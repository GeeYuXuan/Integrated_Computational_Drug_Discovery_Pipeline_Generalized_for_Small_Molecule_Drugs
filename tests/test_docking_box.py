"""
Regression test for the docking box centre.

This is the test that the original pipeline needed and did not have. It runs
against the three real PDB structures committed in the *_project_backup_* folders
and needs nothing but the Python standard library — no RDKit, no Vina.

Run from the repository root:

    python tests/test_docking_box.py

The old implementation averaged the coordinates of every non-water HETATM record
into one point. That is correct only when a structure contains exactly one hetero
group. Measured here:

    EGFR 5Y9T   1 hetero group                 0.0 A off   (accidentally correct)
    KRAS 7RT1   2 ligand copies + GDP + Mg     7.3 A off   (pocket clipped)
    mTOR 4JSV   ADP x2 + ions, in a dimer     37.7 A off   (box held no pocket)
"""
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import pipeline_common as pc  # noqa: E402

# Structures live in tests/data/ so this regression does not depend on which
# example-output folders happen to be committed. They are unmodified RCSB
# downloads — inputs to the pipeline, not results from it — which is why they
# survive when a stale run is deleted.
STRUCTURES = Path(__file__).resolve().parent / "data"

# (target, structure, true pocket centre = centroid of ONE reference ligand copy)
CASES = [
    ("kras", "7RT1.pdb", (1.58, 4.75, -23.51)),
    ("egfr", "5Y9T.pdb", (-0.46, 49.56, -17.66)),
    ("mtor", "4JSV.pdb", (50.12, -1.80, -45.90)),
]


def old_box_center(pdb_path):
    """Reproduce the original implementation, for comparison."""
    pts = []
    with open(pdb_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("HETATM") and "HOH" not in line:
                pts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            elif line.startswith("END"):
                break
    return tuple(sum(p[i] for p in pts) / len(pts) for i in range(3))


def check_pose_parser():
    """parse_vina_poses() must read exactly what the Python bindings reported.

    Docking now runs through the `vina` executable rather than the `vina` Python
    package, which has no wheel for Python 3.13 (Colab's current interpreter).
    The affinity this parser reads from `REMARK VINA RESULT:` must equal the
    affinity recorded in the matching docking_results.csv.

    Only the KRAS deposit is checked. The parser is target-agnostic, and the
    EGFR and mTOR example outputs were removed because they predated the box,
    cutoff and training-range fixes.
    """
    import csv
    checked = mismatched = 0
    for target, backup in (("kras", "kras_project_backup_16Nov26"),):
        results = ROOT / backup / "results" / f"{target}_docking_results.csv"
        if not results.exists():
            continue
        with open(results, encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                pose = results.parent / f"{row['Candidate_ID']}_docked.pdbqt"
                if not pose.exists():
                    continue
                parsed = pc.parse_vina_poses(pose)
                recorded = float(row["Affinity_kcal_mol"])
                checked += 1
                if not parsed or abs(parsed[0] - recorded) > 1e-6:
                    mismatched += 1
                    print(f"  MISMATCH {target} {row['Candidate_ID']}: "
                          f"csv {recorded}, parsed {parsed[:1]}")
    print(f"\n  parse_vina_poses: {checked} pose files checked, "
          f"{mismatched} mismatched")
    return [] if not mismatched else [f"{mismatched} pose files mis-parsed"]


def main():
    failures = []
    for target, fname, expected in CASES:
        path = STRUCTURES / fname
        if not path.exists():
            print(f"SKIP {target}: tests/data/{fname} not present")
            continue

        print("\n" + "=" * 72)
        print(f"{target.upper()}  ({path.name})")
        print("=" * 72)

        cfg = pc.get_config(target, base_dir=str(ROOT / ".pytest_scratch"),
                            persist=False, quiet=True)
        center, info = pc.resolve_box_center(path, cfg, verbose=True)
        size = pc.fit_box_size(info, cfg, verbose=True)

        d_new = math.dist(center, expected)
        d_old = math.dist(old_box_center(path), expected)
        coords = info["reference_coords"]
        half = size / 2.0
        inside_new = sum(1 for p in coords
                         if all(abs(p[i] - center[i]) <= half for i in range(3)))
        old_c = old_box_center(path)
        inside_old = sum(1 for p in coords
                         if all(abs(p[i] - old_c[i]) <= 10.0 for i in range(3)))

        print(f"\n  old centre is {d_old:6.2f} A from the pocket; "
              f"{inside_old}/{len(coords)} reference atoms inside a 20 A box")
        print(f"  new centre is {d_new:6.2f} A from the pocket; "
              f"{inside_new}/{len(coords)} reference atoms inside a {size:.1f} A box")

        if d_new > 0.05:
            failures.append(f"{target}: box centre {d_new:.2f} A from the pocket")
        if inside_new != len(coords):
            failures.append(f"{target}: box holds only {inside_new}/{len(coords)} "
                            f"reference atoms")
        if info["source"] != "config:ligand_resname":
            failures.append(f"{target}: expected the pinned reference ligand, "
                            f"got {info['source']}")

    print("\n" + "=" * 72)
    print("POSE PARSER")
    print("=" * 72)
    failures += check_pose_parser()

    print("\n" + "=" * 72)
    if failures:
        print("FAILED")
        for f in failures:
            print("  - " + f)
        return 1
    print(f"PASSED — box centred on the pocket for all {len(CASES)} structures, "
          f"and the pose parser agrees with every recorded affinity")
    return 0


def test_box_regression():
    """Entry point for `pytest`; main() returns 0 only if every case passed."""
    assert main() == 0


if __name__ == "__main__":
    sys.exit(main())
