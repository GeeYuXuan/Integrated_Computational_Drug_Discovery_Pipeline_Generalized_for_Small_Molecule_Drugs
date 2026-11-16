"""
Regression tests for the bugs an external peer review of the KRAS report found.

Each test pins a number that used to be wrong, so the specific mistake cannot
come back unnoticed:

  1. The PDB resolution parser returned the REMARK number instead of the
     resolution, so every structure was reported as a 2.0 A structure.
  2. The pocket-residue count used whichever hetero group the file listed first
     (GDP in 7RT1, not the inhibitor) and counted ordered waters as protein
     residues, which turned 23 protein residues into "37 pocket residues".
  3. Step 1 kept only the first assay id per aggregated compound, so every
     "unique assays" figure downstream was an undercount.

Run from the repository root:

    python tests/test_review_fixes.py
"""
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = Path(__file__).resolve().parent / "data"
NOTEBOOK = ROOT / (
    "Integrated_Computational_Drug_Discovery_Pipeline_Generalized_for_Small_Molecule_Drugs.ipynb"
)
BACKUP = ROOT / "kras_project_backup_16Nov26"


def _notebook_sources():
    nb = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return ["".join(c["source"]) for c in nb["cells"]]


def _extract_function(name):
    """Pull one top-level function out of the notebook and compile it.

    The pipeline's logic lives in notebook cells rather than in importable
    modules, so a test either duplicates the code -- and then tests its own
    copy -- or lifts the real source. This lifts the real source.
    """
    for src in _notebook_sources():
        m = re.search(rf"^def {re.escape(name)}\(.*?(?=^\S|\Z)", src,
                      re.M | re.S)
        if m:
            ns = {}
            exec(compile(m.group(0), f"<notebook:{name}>", "exec"), ns)
            return ns[name]
    raise AssertionError(f"{name}() not found in the notebook")


# ---------------------------------------------------------------- resolution
def check_resolution():
    """"REMARK   2 RESOLUTION.    1.27 ANGSTROMS." must read as 1.27.

    The old loop returned the first float-convertible token on the line, which
    is the remark number 2.
    """
    fn = _extract_function("get_pdb_resolution_from_file")
    cases = [("7RT1.pdb", 1.27), ("5Y9T.pdb", 3.25), ("4JSV.pdb", 3.50)]
    failures = []
    for fname, expected in cases:
        path = DATA / fname
        if not path.exists():
            print(f"  SKIP {fname}: not present")
            continue
        got = fn(str(path))
        print(f"  {fname}: parsed {got!r}, expected {expected}")
        if not isinstance(got, float) or abs(got - expected) > 1e-9:
            failures.append(f"{fname}: resolution parsed as {got!r}, "
                            f"expected {expected}")
    return failures


# -------------------------------------------------------------------- pocket
def check_pocket():
    """Waters are not pocket residues, and the reference ligand is 7L8.

    Biopython marks waters with the hetflag "W", not "H_", so a test for
    "H_" lets them through. Around 7L8 in 7RT1 chain A there are 23 protein
    residues and 16 ordered waters within 5.0 A; the published figure of 37
    was 24 protein residues plus 13 waters around GDP, a different site.
    """
    try:
        from Bio import PDB
    except ImportError:
        print("  SKIP: biopython not installed")
        return []
    path = DATA / "7RT1.pdb"
    if not path.exists():
        print("  SKIP: tests/data/7RT1.pdb not present")
        return []

    import warnings
    warnings.filterwarnings("ignore")
    model = PDB.PDBParser(QUIET=True).get_structure("x", str(path))[0]

    excluded = ["HOH", "WAT", "SO4", "PO4", "EDO", "DMS", "MG", "NA", "CL"]
    ligands = [r for r in model.get_residues()
               if r.get_id()[0].startswith("H_")
               and r.get_resname() not in excluded]
    first = ligands[0].get_resname()
    chosen = next((l for l in ligands if l.get_resname() == "7L8"), ligands[0])
    print(f"  hetero groups: {[l.get_resname() for l in ligands]}")
    print(f"  file lists {first} first; the reference ligand is "
          f"{chosen.get_resname()}")

    protein = [a for c in model for r in c for a in r.get_atoms()
               if r.get_id()[0] == " "]
    waters = [a for c in model for r in c for a in r.get_atoms()
              if r.get_id()[0] == "W"]
    ns_p, ns_w = PDB.NeighborSearch(protein), PDB.NeighborSearch(waters)
    near_p, near_w = set(), set()
    for atom in chosen.get_atoms():
        near_p.update(ns_p.search(atom.get_coord(), 5.0, level="R"))
        near_w.update(ns_w.search(atom.get_coord(), 5.0, level="R"))
    print(f"  within 5.0 A of {chosen.get_resname()}: {len(near_p)} protein "
          f"residues, {len(near_w)} ordered waters")

    failures = []
    if first == "7L8":
        failures.append("7RT1 no longer lists GDP first; this test's premise "
                        "has changed and the selection logic needs rechecking")
    if len(near_p) != 23:
        failures.append(f"pocket: {len(near_p)} protein residues, expected 23")
    if len(near_w) != 16:
        failures.append(f"pocket: {len(near_w)} waters, expected 16")
    return failures


# --------------------------------------------------------------------- assays
def check_assay_counts():
    """The measurement table must be the source of any assay count.

    Aggregating with assay_chembl_id=('assay_chembl_id', 'first') keeps one id
    per compound, so counting distinct ids in the aggregate table undercounts
    the assays behind the data: 16 instead of 25 for the curated set.
    """
    meas = BACKUP / "results" / "kras_measurements.csv"
    agg = BACKUP / "results" / "kras_bioactivity_full.csv"
    if not meas.exists() or not agg.exists():
        print("  SKIP: deposited result tables not present")
        return []

    def read(p):
        with open(p, encoding="utf-8") as fh:
            return list(csv.DictReader(fh))

    rows_m, rows_a = read(meas), read(agg)
    n_true = len({r["assay_chembl_id"] for r in rows_m})
    n_agg = len({r["assay_chembl_id"] for r in rows_a})
    counts = {}
    for r in rows_m:
        counts[r["assay_chembl_id"]] = counts.get(r["assay_chembl_id"], 0) + 1
    dominant = max(counts.values())
    print(f"  {len(rows_m)} measurements over {n_true} assays; "
          f"aggregating by first id shows only {n_agg}")
    print(f"  largest assay contributes {dominant}/{len(rows_m)} measurements "
          f"({100 * dominant / len(rows_m):.1f}%)")

    failures = []
    if n_true != 25:
        failures.append(f"measurement table has {n_true} assays, expected 25")
    if n_agg >= n_true:
        failures.append("the aggregate table no longer undercounts assays; "
                        "if that column now holds every id, update this test")
    if dominant != 675:
        failures.append(f"dominant assay has {dominant} measurements, "
                        f"expected 675")
    return failures


CHECKS = [
    ("PDB RESOLUTION PARSER", check_resolution),
    ("POCKET RESIDUES AND WATERS", check_pocket),
    ("ASSAY COUNTS", check_assay_counts),
]


def main():
    failures = []
    for title, fn in CHECKS:
        print("\n" + "=" * 72)
        print(title)
        print("=" * 72)
        failures += fn()

    print("\n" + "=" * 72)
    if failures:
        print("FAILED")
        for f in failures:
            print("  - " + f)
        return 1
    print("PASSED — every figure the peer review corrected still holds")
    return 0


def test_review_fixes():
    """Entry point for `pytest`."""
    assert main() == 0


if __name__ == "__main__":
    sys.exit(main())
