"""
Unit tests for the shared configuration and the stdlib-only helpers.

Run from the repository root:

    python tests/test_pipeline_common.py

The first group pins the invariant this project depends on: every target runs one
code path and differs only in its TARGETS entry. The scoring groups pin the
guards using the real numbers measured in the KRAS, EGFR and mTOR runs.
"""
import json
import math
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import pipeline_common as pc

passed = failed = 0


def _try(fn):
    """True when fn() raises (i.e. the expected failure happened)."""
    try:
        fn()
        return False
    except Exception:
        return True


def check(label, cond, detail=""):
    global passed, failed
    if cond:
        passed += 1
        print(f"  PASS  {label}")
    else:
        failed += 1
        print(f"  FAIL  {label}  {detail}")


tmp = Path(tempfile.mkdtemp())

print("\n--- all three targets load, and differ only in config ---")
cfgs = {}
for name in ("kras", "egfr", "mtor"):
    cfgs[name] = pc.get_config(name, base_dir=str(tmp / name), persist=True, quiet=True)
    check(f"{name} config loads", cfgs[name]["target_name"] == name)

shared_keys = [k for k in pc.SHARED if k != "base_dir"]
differing = [k for k in shared_keys
             if len({json.dumps(cfgs[t].get(k), sort_keys=True) for t in cfgs}) > 1]
check("shared settings identical across targets", not differing, f"differ: {differing}")

target_keys = set(pc.TARGETS["kras"]) | set(pc.TARGETS["egfr"]) | set(pc.TARGETS["mtor"])
check("every target defines the same keys",
      all(set(pc.TARGETS[t]) == target_keys for t in pc.TARGETS),
      f"{ {t: sorted(set(target_keys) - set(pc.TARGETS[t])) for t in pc.TARGETS} }")

check("unknown target raises a clear error",
      _try(lambda: pc.get_config("notatarget", base_dir=str(tmp), persist=False, quiet=True)))

print("\n--- config persists and restores (standalone step after kernel restart) ---")
restored = pc.get_config(base_dir=str(tmp / "egfr"), persist=False, quiet=True)
check("restores egfr from pipeline_config.json", restored["target_name"] == "egfr",
      f"got {restored['target_name']}")

print("\n--- export_globals publishes the legacy names ---")
ns = {}
pc.export_globals(cfgs["kras"], ns)
for name, expect in (("TARGET_NAME", "kras"), ("MUTANT_FILTER", "G12D"),
                     ("TARGET_CHEMBL_ID", "CHEMBL2189121"), ("MANUAL_PDB_ID", "7RT1"),
                     ("POTENCY_CUTOFF_NM", 50.0)):
    check(f"{name} == {expect!r}", ns.get(name) == expect, f"got {ns.get(name)!r}")
check("RESULTS_DIR is a Path", isinstance(ns.get("RESULTS_DIR"), Path))

print("\n--- normalize_fixed uses fixed bounds (not batch min/max) ---")
lo, hi = 5.0, 11.0
tight = [8.00, 8.05, 8.10]
scores = pc.normalize_fixed(tight, lo, hi)
check("a tight cluster stays tightly scored",
      max(scores) - min(scores) < 0.02, f"spread {max(scores) - min(scores):.4f}")
check("values are not stretched to 0 and 1",
      min(scores) > 0.4 and max(scores) < 0.6, f"{scores}")
check("clamps below range", pc.normalize_fixed([3.0], lo, hi)[0] == 0.0)
check("clamps above range", pc.normalize_fixed([99.0], lo, hi)[0] == 1.0)
check("reverse maps the most negative affinity to 1.0",
      pc.normalize_fixed([-12.0], -12.0, -4.0, reverse=True)[0] == 1.0)
check("reverse maps the weakest affinity to 0.0",
      pc.normalize_fixed([-4.0], -12.0, -4.0, reverse=True)[0] == 0.0)
check("NaN scores 0 rather than crashing", pc.normalize_fixed([float('nan')], lo, hi)[0] == 0.0)

print("\n--- spread_is_meaningful rejects noise-level spreads ---")
check("EGFR-like spread (0.06) vs RMSE 0.47 -> not rankable",
      not pc.spread_is_meaningful([7.80, 7.83, 7.86], 0.473, factor=2.0))
check("KRAS-like spread (0.60) vs RMSE 0.48 -> not rankable",
      not pc.spread_is_meaningful([7.74, 8.05, 8.34], 0.479, factor=2.0))
check("a real 2-log spread vs RMSE 0.5 -> rankable",
      pc.spread_is_meaningful([6.0, 7.0, 8.0], 0.5, factor=2.0))
check("mTOR docking spread (0.97) vs Vina 2.85 -> not rankable",
      not pc.spread_is_meaningful([-5.86, -5.30, -4.89], 2.85, factor=1.0))
check("single value -> not rankable", not pc.spread_is_meaningful([1.0], 0.1))

print("\n--- json_safe produces parseable JSON ---")
payload = {"a": float("nan"), "b": [1.0, float("inf")], "c": {"d": float("-inf")}, "e": "ok", "f": 3}
safe = pc.json_safe(payload)
text = json.dumps(safe, allow_nan=False)
check("NaN -> null", safe["a"] is None)
check("inf inside a list -> null", safe["b"][1] is None)
check("-inf nested -> null", safe["c"]["d"] is None)
check("strings and ints survive", safe["e"] == "ok" and safe["f"] == 3)
check("strict parser accepts the output", json.loads(text)["a"] is None)

print("\n--- ligand_efficiency ---")
check("LE = -affinity/heavy atoms",
      abs(pc.ligand_efficiency(-8.0, 32) - 0.25) < 1e-9,
      pc.ligand_efficiency(-8.0, 32))
check("a bigger molecule with the same score has lower LE",
      pc.ligand_efficiency(-8.0, 40) < pc.ligand_efficiency(-8.0, 20))
check("zero atoms -> None", pc.ligand_efficiency(-8.0, 0) is None)

print("\n--- fit_box_size ---")
info = {"extent": 13.7, "reference_coords": [(0, 0, 0), (13.7, 0, 0)]}
size = pc.fit_box_size({"extent": 13.7, "reference_coords": []},
                       {"box_size": 22.0, "box_padding": 5.0, "box_size_max": 40.0},
                       verbose=False)
check("grows to extent + 2*padding when that exceeds the minimum",
      abs(size - 23.7) < 1e-9, size)
small = pc.fit_box_size({"extent": 6.0, "reference_coords": []},
                        {"box_size": 22.0, "box_padding": 5.0, "box_size_max": 40.0},
                        verbose=False)
check("keeps the configured minimum for a small ligand", abs(small - 22.0) < 1e-9, small)
capped = pc.fit_box_size({"extent": 60.0, "reference_coords": []},
                         {"box_size": 22.0, "box_padding": 5.0, "box_size_max": 40.0},
                         verbose=False)
check("respects the cap", abs(capped - 40.0) < 1e-9, capped)

print("\n--- centroid / extent ---")
check("centroid of a cube corner pair",
      pc.centroid([(0, 0, 0), (2, 2, 2)]) == (1.0, 1.0, 1.0))
check("extent is the longest axis span",
      abs(pc.group_extent([(0, 0, 0), (4, 1, 1)], (2, 0.5, 0.5)) - 4.0) < 1e-9)
check("empty centroid raises", _try(lambda: pc.centroid([])))

print("\n" + "=" * 60)
print(f"PASSED {passed}   FAILED {failed}")
sys.exit(1 if failed else 0)
