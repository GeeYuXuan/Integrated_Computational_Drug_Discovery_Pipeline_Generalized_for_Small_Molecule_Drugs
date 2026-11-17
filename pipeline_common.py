"""
Shared configuration and helpers for the Integrated Computational Drug Discovery Pipeline.

DESIGN RULE
-----------
Every target runs **exactly the same code**. The only thing that differs between a
KRAS run, an EGFR run and an mTOR run is one entry in ``TARGETS`` below, selected by
``ACTIVE_TARGET``. Nothing target-specific may be hard-coded anywhere else in the
notebook -- if you find yourself editing a notebook cell to switch targets, that is a
bug in this file, not a reason to edit the cell.

The helpers in this module are deliberately split into two groups:

* ``PDB / geometry`` and ``scoring`` helpers depend only on the Python standard
  library, so they can be unit-tested without RDKit, pandas or Vina installed.
* ``chemistry`` helpers import RDKit lazily, so importing this module never fails
  just because RDKit is missing.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
from datetime import date
from pathlib import Path

__all__ = [
    "ACTIVE_TARGET", "TARGETS", "SHARED", "SOLVENT_AND_IONS", "COFACTORS",
    "get_config", "export_globals", "describe_config",
    "parse_hetero_groups", "centroid", "group_extent", "resolve_box_center",
    "fit_box_size", "write_receptor",
    "find_vina_binary", "parse_vina_poses", "dock_with_vina_cli",
    "normalize_fixed", "spread_is_meaningful", "ligand_efficiency",
    "summarise_scores", "estimate_is_stable",
    "inchikey", "is_peptide_like", "morgan_fp", "max_tanimoto_to_training",
    "json_safe",
]

# =============================================================================
# 1. THE ONLY PLACE A TARGET IS DEFINED
# =============================================================================

#: Which entry of ``TARGETS`` to run. This is the single switch for the whole pipeline.
ACTIVE_TARGET = "kras"

TARGETS = {
    # -------------------------------------------------------------------------
    "kras": {
        "target_name":    "kras",          # used for every output filename
        "search_term":    "KRAS",          # ChEMBL / UniProt / PubMed search term
        "chembl_id":      "CHEMBL2189121",  # GTPase KRas, UniProt P01116
        "mutant_filter":  "G12D",          # None for wild type
        # --- structure ---
        "manual_pdb_id":  "7RT1",          # None => automatic RCSB search
        "ligand_resname": "7L8",           # reference ligand that defines the pocket
        "ligand_chain":   "A",             # pin one copy: 7RT1 holds two (res 203, 204)
        "active_site_residue": 12,         # G12 -- fallback only, never silently used
        "box_center":     None,            # explicit [x, y, z] overrides everything
    },
    # -------------------------------------------------------------------------
    # BRD4 stands in for MYC, which cannot be docked: no drug-like ligand has
    # ever been co-crystallised with the MYC protein, so there is nothing to
    # centre a box on and resolve_box_center refuses rather than guess. BET
    # inhibition is the best-established route to suppressing MYC -- JQ1's
    # original result was MYC suppression -- and 3MXF is the canonical BRD4
    # bromodomain with JQ1 bound. Deposited run: brd4_project_backup_17Nov26/.
    "brd4": {
        "target_name":    "brd4",
        "search_term":    "BRD4",
        "chembl_id":      "CHEMBL1163125",  # Bromodomain-containing protein 4
        "mutant_filter":  None,            # wild type
        # --- structure ---
        "manual_pdb_id":  "3MXF",          # BRD4 bromodomain 1, 1.6 A
        "ligand_resname": "JQ1",           # the co-crystal BET inhibitor
        "ligand_chain":   "A",             # one protein chain, one JQ1 copy
        "active_site_residue": 140,        # Asn140, the acetyl-lysine anchor
        "box_center":     None,
    },
    # -------------------------------------------------------------------------
    "egfr": {
        "target_name":    "egfr",
        "search_term":    "EGFR",
        "chembl_id":      "CHEMBL203",
        "mutant_filter":  None,
        "manual_pdb_id":  "5Y9T",
        "ligand_resname": "8RC",
        "ligand_chain":   "A",
        "active_site_residue": 790,        # T790 gatekeeper
        "box_center":     None,
    },
    # -------------------------------------------------------------------------
    "mtor": {
        "target_name":    "mtor",
        "search_term":    "MTOR",
        "chembl_id":      "CHEMBL2842",
        "mutant_filter":  None,
        "manual_pdb_id":  "4JSV",
        # 4JSV carries no drug-like inhibitor: it is an mTOR-mLST8 dimer with ADP in
        # each of two chains. ATP-competitive inhibitors belong in the ADP site, so we
        # pin ONE chain's ADP. Averaging across both chains puts the box ~37 A away
        # from either site, in empty space -- which is what the old code did.
        "ligand_resname": "ADP",
        "ligand_chain":   "A",
        "active_site_residue": None,
        "box_center":     None,
    },
}

# =============================================================================
# 2. SETTINGS SHARED BY EVERY TARGET
# =============================================================================

SHARED = {
    # --- environment ---
    "base_dir": "/content",             # Colab; set PIPELINE_BASE_DIR to override
    "random_seed": 42,
    # NCBI's E-utilities policy asks automated clients to identify themselves with a
    # working contact address, so they can email you rather than block your IP. It is
    # not a credential and needs no registration. Override with the PIPELINE_ENTREZ_EMAIL
    # environment variable if you are running this as someone else.
    "entrez_email": os.environ.get("PIPELINE_ENTREZ_EMAIL", "yuxuan.gee@u.nus.edu"),

    # --- Step 1: data acquisition ---------------------------------------------
    # Compounds kept as "potent hits" for scaffold analysis and seed selection.
    "potency_cutoff_nm": 50.0,
    # Compounds kept for TRAINING the ML model. Training only on sub-50 nM compounds
    # truncates the label range and the model cannot then rank anything -- so the
    # training set deliberately spans the full measured range.
    "training_max_ic50_nm": 1.0e6,
    # How to collapse repeat measurements of the same compound. "median" is honest;
    # "max" (the old behaviour) biases every compound towards its best-ever number.
    "duplicate_policy": "median",
    # Drop peptides/macrocycles: this is a small-molecule pipeline.
    "exclude_peptides": True,
    "peptide_max_amide_bonds": 4,
    "peptide_max_ring_size": 12,

    # --- Step 2: scaffolds -----------------------------------------------------
    "min_scaffold_count": 5,            # a scaffold needs this many compounds to count

    # --- Step 4: generation ----------------------------------------------------
    "generated_mw_range": [300.0, 600.0],
    "generated_max_logp": 5.5,
    "generated_min_qed": 0.5,
    "n_seeds_proven": 10,
    "n_seeds_exploration": 10,
    "n_top_candidates": 20,

    # --- Step 5: docking -------------------------------------------------------
    "box_size": 22.0,                   # minimum box edge, grown to fit the ligand
    "box_padding": 5.0,                 # A of clearance per side around the reference
    "box_size_max": 40.0,               # refuse to grow beyond this
    # Vina's own default (its floor is 1, not 8). Doubling it to 16 doubles the
    # docking time and buys nothing you can act on here: the 20 candidates land
    # within ~0.5 kcal/mol of each other, far inside Vina's ~2.85 kcal/mol scoring
    # error, so the ordering is noise at either setting. The shorter cell also
    # matters operationally -- a 60-minute docking cell was killed outright by a
    # Colab idle disconnect, losing the whole run.
    "exhaustiveness": 8,
    "num_modes": 9,                     # >1 so pose convergence can be inspected
    # Docking runs through the `vina` executable, not the Python bindings, which
    # have no wheel for Python 3.13 (what Colab now runs). None = find it on PATH.
    "vina_binary": None,
    "vina_cpu": None,                   # None = let Vina use every core
    # Published RMSE of the Vina scoring function. Differences smaller than this
    # between two compounds are not real differences.
    "vina_score_error_kcal": 2.85,
    "n_conformers": 50,
    "min_ligand_atoms": 10,             # a hetero group smaller than this is not a ligand
    "keep_cofactors_in_receptor": True, # keep GDP/Mg etc.: removing them opens fake cavities
    "receptor_chains": None,            # None = keep every chain

    # --- Step 6: machine learning ----------------------------------------------
    "fp_radius": 2,
    "fp_bits": 2048,
    "test_size": 0.2,
    "n_estimators": 500,
    "cv_folds": 5,
    # Scaffold splits are repeated, not drawn once. A single draw gave R2 = 0.794
    # against a random-split 0.614 on the first live KRAS run -- the wrong way
    # round, and a sign the split choice was doing the work.
    "scaffold_split_repeats": 10,
    # A candidate less similar than this to anything in the training set is outside
    # the model's applicability domain and its prediction is not trustworthy.
    "applicability_min_tanimoto": 0.40,

    # --- Step 7: consensus scoring ---------------------------------------------
    # FIXED bounds, so a score means the same thing in every run. Min-max scaling
    # within a 20-compound batch stretches pure noise across the full 0-1 range.
    "score_bounds": {
        "pic50":    [5.0, 11.0],        # pIC50 units
        "affinity": [-12.0, -4.0],      # kcal/mol (more negative = better)
    },
    "weights": {
        "Potency": 0.30,
        "Docking": 0.30,
        "QED": 0.20,
        "PAINS": 0.10,
        "DrugLike": 0.10,
    },
    "n_final_candidates": 5,
}

#: Residues that are solvent, cryoprotectant, buffer or ions -- never a reference ligand.
SOLVENT_AND_IONS = {
    "HOH", "WAT", "DOD", "SO4", "PO4", "EDO", "GOL", "PGE", "PEG", "DMS", "MES",
    "EPE", "TRS", "ACT", "ACY", "FMT", "CIT", "NO3", "IOD", "BME", "MPD", "IPA",
    "MG", "MN", "ZN", "NA", "CL", "K", "CA", "FE", "CU", "NI", "CO", "CD", "BR",
    "F", "MGF", "CO3", "AZI", "SCN", "OXL", "TAR", "SIN", "GLC", "NAG", "BMA",
}

#: Biological cofactors. Usable as a pocket reference, but only with a loud warning,
#: because a cofactor site is not necessarily where an inhibitor binds.
COFACTORS = {
    "GDP", "GTP", "GNP", "GSP", "GCP", "ADP", "ATP", "ANP", "AGS", "ACP", "ADN",
    "NAD", "NAI", "NAP", "NDP", "FAD", "FMN", "SAM", "SAH", "COA", "HEM", "PLP",
    "TPP", "UDP", "UTP", "CDP", "CTP",
}


# =============================================================================
# 3. CONFIG ASSEMBLY  (one dict, persisted so each step can run standalone)
# =============================================================================

def _config_path(base_dir) -> Path:
    return Path(base_dir) / "pipeline_config.json"


def get_config(target: str | None = None, base_dir: str | None = None,
               persist: bool = True, quiet: bool = False) -> dict:
    """Build (or restore) the one configuration dict for this run.

    Resolution order for the target:

    1. the ``target`` argument,
    2. the ``PIPELINE_TARGET`` environment variable,
    3. a ``pipeline_config.json`` written by an earlier step in this session,
    4. :data:`ACTIVE_TARGET`.

    Step 3 is what lets you run Step 5 on its own after a kernel restart without
    editing anything.
    """
    base_dir = base_dir or os.environ.get("PIPELINE_BASE_DIR") or SHARED["base_dir"]
    base = Path(base_dir)

    name = target or os.environ.get("PIPELINE_TARGET")
    if name is None:
        saved = _config_path(base)
        if saved.exists():
            try:
                name = json.loads(saved.read_text(encoding="utf-8"))["target_name"]
            except Exception:
                name = None
    name = (name or ACTIVE_TARGET).lower()

    if name not in TARGETS:
        raise KeyError(
            f"Unknown target {name!r}. Known targets: {sorted(TARGETS)}. "
            f"Add a new entry to TARGETS in pipeline_common.py -- do not edit the "
            f"notebook cells."
        )

    cfg = dict(SHARED)
    cfg.update(TARGETS[name])
    cfg["base_dir"] = str(base)

    # Directory layout (identical for every target).
    cfg["results_dir"] = str(base / "results")
    cfg["figures_dir"] = str(base / "figures")
    cfg["data_dir"] = str(base / "data")
    cfg["workflow_data_dir"] = str(base / "workflow" / "data")
    cfg["ligand_prep_dir"] = str(base / "ligands_prep")
    for key in ("results_dir", "figures_dir", "data_dir",
                "workflow_data_dir", "ligand_prep_dir"):
        Path(cfg[key]).mkdir(parents=True, exist_ok=True)

    if persist:
        _config_path(base).write_text(json.dumps(cfg, indent=2), encoding="utf-8")

    if not quiet:
        print(describe_config(cfg))
    return cfg


def describe_config(cfg: dict) -> str:
    mut = cfg.get("mutant_filter") or "wild type"
    return (
        "=" * 78 + "\n"
        f"  TARGET      : {cfg['target_name'].upper()}  ({cfg['search_term']}, {mut})\n"
        f"  ChEMBL ID   : {cfg['chembl_id']}\n"
        f"  Structure   : {cfg['manual_pdb_id'] or 'automatic RCSB search'}"
        f"   reference ligand: {cfg.get('ligand_resname') or 'auto-detect'}\n"
        f"  Hit cutoff  : IC50 < {cfg['potency_cutoff_nm']} nM"
        f"   |  ML trains on the full measured range\n"
        f"  Output      : {cfg['results_dir']}\n"
        + "=" * 78
    )


#: Legacy global names the notebook cells already use, mapped to config keys. Keeping
#: these means a cell's config block can simply be deleted rather than rewritten.
_GLOBAL_ALIASES = {
    "TARGET_NAME": "target_name",
    "TARGET_SEARCH_TERM": "search_term",
    "PROTEIN_SEARCH_TERM": "search_term",
    "TARGET_CHEMBL_ID": "chembl_id",
    "MUTANT_FILTER": "mutant_filter",
    "POTENCY_CUTOFF_NM": "potency_cutoff_nm",
    "MANUAL_PDB_ID": "manual_pdb_id",
    "TARGET_RESIDUE_ID": "active_site_residue",
    "BOX_SIZE": "box_size",
    "EXHAUSTIVENESS": "exhaustiveness",
    "NUM_MODES": "num_modes",
    "RANDOM_SEED": "random_seed",
}


def export_globals(cfg: dict, namespace: dict) -> dict:
    """Publish ``cfg`` into a notebook namespace under its historical global names.

    Also sets ``BASE_DIR``, ``RESULTS_DIR``, ``FIGURES_DIR``, ``DATA_DIR``,
    ``WORKFLOW_DATA_DIR`` and ``LIGAND_PREP_DIR`` as :class:`~pathlib.Path` objects,
    and seeds every random number generator that is importable.
    """
    for gname, key in _GLOBAL_ALIASES.items():
        if key in cfg:
            namespace[gname] = cfg[key]

    namespace["CFG"] = cfg
    namespace["BASE_DIR"] = Path(cfg["base_dir"])
    namespace["RESULTS_DIR"] = Path(cfg["results_dir"])
    namespace["FIGURES_DIR"] = Path(cfg["figures_dir"])
    namespace["DATA_DIR"] = Path(cfg["data_dir"])
    namespace["WORKFLOW_DATA_DIR"] = Path(cfg["workflow_data_dir"])
    namespace["LIGAND_PREP_DIR"] = Path(cfg["ligand_prep_dir"])

    seed = cfg["random_seed"]
    import random as _random
    _random.seed(seed)
    try:
        import numpy as _np
        _np.random.seed(seed)
    except ImportError:
        pass
    return cfg


def json_safe(obj):
    """Recursively replace NaN/Infinity with ``None`` so output is *valid* JSON.

    ``json.dump`` emits bare ``NaN``, which Python reads back but ``JSON.parse``,
    ``jq`` and most other parsers reject.
    """
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_safe(v) for v in obj]
    if isinstance(obj, float):
        return None if (math.isnan(obj) or math.isinf(obj)) else obj
    # numpy scalars, pandas NA, etc.
    for attr in ("item",):
        if hasattr(obj, attr) and not isinstance(obj, (str, bytes)):
            try:
                return json_safe(obj.item())
            except Exception:
                break
    if obj is not None and not isinstance(obj, (str, bytes, int, bool)):
        try:
            if obj != obj:  # NaN-like
                return None
        except Exception:
            pass
    return obj


def today() -> str:
    return date.today().isoformat()


# =============================================================================
# 4. PDB GEOMETRY  -- stdlib only, so it is unit-testable
# =============================================================================

def parse_hetero_groups(pdb_path) -> dict:
    """Return ``{(resname, chain, resseq): [(x, y, z), ...]}`` for heavy HETATM atoms.

    Only the first model is read, only the blank/``A`` alternate location is kept, and
    water and hydrogens are skipped. Grouping by *residue* is the whole point: the
    previous implementation averaged every hetero atom in the file into a single
    point, which silently produced a box centre in empty space whenever a structure
    held more than one hetero group (a dimer, a cofactor, an ion, a second ligand copy).
    """
    groups: dict = {}
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            record = line[:6]
            if record.startswith("ENDMDL"):
                break
            if not record.startswith("HETATM"):
                continue
            resname = line[17:20].strip().upper()
            if resname in ("HOH", "WAT", "DOD"):
                continue
            if line[16] not in (" ", "A"):
                continue
            element = line[76:78].strip().upper()
            if element == "H" or (not element and line[12:16].strip()[:1] == "H"):
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            groups.setdefault((resname, line[21], line[22:27].strip()), []).append(xyz)
    return groups


def parse_residue_atoms(pdb_path, resseq: int, chain: str | None = None) -> list:
    """Heavy ``ATOM`` coordinates of one protein residue (fallback pocket reference)."""
    coords = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line[:6].startswith("ENDMDL"):
                break
            if not line[:6].startswith("ATOM"):
                continue
            if line[16] not in (" ", "A"):
                continue
            if chain is not None and line[21] != chain:
                continue
            try:
                if int(line[22:26]) != int(resseq):
                    continue
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    return coords


def centroid(coords) -> tuple:
    n = len(coords)
    if not n:
        raise ValueError("centroid() of an empty coordinate list")
    return tuple(sum(c[i] for c in coords) / n for i in range(3))


def group_extent(coords, center=None) -> float:
    """Longest edge of the axis-aligned bounding box needed to hold ``coords``."""
    c = center or centroid(coords)
    return 2.0 * max(max(abs(p[i] - c[i]) for p in coords) for i in range(3))


def resolve_box_center(pdb_path, cfg: dict, verbose: bool = True) -> tuple:
    """Decide where the docking grid belongs. Returns ``(center, info)``.

    Resolution order, strictest first -- and it **raises** rather than quietly
    falling back to the geometric centre of the protein, which is never a pocket:

    1. ``cfg['box_center']``                 -- explicit coordinates.
    2. ``cfg['ligand_resname']`` (+ chain)   -- the pinned reference ligand.
    3. auto: the single largest drug-like hetero group.
    4. auto: the single largest cofactor group, with a warning.
    5. ``cfg['active_site_residue']``        -- a named residue, with a warning.

    At no point are two hetero groups averaged together.
    """
    info: dict = {"pdb": str(pdb_path)}

    if cfg.get("box_center"):
        center = tuple(float(v) for v in cfg["box_center"])
        info.update(source="config:box_center", reference="explicit", n_atoms=0,
                    extent=0.0, reference_coords=[])
        if verbose:
            print(f"  Box centre from config: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
        return center, info

    groups = parse_hetero_groups(pdb_path)
    if verbose:
        print(f"  Hetero groups in {Path(pdb_path).name}: {len(groups)}")
        for (rn, ch, rs), atoms in sorted(groups.items(), key=lambda kv: -len(kv[1]))[:8]:
            tag = ("solvent/ion" if rn in SOLVENT_AND_IONS
                   else "cofactor" if rn in COFACTORS else "drug-like")
            print(f"    {rn:>4} chain {ch or '-'} res {rs:>5}  {len(atoms):>3} atoms  [{tag}]")

    def _pick(pool, source, warn=None):
        # Deterministic: most atoms first, then chain, then residue number.
        key, coords = sorted(pool.items(), key=lambda kv: (-len(kv[1]), kv[0][1], kv[0][2]))[0]
        center = centroid(coords)
        copies = [k for k in pool if k[0] == key[0]]
        if warn and verbose:
            print(f"  !! {warn}")
        if len(copies) > 1 and verbose:
            print(f"  !! {len(copies)} copies of {key[0]} present; using chain "
                  f"{key[1]} res {key[2]} only (never the average of all copies)")
        if verbose:
            print(f"  Pocket reference: {key[0]} chain {key[1]} res {key[2]} "
                  f"({len(coords)} atoms)")
            print(f"  Box centre: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
        info.update(source=source, reference=f"{key[0]}/{key[1]}/{key[2]}",
                    n_atoms=len(coords), extent=group_extent(coords, center),
                    reference_coords=coords, n_copies=len(copies))
        return center, info

    want = cfg.get("ligand_resname")
    if want:
        chain = cfg.get("ligand_chain")
        pool = {k: v for k, v in groups.items()
                if k[0] == want.upper() and (chain is None or k[1] == chain)}
        if not pool:
            present = sorted({k[0] for k in groups})
            raise LookupError(
                f"Reference ligand {want!r}"
                + (f" in chain {chain!r}" if chain else "")
                + f" is not present in {Path(pdb_path).name}. Hetero residues found: "
                f"{present}. Fix 'ligand_resname'/'ligand_chain' for target "
                f"{cfg['target_name']!r} in pipeline_common.py."
            )
        return _pick(pool, "config:ligand_resname")

    min_atoms = cfg.get("min_ligand_atoms", 10)
    drug_like = {k: v for k, v in groups.items()
                 if k[0] not in SOLVENT_AND_IONS and k[0] not in COFACTORS
                 and len(v) >= min_atoms}
    if drug_like:
        return _pick(drug_like, "auto:drug-like-ligand")

    cofactor = {k: v for k, v in groups.items()
                if k[0] in COFACTORS and len(v) >= min_atoms}
    if cofactor:
        return _pick(cofactor, "auto:cofactor", warn=(
            "No drug-like ligand in this structure; centring on a COFACTOR site. "
            "That is right for an ATP/GTP-competitive inhibitor and wrong for an "
            "allosteric one -- check this before trusting the scores."))

    residue = cfg.get("active_site_residue")
    if residue is not None:
        coords = parse_residue_atoms(pdb_path, residue, cfg.get("ligand_chain"))
        if coords:
            center = centroid(coords)
            if verbose:
                print(f"  !! No ligand found; centring on residue {residue} "
                      f"({len(coords)} atoms). A residue centroid is a weaker "
                      f"reference than a bound ligand -- treat results as indicative.")
                print(f"  Box centre: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
            info.update(source="config:active_site_residue", reference=f"residue {residue}",
                        n_atoms=len(coords), extent=group_extent(coords, center),
                        reference_coords=coords)
            return center, info

    raise LookupError(
        f"Could not locate a binding site in {Path(pdb_path).name}. No drug-like "
        f"ligand, no cofactor, and no usable 'active_site_residue' for target "
        f"{cfg['target_name']!r}. Set 'ligand_resname' or 'box_center' in "
        f"pipeline_common.py rather than letting the pipeline guess -- the geometric "
        f"centre of a protein is not a pocket."
    )


def fit_box_size(info: dict, cfg: dict, verbose: bool = True) -> float:
    """Return a box edge that actually encloses the pocket reference, plus padding.

    The old code used a fixed 20 A edge regardless of the reference ligand's size, so
    a large ligand was clipped by the search space without any warning.
    """
    requested = float(cfg.get("box_size", 24.0))
    padding = float(cfg.get("box_padding", 8.0))
    cap = float(cfg.get("box_size_max", 40.0))
    extent = float(info.get("extent") or 0.0)

    needed = extent + 2.0 * padding
    size = max(requested, needed)
    if size > cap:
        if verbose:
            print(f"  !! Pocket reference needs a {size:.1f} A box but the cap is "
                  f"{cap:.1f} A; using the cap. Poses may be clipped.")
        size = cap
    if verbose:
        if size > requested:
            print(f"  Box grown {requested:.1f} -> {size:.1f} A to fit a "
                  f"{extent:.1f} A reference ligand with {padding:.1f} A padding")
        else:
            print(f"  Box size {size:.1f} A (reference extent {extent:.1f} A, "
                  f"clearance {(size - extent) / 2:.1f} A per side)")
    coords = info.get("reference_coords") or []
    if coords:
        c = centroid(coords)
        worst = max(max(abs(p[i] - c[i]) for p in coords) for i in range(3))
        if worst > size / 2.0:
            raise ValueError(
                f"Box of {size:.1f} A does not enclose the pocket reference "
                f"(needs {2 * worst:.1f} A). Raise 'box_size_max'."
            )
    return size


def write_receptor(pdb_path, out_path, cfg: dict, exclude_reference: dict | None = None,
                   verbose: bool = True) -> Path:
    """Write a docking receptor: protein, plus structurally integral cofactors.

    Stripping every ``HETATM`` -- the old behaviour -- removes GDP/Mg from KRAS and
    ADP/Mg from mTOR, opening a cavity that does not exist in the real protein and
    changing the electrostatics of the pocket next door. The reference ligand itself
    *is* removed, since that is what we are docking against.
    """
    keep_cofactors = bool(cfg.get("keep_cofactors_in_receptor", True))
    chains = cfg.get("receptor_chains")
    skip = set()
    if exclude_reference and exclude_reference.get("reference"):
        parts = str(exclude_reference["reference"]).split("/")
        if len(parts) == 3:
            skip.add((parts[0], parts[1], parts[2]))

    kept_atoms = kept_cofactors = 0
    out_lines = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            record = line[:6]
            if record.startswith("ENDMDL"):
                break
            if line[16] not in (" ", "A"):
                continue
            if chains and line[21] not in chains:
                continue
            if record.startswith("ATOM"):
                out_lines.append(line)
                kept_atoms += 1
            elif record.startswith("HETATM") and keep_cofactors:
                resname = line[17:20].strip().upper()
                if resname in ("HOH", "WAT", "DOD"):
                    continue
                if (resname, line[21], line[22:27].strip()) in skip:
                    continue
                if resname in COFACTORS or resname in ("MG", "MN", "ZN", "CA", "FE"):
                    out_lines.append(line)
                    kept_cofactors += 1
    out_lines.append("END\n")

    out_path = Path(out_path)
    out_path.write_text("".join(out_lines), encoding="utf-8")
    if verbose:
        print(f"  Receptor written: {out_path.name} "
              f"({kept_atoms} protein atoms, {kept_cofactors} cofactor/metal atoms kept)")
        if skip:
            print(f"  Reference ligand {sorted(skip)[0][0]} removed from the receptor")
    return out_path


# =============================================================================
# 4b. RUNNING AUTODOCK VINA  -- through the command line, not the Python bindings
# =============================================================================
#
# The `vina` Python package publishes no wheel for Python 3.13, which is what Colab
# now runs, and building it from source there fails:
#
#     Building wheel for vina (pyproject.toml) did not run successfully.
#     ERROR: Failed building wheel for vina
#
# Worse, putting it in the same `pip install` as everything else made its failure
# abort the whole command, so RDKit and Biopython were not installed either.
#
# The apt package `autodock-vina` provides /usr/bin/vina (v1.2.5) -- the same engine
# the bindings wrap. Driving it through the command line removes the Python-version
# dependency entirely and gives the same numbers, which the pose files prove: the
# affinity written as `REMARK VINA RESULT:` is exactly what the bindings return.

def find_vina_binary(explicit: str | None = None):
    """Locate the AutoDock Vina executable, or return None."""
    candidates = [explicit] if explicit else []
    candidates += ["vina", "autodock_vina", "vina_1.2.5_linux_x86_64",
                   "/usr/bin/vina", "/usr/local/bin/vina"]
    for name in candidates:
        if not name:
            continue
        path = shutil.which(name) or (name if Path(name).is_file() else None)
        if path:
            return str(path)
    return None


def parse_vina_poses(pose_path) -> list:
    """Affinities (kcal/mol) from a Vina output PDBQT, best pose first.

    Vina writes one ``REMARK VINA RESULT:`` line per pose:

        REMARK VINA RESULT:    -8.186      0.000      0.000
                               ^ affinity  ^ rmsd l.b ^ rmsd u.b
    """
    energies = []
    try:
        with open(pose_path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            energies.append(float(parts[3]))
                        except ValueError:
                            pass
    except OSError:
        return []
    return energies


def dock_with_vina_cli(receptor_pdbqt, ligand_pdbqt, center, box_size, out_pose,
                       cfg: dict | None = None, vina_bin: str | None = None,
                       timeout: int = 3600):
    """Dock one ligand. Returns ``(best_affinity, all_affinities, error)``.

    ``error`` is None on success, otherwise a short string explaining the failure.
    """
    cfg = cfg or SHARED
    vina_bin = vina_bin or find_vina_binary(cfg.get("vina_binary"))
    if not vina_bin:
        return None, None, (
            "AutoDock Vina executable not found. Install it with "
            "`apt-get install -y autodock-vina` (Colab/Debian) or "
            "`brew install autodock-vina` (macOS), or set 'vina_binary' in "
            "pipeline_common.py to its full path.")

    cmd = [
        vina_bin,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{box_size:.3f}",
        "--size_y", f"{box_size:.3f}",
        "--size_z", f"{box_size:.3f}",
        "--exhaustiveness", str(int(cfg.get("exhaustiveness", 8))),
        "--num_modes", str(int(cfg.get("num_modes", 9))),
        "--seed", str(int(cfg.get("random_seed", 42))),
        "--out", str(out_pose),
    ]
    if cfg.get("vina_cpu"):
        cmd += ["--cpu", str(int(cfg["vina_cpu"]))]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, None, f"timed out after {timeout}s"
    except OSError as exc:
        return None, None, f"could not run {vina_bin}: {exc}"

    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()
        return None, None, (tail[-1][:200] if tail else f"exit code {proc.returncode}")

    energies = parse_vina_poses(out_pose)
    if not energies:
        return None, None, "docking produced no poses"
    return min(energies), energies, None


# =============================================================================
# 5. SCORING  -- stdlib only
# =============================================================================

def normalize_fixed(values, lower: float, upper: float, reverse: bool = False):
    """Scale to 0-1 against **fixed** bounds, so scores compare across runs.

    Min-max scaling inside one batch of 20 candidates always awards 1.0 to the best
    and 0.0 to the worst however small the real difference is, which turns prediction
    noise into an apparently decisive ranking.
    """
    if upper == lower:
        raise ValueError("normalize_fixed() needs lower != upper")
    span = upper - lower

    def one(v):
        if v is None:
            return 0.0
        try:
            if v != v:  # NaN
                return 0.0
        except Exception:
            return 0.0
        x = (float(v) - lower) / span
        x = min(1.0, max(0.0, x))
        return 1.0 - x if reverse else x

    if hasattr(values, "apply"):          # pandas Series
        return values.apply(one)
    if isinstance(values, (list, tuple)):
        return [one(v) for v in values]
    return one(values)


def spread_is_meaningful(values, error, factor: float = 2.0) -> bool:
    """True when a set of values varies by enough to be worth ranking.

    ``error`` is the measurement/model error in the same units (e.g. the model's test
    RMSE, or Vina's ~2.85 kcal/mol). Ranking a spread smaller than the error is
    ranking noise.
    """
    clean = []
    for v in values:
        try:
            if v is None or v != v:
                continue
            clean.append(float(v))
        except Exception:
            continue
    if len(clean) < 2 or not error:
        return False
    return (max(clean) - min(clean)) >= factor * float(error)


def summarise_scores(values) -> dict | None:
    """mean / sd / min / max over repeated train-test splits.

    One split is a single draw from a wide distribution, and quoting it invites
    quoting a lucky number. The first live KRAS run made the point: a single
    scaffold split returned R2 = 0.794 while the random split returned 0.614 --
    backwards, since scaffold splits are normally the harder test. Repeating the
    split and reporting the spread says how much of that was the draw.
    """
    clean = []
    for v in values:
        try:
            if v is None or v != v:          # None or NaN
                continue
            clean.append(float(v))
        except (TypeError, ValueError):
            continue
    if not clean:
        return None
    n = len(clean)
    mean = sum(clean) / n
    sd = (sum((v - mean) ** 2 for v in clean) / (n - 1)) ** 0.5 if n > 1 else 0.0
    return {
        "mean": mean,
        "sd": sd,
        "min": min(clean),
        "max": max(clean),
        "n": n,
        "values": [round(v, 4) for v in clean],
    }


def estimate_is_stable(summary: dict | None, max_sd: float = 0.10,
                       max_range: float = 0.25) -> bool:
    """True when repeated splits agree closely enough to quote one number.

    Defaults are chosen for R2: an sd above 0.10, or a best-to-worst range above
    0.25, means the split choice moves the answer more than most people would
    accept in a reported figure.
    """
    if not summary or summary.get("n", 0) < 2:
        return False
    return (summary["sd"] <= max_sd
            and (summary["max"] - summary["min"]) <= max_range)


def ligand_efficiency(affinity: float, heavy_atoms: int):
    """-affinity / heavy atom count. Raw Vina scores reward size, not fit."""
    if not heavy_atoms or affinity is None:
        return None
    try:
        return -float(affinity) / float(heavy_atoms)
    except Exception:
        return None


# =============================================================================
# 6. CHEMISTRY  -- RDKit imported lazily
# =============================================================================

def _rdkit():
    from rdkit import Chem
    return Chem


def inchikey(smiles: str):
    """Structure key for identity comparison. SMILES strings compare unreliably."""
    Chem = _rdkit()
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToInchiKey(mol) if mol is not None else None
    except Exception:
        return None


def is_peptide_like(mol, cfg: dict | None = None) -> bool:
    """True for peptides and large macrocycles, which are not small-molecule drugs.

    The KRAS dataset contains disulfide-bridged cyclic peptides; left in, they become
    "privileged scaffolds" and seeds for a small-molecule generator.
    """
    cfg = cfg or SHARED
    Chem = _rdkit()
    if mol is None:
        return False
    try:
        amide = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
        n_amide = len(mol.GetSubstructMatches(amide)) if amide is not None else 0
        if n_amide > int(cfg.get("peptide_max_amide_bonds", 4)):
            return True
        biggest = max((len(r) for r in mol.GetRingInfo().AtomRings()), default=0)
        if biggest > int(cfg.get("peptide_max_ring_size", 12)):
            return True
    except Exception:
        return False
    return False


def morgan_fp(smiles: str, cfg: dict | None = None):
    """Morgan/ECFP4 bit vector, or ``None``. One definition, used by every step."""
    cfg = cfg or SHARED
    Chem = _rdkit()
    from rdkit.Chem import rdFingerprintGenerator
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        gen = rdFingerprintGenerator.GetMorganGenerator(
            radius=int(cfg.get("fp_radius", 2)), fpSize=int(cfg.get("fp_bits", 2048)))
        return gen.GetFingerprint(mol)
    except Exception:
        return None


def max_tanimoto_to_training(candidate_smiles, training_smiles, cfg: dict | None = None):
    """Nearest-neighbour Tanimoto of each candidate to the training set.

    This is the applicability-domain check the pipeline was missing. A random forest
    given an out-of-domain input does not fail -- it returns the training mean, which
    then looks like a confident prediction.
    """
    cfg = cfg or SHARED
    from rdkit import DataStructs
    train_fps = [fp for fp in (morgan_fp(s, cfg) for s in training_smiles) if fp is not None]
    out = []
    for smi in candidate_smiles:
        fp = morgan_fp(smi, cfg)
        if fp is None or not train_fps:
            out.append(None)
            continue
        out.append(max(DataStructs.BulkTanimotoSimilarity(fp, train_fps)))
    return out
