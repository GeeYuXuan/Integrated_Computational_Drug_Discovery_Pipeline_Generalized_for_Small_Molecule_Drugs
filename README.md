# 🧬 Integrated Computational Drug Discovery Pipeline
### Generalized for small-molecule inhibitors

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GeeYuXuan/Integrated_Computational_Drug_Discovery_Pipeline_Generalized_for_Small_Molecule_Drugs/blob/main/Integrated_Computational_Drug_Discovery_Pipeline_Generalized_for_Small_Molecule_Drugs.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An end-to-end *in silico* pipeline that retrieves bioactivity data for a chosen
protein target, analyses structure–activity relationships, generates novel analogs,
and evaluates them by molecular docking, machine learning and multi-parameter
optimisation.

> ### ⚠️ Read this before using any number from this pipeline
> Everything here is a **computational hypothesis**, not a measurement. On the KRAS
> re-run the potency model reaches an R² of 0.61 on a random split and 0.61 ± 0.10
> across ten scaffold splits; the earlier EGFR and mTOR runs, trained on a truncated
> potency range, reached 0.3–0.5. Docking scores carry an error of about
> 2.85 kcal/mol — larger than the spread between most candidates a run produces, in
> which case the pipeline withholds the term rather than ranking on it. Candidates
> are starting points for synthesis and assay — they are not predicted drugs. See
> [Limitations](#️-limitations) for the specific things this pipeline cannot tell you.

---

## 🚀 Quick start

1. Open the notebook in Colab with the badge above.
2. Run **Step 0** (installs dependencies, once per session).
3. Run **Step 0b** and set one line:
   ```python
   TARGET = "kras"      # "kras" | "egfr" | "mtor"
   ```
4. Run the remaining cells in order.

That is the whole configuration. **Switching target means editing that one line** —
nothing else. All three published runs (KRAS G12D, EGFR, mTOR) execute identical code.

---

## ⚙️ Configuration

Every target-specific value lives in one dictionary in
[`pipeline_common.py`](pipeline_common.py). No notebook cell contains a target name,
a PDB ID, a mutation or a threshold.

```python
TARGETS = {
    "kras": {
        "target_name":    "kras",           # prefix for every output file
        "search_term":    "KRAS",           # ChEMBL / UniProt / PubMed search term
        "chembl_id":      "CHEMBL2189121",  # GTPase KRas, UniProt P01116
        "mutant_filter":  "G12D",           # None for wild type
        "manual_pdb_id":  "7RT1",           # None => automatic RCSB search
        "ligand_resname": "7L8",            # reference ligand that defines the pocket
        "ligand_chain":   "A",              # pin one copy (7RT1 has two)
        "active_site_residue": 12,          # fallback only, never silently used
        "box_center":     None,             # explicit [x, y, z] overrides everything
    },
    ...
}
```

### Adding a new target

Add an entry to `TARGETS` and select it. The fields you need:

| Field | How to find it | Required |
|---|---|---|
| `chembl_id` | Search [ChEMBL](https://www.ebi.ac.uk/chembl/) for your protein, take the **single protein** ID | Recommended — avoids ambiguous name matches |
| `search_term` | The gene or protein name | Yes |
| `target_name` | A short, filesystem-safe label | Yes |
| `mutant_filter` | Mutation code (`"G12D"`, `"T790M"`) or `None` | Yes |
| `manual_pdb_id` | Pick from the structures Step 3 reports | Recommended |
| `ligand_resname` | The 3-letter ligand code in that structure, from [RCSB](https://www.rcsb.org/) | Strongly recommended — see below |
| `active_site_residue` | Residue number of the active site | Fallback only |

Settings shared by every target — potency cutoff, box size, exhaustiveness, scoring
weights, random seed — are in the `SHARED` dictionary in the same file.

### Why `ligand_resname` matters

The docking box is centred on **one named reference ligand**. If you leave it as
`None`, the pipeline auto-detects the largest drug-like hetero group and tells you
what it chose. It will never average several hetero groups together, and it raises
an error rather than falling back to the geometric centre of the protein.

If your structure has no drug-like ligand (only a cofactor such as ADP or GDP), the
pipeline centres on the cofactor site **and says so loudly**. That is correct for an
ATP- or GTP-competitive inhibitor and wrong for an allosteric one — check before
trusting the scores.

---

## 📁 Repository layout

```
├── Integrated_..._Small_Molecule_Drugs.ipynb   the pipeline (12 cells)
├── pipeline_common.py                          shared config + helpers  ← edit this
├── requirements.txt                             pinned dependency floors
├── tests/
│   ├── data/                                    KRAS, EGFR and mTOR structures
│   ├── test_docking_box.py                      box-centring regression test
│   └── test_pipeline_common.py                  config + scoring guard tests
├── Computational Discovery of Novel KRAS Inhibitors.pdf   write-up of the KRAS run
└── kras_project_backup_16Nov26/                 example outputs — KRAS G12D
```

The backup folder holds one complete run:

```
├── workflow/data/    the PDB structure and the prepared receptor
├── results/          CSVs, the trained model, the summary JSON, docked poses
└── figures/          every plot
```

`kras_project_backup_16Nov26/` is the reference run for the current code. It was
produced after every fix listed in [What changed](#-what-changed-and-why), and the
figures and tables in the PDF write-up come from it. It also carries
`results/kras_measurements.csv`, the measurement-level table behind the provenance
figures, so the inter-assay spread can be recomputed from the deposit alone.

The earlier EGFR and mTOR example outputs have been **removed**: they predated the
box, potency-cutoff and training-range fixes, so every number in them was
superseded — the mTOR search box alone sat 37.7 Å from the binding site. Their
crystal structures are kept in `tests/data/`, because those are RCSB inputs rather
than pipeline results and the box regression still checks all three. To regenerate
either target, set `ACTIVE_TARGET` in `pipeline_common.py` and run the notebook.

### Running the tests

No RDKit or Vina needed — the tests use only the standard library, the structures in
`tests/data/` and the committed KRAS poses. Run both suites at once:

```bash
python -m pytest tests/ -q
```

Or run either on its own, which prints every individual check:

```bash
python tests/test_docking_box.py
```

```bash
python tests/test_pipeline_common.py
```

```bash
python tests/test_review_fixes.py
```

`test_review_fixes.py` needs Biopython for its pocket check and skips that one
check without it.

---

## 🔬 Pipeline steps

### Step 0 · Environment
Installs OpenBabel, AutoDock Vina, RDKit, scikit-learn and the ChEMBL client.

### Step 0b · Configuration
Loads `pipeline_common.py` and selects the target. Writes `pipeline_config.json`, so
any later step can be re-run on its own after a kernel restart.

### Step 1 · Data acquisition
Retrieves IC₅₀ data from ChEMBL, converts units, filters to binding/functional assays
with exact relations, applies the mutant filter, standardises structures, removes
salts and peptides, and collapses repeat measurements by **median**.

Reports assay provenance and warns when one assay dominates the dataset.

**Outputs**
- `{target}_measurements.csv` — the measurement-level table, written **before** the medians
  are taken. One row per ChEMBL measurement, so the inter-assay disagreement the step reports
  can be recomputed by anyone holding only the results. Without it the aggregation is a
  one-way door.
- `{target}_bioactivity_full.csv` — all curated compounds, full potency range → used for ML training
- `{target}_inhibitors_cleaned.csv` — the potent subset (IC₅₀ < `potency_cutoff_nm`) → used for SAR and seeds
- `{target}_pic50_distribution.png`

Both aggregate files carry `ic50_min_nm` and `ic50_max_nm` alongside the median, so the
per-compound spread is available without loading the measurement table.

### Step 2 · SAR and scaffold analysis
Computes descriptors (MW, LogP, TPSA, HBD, HBA, QED), performs Bemis–Murcko scaffold
decomposition, and visualises chemical space by PCA on Morgan fingerprints.

**Outputs** `{target}_scaffold_analysis.csv`, `{target}_sar_analysis.csv`,
`{target}_chemical_space_pca.png`, `{target}_physicochemical_properties.png`,
`{target}_top_scaffolds_potency.png`

### Step 3 · Genomic and structural context
Resolves the UniProt ID (ChEMBL first, UniProt API as fallback), searches RCSB for
structures, analyses the chosen structure's binding pocket, and mines PubMed.

**Outputs** `{target}_structural_analysis.json`, `{target}_literature_findings.txt`

> Set `entrez_email` in `pipeline_common.py` to your own address. NCBI throttles or
> blocks queries with an invalid contact address.

### Step 4 · Generative design
Selects seed scaffolds using a hybrid strategy — recurring scaffolds (≥
`min_scaffold_count` compounds) plus high-potency rare ones — then generates analogs
by atom substitution, functional-group addition, ring modification and substituent
addition. Filters on MW, LogP, QED and synthetic accessibility.

Novelty is checked against known molecules **and** known scaffolds, by InChIKey.

**Outputs** `{target}_selected_seeds.csv`, `{target}_generated_candidates.csv`,
`{target}_top20_generated_candidates.csv`, `{target}_generation_pca.png`,
`{target}_generation_tsne.png`

### Step 5 · Virtual screening (docking)
Locates the pocket, builds a receptor that **keeps structurally integral cofactors**,
grows the box to fit the reference ligand, and docks with AutoDock Vina.

Reports **ligand efficiency** alongside raw affinity, and states plainly when the
spread across candidates is smaller than Vina's own error.

**Outputs** `{target}_docking_results.csv`, `{target}_docking_scores.png`,
`{target}_docking_site.json` (an audit record of how the box was chosen),
`candidate_*_docked.pdbqt`

### Step 6 · Machine-learning prediction
Trains a Random Forest on Morgan fingerprints over the **full measured potency
range**. Reports a random split, **repeated scaffold splits** (mean ± sd, the
honest estimate) and a mean-only baseline, and scores every candidate for
**applicability domain**.

**Outputs** `{target}_candidate_predictions.csv`, `{target}_pIC50_model.pkl`,
`{target}_model_metrics.json`, `{target}_model_performance.png`,
`{target}_ml_docking_comparison.png`

### Step 7 · Physicochemical profiling and final selection
Computes physicochemical descriptors, Lipinski and Veber compliance, and PAINS
alerts, then combines them into a consensus score using **fixed** bounds.

> This step does **not** predict ADME or toxicity. There is no BBB, Caco-2, CYP,
> hERG or hepatotoxicity model here. For real ADMET endpoints, run the final
> candidates through [ADMETlab](https://admetlab3.scbdd.com/),
> [pkCSM](https://biosig.lab.uq.edu.au/pkcsm/) or
> [SwissADME](http://www.swissadme.ch/).

**Outputs** `{target}_admet_analysis.csv`, `{target}_final_candidates.csv`,
`{target}_candidate_radar_plot.png`

### Step 8 · Consolidation
Aggregates every step's metrics into `{target}_project_summary_metrics.json`, and
verifies that each expected artifact exists.

---

## ⚠️ Limitations

**The potency model is modest, and a single split does not pin it down.** The
January runs reached R² of 0.29 (KRAS), 0.41 (EGFR) and 0.46 (mTOR) on random
splits. Training across the full potency range lifted KRAS to **0.61** on a live
re-run. Treat a predicted IC₅₀ as "plausibly active", never as a number.

Scaffold splits — holding out whole Bemis–Murcko scaffolds — are the harder and
more honest test, but **one draw is not an estimate**. The first live KRAS run
drew a scaffold split of 0.794 against a random split of 0.614, which is the wrong
way round and was mostly the draw. Step 6 therefore repeats the split
`scaffold_split_repeats` times (default 10) and reports mean ± sd, flagging the
estimate as unstable when the repeats disagree. Quote the mean with its spread;
never quote the best draw.

**Predictions on generated candidates may carry no signal.** Generation starts from
bare scaffolds, which sit outside the chemistry the model was trained on. A random
forest given an out-of-domain molecule returns the training mean, which looks like a
confident prediction. Step 6 measures this and refuses to rank when the prediction
spread falls below the model's own error. In the original runs it always did.

**Docking is a single rigid-receptor scoring function.** No protein flexibility, no
explicit water, no free-energy calculation. Vina's error (~2.85 kcal/mol) is larger
than the differences between most candidates here.

**Scaffold conclusions reflect dataset composition.** The KRAS G12D set draws 68% of
its measurements from a single assay, so its "privileged scaffolds" largely describe
one compound series. Step 1 reports this.

**The mutant filter is a text match** on the assay description, so an assay comparing
your mutation against another variant also matches.

**Generated analogs are not synthesis routes.** The SA score estimates difficulty; it
does not propose a route, and no retrosynthetic check is performed.

**This is not a validated drug-discovery workflow.** It is a learning and
hypothesis-generation project.

---

## 🔧 What changed, and why

### Second review, 17 November 2026

A peer review of the KRAS report checked every reported figure against the
deposited outputs and the code behind them. Seven code defects came out of it.
None changes a docking affinity or a model metric, but five change numbers that
were reported, and the report has been corrected to match.

| Issue | Effect | Fix |
|---|---|---|
| **Only the first assay id survived aggregation** | Every "unique assays" figure was an undercount — 16 reported against 25 real for the curated set, 11 against 17 for the potent subset — and the dominant-assay share was computed over compounds rather than measurements | Carry `n_assays` per compound and take any dataset-level count from `{target}_measurements.csv` ([`test_review_fixes.py`](tests/test_review_fixes.py)) |
| **Pocket residues counted waters, around the wrong ligand** | Biopython marks waters `W`, not `H_`, so they passed the protein filter; and `ligands[0]` is GDP in 7RT1, not the inhibitor. "37 pocket residues" was 24 protein residues + 13 waters around a different site | Pin the reference to `ligand_resname`, exclude waters, report both counts |
| **Resolution parser returned the remark number** | `REMARK   2 RESOLUTION.    1.27 ANGSTROMS.` was read as 2.0 Å, for every structure | Take the token after `RESOLUTION`, never the first float on the line |
| **Receptor converted with no hydrogens and no charges** | The deposited `7RT1_receptor.pdbqt` has 1,367 atoms, 0 hydrogens and every partial charge 0.000, while the Methods claimed protonation at pH 7.4 with assigned charges | Pass `-p 7.4`, then verify the written file and say so loudly if it still has no hydrogens. Charges are a separate problem — see below |
| **Novelty measured against the potent subset only** | A generated molecule identical to one of the weaker curated compounds was reported as novel | Build the reference from `{target}_bioactivity_full.csv` |
| **Applicability domain measured against all 552** | Includes the 20 % the model never saw, so the domain looks wider than it is | Compare against the training partition |
| **Effective weights and dropped terms never exported** | `Scoring_Dropped_Terms` was computed and printed but written to no file, so the CSVs gave no sign that two of five terms had been withheld | Export the column with both tables and write `{target}_scoring_weights.json` |

`TSNE(n_iter=)` was also replaced with `max_iter` (renamed in scikit-learn 1.5,
removed in 1.7).

### What the verification re-run found

The whole pipeline was then re-run on the fixed code to check the fixes in
place. Three portability bugs surfaced that Colab had hidden, all now fixed:
Step 4 shelled out to `wget` (absent on Windows and macOS) for the SA-score
table; eight `open()` calls had no explicit `encoding`, so Step 3 died writing
PubMed abstracts under a cp1252 locale; and `--partialcharge gasteiger` turned
out to make Open Babel emit **zero molecules**.

Measured against Open Babel 3.1.0 on `7RT1_clean.pdb`:

| obabel arguments | atoms | hydrogens | non-zero charges |
|---|---|---|---|
| `-xr` (as deposited) | 1367 | 0 | 0 |
| `-xr -p 7.4` (now) | 1662 | 295 | 0 |
| `-xr --partialcharge gasteiger` | — | — | conversion fails |

So `-p 7.4` does add the polar hydrogens, but **this Open Babel's PDBQT writer
emits `+0.000` for every atom whatever it is asked**. Real partial charges need
`prepare_receptor4.py` (ADFR), which this pipeline does not run. The report
says so rather than implying otherwise.

**Docking reproducibility.** Re-docking the same 20 ligands against the *same*
unprotonated receptor with a different Vina build (1.2.7 vs Colab's apt
package) moved the mean by −0.010 kcal/mol but individual scores by up to
**0.45 kcal/mol**, reshuffling **25% of all pairwise rankings**. The candidate
set's entire spread is 0.83 kcal/mol. Adding the 295 polar hydrogens, by
contrast, was a near-uniform −0.364 kcal/mol offset that left 85% of pairwise
rankings intact. Run-to-run search noise, not the receptor, is what makes these
affinities unrankable — measured on this system rather than argued from Vina's
published ~2.85 kcal/mol error.

**Everything that matters reproduced.** Steps 1-4 are identical to the deposit
to within 1e-14 (897 measurements, 552 compounds, spread 1.92/3.79, 41
scaffolds, 15 seeds, 70 candidates, the same top 20). The final five
candidates, their order and their consensus scores are unchanged;
`Consensus_Score = 0.5 x QED + 0.5` still holds exactly; Potency and Docking
are still both dropped as noise; and all 20 candidates remain inside the
applicability domain even under the stricter training-partition test.

**The deposited `kras_project_backup_16Nov26/` predates these fixes.** Its
docking affinities came from the unprotonated, uncharged receptor; its
`structural_analysis.json` holds the 37-residue count; its CSVs carry no
dropped-terms column. The figures the report quotes are the deposited ones,
which is why the report states the receptor limitation rather than hiding it. A
re-run on the fixed code will differ in the pocket count, the assay counts, the
applicability-domain margins and the exported columns, and may differ slightly
in affinity.

### First review, November 2026

A review in November 2026 found four issues serious enough to invalidate results.
All are fixed; each fix is covered by a test.

| Issue | Effect | Fix |
|---|---|---|
| **Docking box centred on the mean of all hetero atoms** | Worked for EGFR (one ligand). KRAS was 7.3 Å off; mTOR was **37.7 Å off, with zero pocket inside the box** — every mTOR docking score was meaningless | Group HETATM by residue, centre on one named reference ligand, never average; raise rather than guess ([`test_docking_box.py`](tests/test_docking_box.py)) |
| **`POTENCY_CUTOFF_NM` was ignored** | The documented knob did nothing; the filter was hard-coded to 50 nM while the plot title showed the configured value | The config value drives the filter |
| **Novelty checked only against whole molecules** | Generated structures reproducing a known scaffold passed as "novel" — the rank-1 mTOR candidate was an existing scaffold | Check molecules *and* scaffolds, by InChIKey |
| **Model trained only on sub-50 nM compounds** | Truncated label range; candidate predictions varied 3–11× *less* than the model's own RMSE, so ranking ranked noise | Train on the full measured range; add an applicability-domain check; refuse to rank a spread below the model error |

Also fixed: cofactors are no longer stripped from the receptor; the box grows to fit
the reference ligand; repeat measurements are aggregated by median rather than best
value; peptides are excluded; consensus scoring uses fixed bounds instead of
per-batch min–max; ligand efficiency is reported; the model and two figures are
target-prefixed (they previously overwrote each other between runs);
`unique_scaffolds` counted compounds rather than scaffolds; summary JSON contained
bare `NaN`, which strict parsers reject; the completion date was hard-coded;
conformer generation is seeded; and all configuration now lives in one place.

---

## 📋 Requirements

See [`requirements.txt`](requirements.txt). The pipeline needs OpenBabel and AutoDock
Vina as system binaries, installed by Step 0 on Colab via `apt-get`.

RDKit must be 2023.03 or newer — earlier versions lack `rdFingerprintGenerator`,
which Steps 2, 4 and 6 use.

---

## 📚 Citing the tools used

If you build on this, please cite the underlying tools:

- **AutoDock Vina** — Eberhardt J., Santos-Martins D., Tillack A.F., Forli S. (2021) *J. Chem. Inf. Model.* 61:3891–3898
- **RDKit** — Landrum G., [rdkit.org](https://www.rdkit.org/)
- **ChEMBL** — Zdrazil B. *et al.* (2024) *Nucleic Acids Res.* 52:D1180–D1192
- **Synthetic accessibility score** — Ertl P. & Schuffenhauer A. (2009) *J. Cheminform.* 1:8
- **QED** — Bickerton G.R. *et al.* (2012) *Nat. Chem.* 4:90–98
- **OpenBabel** — O'Boyle N.M. *et al.* (2011) *J. Cheminform.* 3:33

---

## 🩺 Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `Failed building wheel for vina` | The `vina` **Python package** has no wheel for Python 3.13 (Colab's current interpreter). The pipeline does not use it — docking drives the `vina` **executable** from `apt install autodock-vina`. If you see this, you are running an older copy of Step 0 that pip-installed `vina`; pull the latest notebook |
| `AutoDock Vina executable not found` | Run `!apt-get install -y autodock-vina`, or set `vina_binary` in `pipeline_common.py` to its full path |
| `Install failed for: ['rdkit', 'vina', 'Bio']` | One failing build aborts a whole `pip install`. Step 0 now installs each package separately so this cannot cascade — pull the latest notebook |
| `Could not download pipeline_common.py` | Upload `pipeline_common.py` from this repository into `/content`, then re-run Step 0b |
| `Reference ligand 'XXX' is not present` | The structure does not contain that ligand. The error lists what *is* present — pick one, or set `ligand_resname` to `None` for auto-detection |
| `Could not locate a binding site` | No drug-like ligand, no cofactor, no usable `active_site_residue`. Set `ligand_resname` or `box_center` explicitly rather than letting it guess |
| `No compounds below N nM` | Raise `potency_cutoff_nm` (1000.0 is reasonable for a hard target) |
| "Centring on a COFACTOR site" warning | The structure has no inhibitor bound. Fine for an ATP/GTP-competitive series; wrong for an allosteric one |
| "the candidate ORDER is not meaningful" | Working as intended — the prediction spread is below the model's error. Report the candidates as an unordered set |
| `ModuleNotFoundError` | Re-run Step 0 |
| Step 5 cannot find a PDB | Set `manual_pdb_id`, or place `{PDB_ID}.pdb` in `workflow/data/` by hand |

---

## 📄 License

MIT — see [LICENSE](LICENSE).
