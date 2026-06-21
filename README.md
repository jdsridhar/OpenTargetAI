# OpenTargetAI

**Open Source Ligand-Based Target Prediction Platform**

OpenTargetAI is a desktop application for **ligand-based protein target
prediction**. Given a small molecule, it finds chemically similar compounds in a
local bioactivity database (built from public ChEMBL / BindingDB data),
aggregates their known targets, and produces a ranked, **explainable** list of
candidate targets with transparent confidence scores.

> ⚠️ **Scientific scope.** OpenTargetAI is an independent, open-source tool for
> research and education. It is **not** a reproduction of SwissTargetPrediction
> and does not claim to replicate its model or results. Predictions are
> hypotheses derived from chemical similarity and must not be used for clinical,
> regulatory, or safety decisions.

---

## Key features

- **Modern PyQt6 desktop UI** with a dark scientific theme.
- **Molecule input** by SMILES paste, SDF/MOL file, worked examples, or an
  in-app molecule builder; live 2D depiction and full physicochemical /
  drug-likeness profile (Lipinski, Veber, Ghose, Lead-likeness, QED, PAINS).
- **Five fingerprint types** (Morgan ECFP4/ECFP6, MACCS, Atom Pair, Topological)
  via RDKit's modern generator API, stored as compact binary blobs.
- **Three similarity metrics** (Tanimoto, Dice, Cosine) with fast bulk search.
- **Weighted target-fishing engine**: scores each candidate target from
  similarity strength × bioactivity potency × occurrence frequency, with a
  separate, tunable **confidence model** (High / Medium / Low).
- **Explainable predictions**: per-target evidence breakdown, supporting
  ligands, shared Morgan substructure highlighting, Murcko scaffold and MCS.
- **Interactive network** (query → targets → ligands) with zoom/pan.
- **Advanced charts**: similarity / confidence / activity distributions, target
  class pie chart, plus a standalone interactive Plotly dashboard.
- **PDB integration**: retrieve experimental structures (PDBe/RCSB) for a
  predicted target and download coordinate files.
- **Docking preparation**: build a ready-to-dock AutoDock Vina package (3D
  ligand, receptor, suggested grid box, config).
- **Export** to CSV, Excel, HTML and PDF.
- **Validation module**: leave-one-out benchmark reporting Top-1/5/10 accuracy,
  MRR, Precision@10, Recall@10 and ROC-AUC.
- **Standalone CLI importer** and a SQLite backend designed to scale to millions
  of compounds (WAL mode, indexes, in-memory fingerprint cache, background
  workers).

---

## Project structure

```
OpenTargerAi/
├── app.py                     # GUI entry point
├── import_data.py             # standalone CLI database importer
├── schema.sql                 # SQLite DDL (documentation / manual setup)
├── requirements.txt
├── README.md
│
├── database/                  # storage layer
│   ├── db_manager.py          # schema, CRUD, caching, statistics
│   ├── importer.py            # ChEMBL / BindingDB / demo importers
│   └── demo_data.py           # curated demo dataset
├── similarity/                # cheminformatics
│   ├── fingerprint_engine.py  # fingerprint generation & similarity
│   ├── search_engine.py       # similarity search over the database
│   └── mol_utils.py           # descriptors, drug-likeness, scaffolds, 2D SVG
├── prediction/                # target prediction
│   ├── prediction_engine.py   # weighted target fishing
│   └── confidence.py          # confidence model
├── benchmark/
│   └── validator.py           # leave-one-out benchmark + metrics
├── visualization/
│   ├── network_view.py        # PyQtGraph + NetworkX target network
│   └── charts.py              # charts + Plotly HTML export
├── integration/
│   ├── pdb_client.py          # PDBe / RCSB structure retrieval
│   └── docking_prep.py        # AutoDock Vina package builder
├── reports/
│   └── exporter.py            # CSV / Excel / HTML / PDF exports
├── workers/
│   └── workers.py             # background QThread workers
└── ui/                        # PyQt6 panels & dialogs
    ├── main_window.py
    ├── molecule_input.py      # left panel
    ├── molecule_viewer.py     # center panel
    ├── prediction_controls.py # right panel
    ├── results_view.py
    ├── target_details.py
    ├── explainer_view.py
    ├── charts_view.py
    ├── pdb_view.py
    ├── validation_view.py
    ├── dialogs.py
    └── widgets.py
```

---

## Installation

Requires **Python 3.12+** (tested through 3.14).

```bash
# 1. (recommended) create a virtual environment
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS / Linux:
source .venv/bin/activate

# 2. install dependencies
pip install -r requirements.txt
```

All dependencies install from PyPI, including RDKit (`pip install rdkit`).

---

## Quick start

```bash
# Launch the desktop application
python app.py
```

On first launch the database is empty and you'll be offered the **curated demo
dataset** (37 well-known drugs across 11 targets). Accept it, then:

1. Go to the **Predict** tab and paste a SMILES (or pick an example), e.g.
   ibuprofen `CC(C)Cc1ccc(C(C)C(=O)O)cc1`.
2. Adjust the fingerprint, metric, similarity threshold and Top-N on the right.
3. Click **Run Prediction**.
4. Explore the **Results**, **Network**, **Charts**, **Structures & Docking**
   and **Validation** tabs.

### Loading data without the GUI

```bash
# Load the demo dataset into the default database (data/opentargetai.db)
python import_data.py --demo

# Import a full ChEMBL activities TSV export
python import_data.py --chembl path/to/chembl_activities.tsv

# Import a BindingDB TSV export
python import_data.py --bindingdb path/to/BindingDB_All.tsv

# Inspect the database
python import_data.py --stats
```

### Expected import columns

**ChEMBL TSV** (tab-separated): `canonical_smiles`, `chembl_id`,
`target_chembl_id`, `uniprot_id`, `gene_name`, `target_name`, `organism`,
`standard_type`, `standard_value`, `standard_units`, `pchembl_value`,
`assay_type`. Such exports can be produced from the ChEMBL web interface or the
ChEMBL SQL/Postgres dump.

**BindingDB TSV**: the standard BindingDB download columns, including
`Ligand SMILES`, `Target Name`, `UniProt (SwissProt) Primary ID of Target
Chain`, `Ki (nM)` / `IC50 (nM)` / `Kd (nM)` / `EC50 (nM)`, and
`Target Source Organism …`.

Rows with invalid SMILES or missing target names are skipped; compounds are
deduplicated by InChIKey and all five fingerprints are generated on import.

---

## How it works

1. **Fingerprint the query** with the selected representation.
2. **Similarity search** — compute bulk Tanimoto/Dice/Cosine against every stored
   compound and keep neighbours above the threshold (Top-N).
3. **Target fishing** — collect every known target of those neighbours.
4. **Scoring** — for each target:
   `target_score = similarity_weight × bioactivity_weight × occurrence_frequency`
   where *similarity_weight* is the mean squared neighbour similarity,
   *bioactivity_weight* is derived from mean pChEMBL potency, and
   *occurrence_frequency* is the fraction of neighbours that hit the target.
5. **Confidence** — a weighted blend of similarity strength, potency, recurrence
   and activity consistency, mapped to High / Medium / Low.
6. **Explain** — every prediction exposes its supporting ligands, the shared
   substructure with the query, and its confidence components.

---

## Validation

The **Validation** tab (and `benchmark/validator.py`) runs leave-one-out
cross-validation: each compound is masked from its own neighbour set, its
targets are predicted, and the ranking is compared with its experimentally known
targets. Reported metrics: Top-1/5/10 accuracy, Mean Reciprocal Rank,
Precision@10, Recall@10 and a pooled ROC-AUC. Results depend heavily on the
density and diversity of the imported dataset.

---

## Structures & docking

For a selected target with a UniProt accession, OpenTargetAI retrieves
experimental structures from PDBe/RCSB and can download coordinate files. The
docking preparation step assembles an AutoDock Vina package: a 3D-embedded,
energy-minimised ligand, the chosen receptor, a suggested search box (centred on
a co-crystallised ligand when present), a `vina_config.txt`, and a README with
the remaining preparation steps. **AutoDock Vina / AutoDockTools are not
bundled** — install them separately to run the docking.

---

## Performance notes

- SQLite is opened in WAL mode with tuned pragmas and indexed lookups.
- Fingerprints are stored as compact native binary blobs and cached in memory
  after the first search of a session.
- All long-running work (search, import, validation, network/file IO) runs in
  background `QThread` workers so the UI stays responsive.

---

## Disclaimer

OpenTargetAI is provided "as is" for research and educational purposes. The
bundled demo dataset uses rounded, representative potency values and must not be
used for real decision making. Predictions are similarity-based hypotheses, not
validated facts.

## License

Released under the MIT License. Public data sources (ChEMBL, BindingDB, RCSB
PDB, UniProt) are subject to their own respective licenses and terms of use.
#   O p e n T a r g e t A I  
 