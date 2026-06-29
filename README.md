<div align="center">
  <img src="resources/banner.png" alt="OpenTargetAI Banner" width="100%" style="border-radius: 8px; margin-bottom: 20px;" />
  
  # OpenTargetAI
  
  ### 🔬 Open-Source Ligand-Based Target Prediction & Cheminformatics Platform
  
  [![Python Version](https://img.shields.io/badge/Python-3.12%20%7C%203.13%20%7C%203.14-blue.svg?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
  [![PyQt6](https://img.shields.io/badge/GUI-PyQt6-brightgreen.svg?style=flat-square&logo=qt&logoColor=white)](https://www.riverbankcomputing.com/software/pyqt/)
  [![RDKit](https://img.shields.io/badge/Cheminformatics-RDKit-orange.svg?style=flat-square)](https://www.rdkit.org/)
  [![SQLite](https://img.shields.io/badge/Database-SQLite-003B57.svg?style=flat-square&logo=sqlite&logoColor=white)](https://sqlite.org/)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](https://opensource.org/licenses/MIT)
  [![PRs Welcome](https://img.shields.io/badge/PRs-welcome-violet.svg?style=flat-square)](https://github.com/jdsridhar/OpenTargetAI/pulls)

  <p align="center">
    <strong>An explainable, local target-fishing workstation designed for rapid drug discovery screening and chemical biology workflows.</strong>
  </p>
</div>

---

> [!WARNING]  
> **Scientific Scope Disclaimer**  
> OpenTargetAI is an independent, open-source tool for research and education. It is **not** a reproduction of SwissTargetPrediction and does not claim to replicate its model or results. Predictions are hypotheses derived from chemical similarity and must not be used for clinical, regulatory, or safety decisions.

---

## 📌 Table of Contents
- [🚀 Key Features](#-key-features)
- [📂 Project Structure](#-project-structure)
- [⚙️ Installation](#%EF%B8%8F-installation)
- [⏱️ Quick Start](#%EF%B8%8F-quick-start)
- [🔬 How It Works](#-how-it-works)
- [📊 Validation](#-validation)
- [🧬 Structures & Docking](#-structures--docking)
- [⚡ Performance Notes](#-performance-notes)
- [⚖️ Disclaimer & License](#%EF%B8%8F-disclaimer--license)

---

## 🚀 Key Features

*   **🖥️ PyQt6 Desktop UI** – Premium dark scientific theme featuring an in-app molecule builder, SMILES pasting, SDF/MOL import, and live 2D chemical structure rendering.
*   **🩺 Physicochemical Profiling** – Real-time estimation of drug-likeness rules (Lipinski, Veber, Ghose, Lead-likeness, QED) and structural flags (PAINS).
*   **🧬 Advanced Cheminformatics** – Compute five distinct molecular representations (Morgan ECFP4/ECFP6, MACCS, Atom Pair, Topological) via RDKit’s modern generator API, stored as compact binary database blobs.
*   **⚡ Fast Similarity Engine** – Perform bulk similarity matching over millions of compounds using Tanimoto, Dice, or Cosine metrics.
*   **🎣 Target-Fishing Engine** – Custom scoring algorithm combining similarity coefficients, bioactivity potencies, and target occurrences into a consolidated confidence rating (High / Medium / Low).
*   **💡 Explainable AI (XAI)** – Trace the evidence for any prediction with detailed target breakdowns, supporting reference ligands, maximum common substructures (MCS), and highlighting shared Morgan scaffolds.
*   **🌐 Interactive Network View** – Zoomable and pannable interactome graph showing target-ligand association networks (powered by PyQtGraph & NetworkX).
*   **🔌 Structure & Docking Prep** – Automatically download target experimental 3D structures from PDBe/RCSB and generate complete, ready-to-run AutoDock Vina packages (3D minimized ligand, receptor, suggested grid box, and configuration).
*   **📊 Benchmark Reporting** – Evaluate predictions locally with a leave-one-out cross-validation suite reporting Top-1/5/10 accuracy, MRR, Precision, Recall, and ROC-AUC.
*   **💾 Enterprise Database Layer** – SQLite database running in WAL mode with background `QThread` workers, custom pragmas, and in-memory caches to maintain seamless UI interactivity.

---

## 📂 Project Structure

```
OpenTargerAi/
├── app.py                     # GUI desktop application entry point
├── import_data.py             # CLI database importer (ChEMBL / BindingDB / Demo)
├── schema.sql                 # SQLite DDL schema for documentation & manual setup
├── requirements.txt           # Package dependencies
├── README.md                  # Beautiful platform documentation
│
├── database/                  # Storage Layer
│   ├── db_manager.py          # SQLite schema, CRUD operations, caches & statistics
│   ├── importer.py            # Custom CLI and demo dataset parsers
│   └── demo_data.py           # Preloaded small-scale reference dataset
│
├── similarity/                # Cheminformatics Core
│   ├── fingerprint_engine.py  # RDKit fingerprint generators and binary encoders
│   ├── search_engine.py       # High-performance bulk database search
│   └── mol_utils.py           # Molecular descriptors, scaffolds, & SVG renderings
│
├── prediction/                # Prediction Mechanics
│   ├── prediction_engine.py   # Weighted target fishing logic
│   └── confidence.py          # Multiclass target confidence model
│
├── benchmark/                 # Evaluation
│   └── validator.py           # Cross-validation testing suite & metrics
│
├── visualization/             # Graphics and UI Plots
│   ├── network_view.py        # Interactome graph (PyQtGraph + NetworkX)
│   └── charts.py              # Distribution plots and standalone interactive Plotly export
│
├── integration/               # External APIs
│   ├── pdb_client.py          # RCSB/PDBe structure and metadata downloaders
│   └── docking_prep.py        # AutoDock Vina project package generator
│
├── reports/                   # Export Engine
│   └── exporter.py            # CSV, Excel, HTML, and PDF document compilers
│
├── workers/                   # Async Architecture
│   └── workers.py             # Background QThread workers preventing GUI freeze
│
└── ui/                        # PyQt6 View Controller
    ├── main_window.py         # Main viewport shell
    ├── molecule_input.py      # Query input pane (SMILES/builder/SDF)
    ├── molecule_viewer.py     # Descriptors & 2D rendering pane
    ├── prediction_controls.py # Hyperparameter & threshold settings pane
    └── ...                    # Specialized layouts (results, network, chart, validation, etc.)
```

---

## ⚙️ Installation

OpenTargetAI requires **Python 3.12+** (fully verified up to Python 3.14).

```bash
# 1. Create a clean virtual environment
python -m venv .venv

# 2. Activate the virtual environment
# On Windows (PowerShell):
.venv\Scripts\activate
# On macOS / Linux:
source .venv/bin/activate

# 3. Install core dependencies
pip install -r requirements.txt
```

> [!TIP]
> All primary dependencies (including RDKit) install cleanly via pip on modern platforms without needing complex conda environments.

---

## ⏱️ Quick Start

```bash
# Start the application
python app.py
```

### Running Your First Prediction
On your first startup, the database is blank. OpenTargetAI will prompt you to automatically load the **curated demo dataset** (containing 37 reference drugs spanning 11 unique biological targets). Accept the import, then follow these steps:

1. Navigate to the **Predict** tab on the left panel.
2. Paste a query SMILES (e.g., Ibuprofen: `CC(C)Cc1ccc(C(C)C(=O)O)cc1`) or click an example preset.
3. Configure settings on the right panel (Fingerprint type, Distance metric, Top-N neighbors, and Similarity threshold).
4. Click **Run Prediction**.
5. Explore predicted targets across the **Results**, **Network**, **Charts**, **Structures & Docking**, and **Validation** tabs.

### Database Importer CLI
You can inspect or import larger datasets directly from your terminal:

```bash
# Load the demo dataset into the default SQLite database (data/opentargetai.db)
python import_data.py --demo

# Import a custom ChEMBL activity TSV export
python import_data.py --chembl path/to/chembl_activities.tsv

# Import a BindingDB raw TSV export
python import_data.py --bindingdb path/to/BindingDB_All.tsv

# Print database stats (number of compounds, activities, and targets)
python import_data.py --stats
```

---

## 📋 Expected Input Specifications

To populate the database using the CLI importer, TSV files must match the structures below:

### 1. ChEMBL Format (Tab-Separated)
| Column Header | Description |
| :--- | :--- |
| `canonical_smiles` | RDKit-compatible compound SMILES string (required) |
| `chembl_id` | Unique compound identifier (e.g., `CHEMBL521`) |
| `target_chembl_id` | Unique target identifier (e.g., `CHEMBL220`) |
| `uniprot_id` | Swiss-Prot accession code (e.g., `P35968`) |
| `gene_name` | Primary gene symbol |
| `target_name` | Full biological target description |
| `organism` | Species name (e.g., *Homo sapiens*) |
| `standard_type` | Bioactivity type (e.g., `IC50`, `Ki`, `Kd`) |
| `standard_value` | Quantitative value (e.g., `12.5`) |
| `standard_units` | Measurement units (e.g., `nM`) |
| `pchembl_value` | Negative log concentration (e.g., `7.90`) |
| `assay_type` | Assay classification (e.g., `B`, `F`) |

### 2. BindingDB Format (Tab-Separated)
The standard BindingDB public download requires the following columns:
*   `Ligand SMILES`
*   `Target Name`
*   `UniProt (SwissProt) Primary ID of Target Chain`
*   `Ki (nM)` / `IC50 (nM)` / `Kd (nM)` / `EC50 (nM)`
*   `Target Source Organism...`

---

## 🔬 How It Works

OpenTargetAI implements an explainable, data-driven target fishing model:

```
[Query Molecule] ──> [Fingerprint Generator] ──> [Bulk Similarity Search]
                                                           │
                                                           ▼
[Target Explainer] <── [Scoring Engine] <── [Top-N Supporting Neighbors]
```

1. **Fingerprint Encoding**: The query molecule is converted into the selected binary fingerprint vector.
2. **Similarity Neighborhood Search**: The search engine computes distance metrics against every compound in the database, retaining neighbors above the similarity cutoff.
3. **Evidence Aggregation**: All targets annotated to the matching neighbors are retrieved.
4. **Target Scoring**: Each candidate target is ranked by a composite mathematical model:
   
   $$\text{Target Score} = w_{\text{similarity}} \times w_{\text{bioactivity}} \times f_{\text{occurrence}}$$
   
   *   **Similarity Weight ($w_{\text{similarity}}$)**: The mean squared similarity coefficient of the supporting neighbors.
   *   **Bioactivity Weight ($w_{\text{bioactivity}}$)**: A potency scaling coefficient computed from the mean $p\text{ChEMBL}$ value of supporting activities.
   *   **Occurrence Frequency ($f_{\text{occurrence}}$)**: The proportion of retrieved neighbors associated with this target.
5. **Confidence Classification**: Targets are classified as **High**, **Medium**, or **Low** confidence based on neighborhood size, similarity depth, and potency consistency.
6. **Explainability Details**: For each predicted target, RDKit calculates the **Maximum Common Substructure (MCS)** and highlights shared scaffolds to explain *why* the prediction was made.

---

## 📊 Validation

You can evaluate the accuracy of your local bioactivity database using the **Validation** panel or by executing:

```bash
# Run leave-one-out cross-validation
python -m benchmark.validator
```

The validation suite performs a **leave-one-out cross-validation (LOOCV)** benchmark:
*   Each compound in the database is masked from its own target-fishing neighborhood.
*   Its target profiles are predicted using the remaining compound network.
*   Predicted rankings are matched against experimental ground-truth to calculate accuracy metrics:
    - **Top-1 / Top-5 / Top-10 Accuracy**
    - **Mean Reciprocal Rank (MRR)**
    - **Precision@10 & Recall@10**
    - **Area Under the ROC Curve (ROC-AUC)**

---

## 🧬 Structures & Docking

OpenTargetAI bridges ligand-based screening with structure-based verification:

*   **PDB Fetching**: Query the PDBe API to list and download high-resolution crystallographic structure coordinate files (.pdb / .cif) for predicted target accessions.
*   **Vina Docking Preparation**: Auto-prepare docking folders containing:
    1.  A 3D conformer of the query molecule, energy-minimized using RDKit's MMFF94 force field.
    2.  The target protein structure.
    3.  A suggested grid box configuration centered on the active site (defined by co-crystallized reference ligands).
    4.  A `vina_config.txt` control file and README detailing execution commands.

> [!NOTE]  
> AutoDock Vina and preparation command-line tools (like ADFRSuite or Meeko) must be installed separately on your system to run the prepared docking calculations.

---

## ⚡ Performance Notes

*   **SQLite WAL Mode**: SQLite operates in Write-Ahead Logging mode with customized page sizing, synchronous settings, and query indexing to handle datasets containing millions of datapoints without interface lag.
*   **In-Memory Fingerprint Caching**: Binary fingerprints are decompressed and cached in RAM upon the first query, speeding up subsequent predictions by up to 20x.
*   **Async Threading**: Long tasks (database loads, similarity calculations, structure preparation, and benchmarks) run in background `QThread` workers to keep the desktop frame smooth and responsive.

---

## ⚖️ Disclaimer & License

*   **Disclaimer**: OpenTargetAI is provided "as is" for research and educational purposes. The bundled demo dataset uses rounded, representative potency values and must not be used for real decision making.
*   **License**: Released under the [MIT License](LICENSE). Public data sources (ChEMBL, BindingDB, RCSB PDB, UniProt) are subject to their own respective licensing terms.