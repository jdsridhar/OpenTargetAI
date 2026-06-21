"""
Left panel: molecule input for OpenTargetAI.

Lets the user provide a query molecule by pasting SMILES, loading a SMILES/SDF/
MOL file, picking a worked example, or building one in the sketch dialog. Emits
``molecule_changed`` with a canonical SMILES once a valid molecule is parsed.
"""

from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QGridLayout,
    QGroupBox,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from similarity import mol_utils
from ui.widgets import SectionTitle

# A few well-known drugs for one-click loading.
EXAMPLE_MOLECULES = [
    ("— Examples —", ""),
    ("Ibuprofen (NSAID)", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
    ("Imatinib (kinase inhibitor)", "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"),
    ("Sildenafil (PDE5)", "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12"),
    ("Diphenhydramine (H1)", "CN(C)CCOC(c1ccccc1)c1ccccc1"),
    ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("Estradiol", "CC12CCC3c4ccc(O)cc4CCC3C1CCC2O"),
]


class MoleculeInputPanel(QWidget):
    """Query-molecule input controls."""

    molecule_changed = pyqtSignal(str)   # canonical SMILES of a valid molecule
    status_message = pyqtSignal(str)

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._build_ui()

    # ── UI ────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.addWidget(SectionTitle("Input Molecule"))

        # SMILES entry
        smiles_group = QGroupBox("Paste SMILES")
        sg_layout = QVBoxLayout(smiles_group)
        self.smiles_edit = QPlainTextEdit()
        self.smiles_edit.setPlaceholderText(
            "Paste or type a SMILES string, e.g. CC(=O)Oc1ccccc1C(=O)O"
        )
        self.smiles_edit.setMaximumHeight(90)
        sg_layout.addWidget(self.smiles_edit)

        self.example_combo = QComboBox()
        for label, smiles in EXAMPLE_MOLECULES:
            self.example_combo.addItem(label, smiles)
        self.example_combo.currentIndexChanged.connect(self._on_example_selected)
        sg_layout.addWidget(self.example_combo)
        layout.addWidget(smiles_group)

        # Action buttons
        buttons_group = QGroupBox("Actions")
        grid = QGridLayout(buttons_group)
        self.btn_load_smiles = QPushButton("Load SMILES")
        self.btn_validate = QPushButton("Validate Molecule")
        self.btn_load_sdf = QPushButton("Load SDF")
        self.btn_load_mol = QPushButton("Load MOL")
        self.btn_draw = QPushButton("Draw Molecule")
        self.btn_generate = QPushButton("Generate Structure")

        self.btn_load_smiles.clicked.connect(self._load_from_text)
        self.btn_validate.clicked.connect(self._validate)
        self.btn_load_sdf.clicked.connect(lambda: self._load_file("sdf"))
        self.btn_load_mol.clicked.connect(lambda: self._load_file("mol"))
        self.btn_draw.clicked.connect(self._open_sketch)
        self.btn_generate.clicked.connect(self._generate_structure)

        grid.addWidget(self.btn_load_smiles, 0, 0)
        grid.addWidget(self.btn_validate, 0, 1)
        grid.addWidget(self.btn_load_sdf, 1, 0)
        grid.addWidget(self.btn_load_mol, 1, 1)
        grid.addWidget(self.btn_draw, 2, 0)
        grid.addWidget(self.btn_generate, 2, 1)
        layout.addWidget(buttons_group)

        self.status_label = QLabel("Enter a molecule to begin.")
        self.status_label.setWordWrap(True)
        self.status_label.setStyleSheet("color:#8b949e;padding:4px;")
        layout.addWidget(self.status_label)
        layout.addStretch(1)

    # ── Public API ────────────────────────────────────────────────────────

    def set_smiles(self, smiles: str, emit: bool = True) -> None:
        """Programmatically set the SMILES text and optionally validate it."""
        self.smiles_edit.setPlainText(smiles)
        if emit:
            self._load_from_text()

    def current_text(self) -> str:
        return self.smiles_edit.toPlainText().strip()

    # ── Handlers ──────────────────────────────────────────────────────────

    def _on_example_selected(self, index: int) -> None:
        smiles = self.example_combo.currentData()
        if smiles:
            self.smiles_edit.setPlainText(smiles)
            self._load_from_text()

    def _load_from_text(self) -> None:
        smiles = self.current_text()
        if not smiles:
            self._set_status("Please enter a SMILES string.", ok=False)
            return
        mol = mol_utils.parse_smiles(smiles)
        if mol is None:
            self._set_status("Invalid SMILES — could not parse molecule.", ok=False)
            return
        canonical = mol_utils.canonical_smiles(mol)
        self._set_status(
            f"Valid molecule loaded ({mol.GetNumAtoms()} atoms).", ok=True
        )
        self.molecule_changed.emit(canonical)

    def _validate(self) -> None:
        smiles = self.current_text()
        mol = mol_utils.parse_smiles(smiles)
        if mol is None:
            self._set_status("✗ Invalid SMILES.", ok=False)
            return
        profile = mol_utils.compute_profile(mol)
        msg = (
            f"✓ Valid: {profile.formula}, MW {profile.mol_weight:.1f}, "
            f"QED {profile.qed:.2f}"
            if profile
            else "✓ Valid molecule."
        )
        self._set_status(msg, ok=True)
        self.molecule_changed.emit(mol_utils.canonical_smiles(mol))

    def _load_file(self, kind: str) -> None:
        filters = {
            "sdf": "SDF files (*.sdf *.sd);;All files (*)",
            "mol": "MOL files (*.mol);;All files (*)",
        }
        path, _ = QFileDialog.getOpenFileName(
            self, f"Open {kind.upper()} file", "", filters.get(kind, "All files (*)")
        )
        if not path:
            return
        try:
            mol = mol_utils.mol_from_file(path)
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Error loading file: {exc}", ok=False)
            return
        if mol is None:
            self._set_status("Could not read a valid molecule from file.", ok=False)
            return
        smiles = mol_utils.canonical_smiles(mol)
        self.smiles_edit.setPlainText(smiles)
        self._set_status(f"Loaded molecule from {path.split('/')[-1]}.", ok=True)
        self.molecule_changed.emit(smiles)

    def _generate_structure(self) -> None:
        """Clean up / canonicalise the current SMILES (and report the result)."""
        smiles = self.current_text()
        mol = mol_utils.parse_smiles(smiles)
        if mol is None:
            self._set_status("Nothing to generate — enter a valid SMILES first.", ok=False)
            return
        canonical = mol_utils.canonical_smiles(mol)
        self.smiles_edit.setPlainText(canonical)
        self._set_status("Canonical 2D structure generated.", ok=True)
        self.molecule_changed.emit(canonical)

    def _open_sketch(self) -> None:
        from ui.dialogs import SketchDialog

        dialog = SketchDialog(self.current_text(), self)
        if dialog.exec():
            smiles = dialog.smiles()
            if smiles:
                self.smiles_edit.setPlainText(smiles)
                self._load_from_text()

    def _set_status(self, message: str, ok: bool = True) -> None:
        color = "#3fb950" if ok else "#f85149"
        self.status_label.setStyleSheet(f"color:{color};padding:4px;")
        self.status_label.setText(message)
        self.status_message.emit(message)
