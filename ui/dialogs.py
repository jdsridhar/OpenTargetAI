"""
Dialogs for OpenTargetAI: a molecule builder, an import progress dialog and an
about box.
"""

from __future__ import annotations

from typing import Optional

from typing import Callable

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QProgressBar,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from similarity import mol_utils
from ui.widgets import MoleculeSvgWidget

# Common scaffolds / building blocks offered as starting points in the builder.
TEMPLATES = [
    ("— Choose a scaffold —", ""),
    ("Benzene", "c1ccccc1"),
    ("Pyridine", "c1ccncc1"),
    ("Indole", "c1ccc2[nH]ccc2c1"),
    ("Piperazine", "C1CNCCN1"),
    ("Piperidine", "C1CCNCC1"),
    ("Quinazoline", "c1ccc2ncncc2c1"),
    ("Purine", "c1nc2[nH]cnc2cn1"),
    ("Biphenyl", "c1ccc(-c2ccccc2)cc1"),
    ("Benzimidazole", "c1ccc2[nH]cnc2c1"),
    ("Sulfonamide benzene", "O=S(=O)(N)c1ccccc1"),
    ("Carboxylic acid (propanoic)", "CCC(=O)O"),
    ("Diphenylmethane", "C(c1ccccc1)c1ccccc1"),
]


class SketchDialog(QDialog):
    """A lightweight molecule builder: pick a scaffold or edit SMILES live."""

    def __init__(self, initial: str = "", parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Molecule Builder")
        self.setMinimumWidth(520)
        self._smiles = initial

        layout = QVBoxLayout(self)
        layout.addWidget(
            QLabel(
                "Pick a scaffold to start from, then refine the SMILES below. "
                "The preview updates live."
            )
        )

        self.template_combo = QComboBox()
        for label, smiles in TEMPLATES:
            self.template_combo.addItem(label, smiles)
        self.template_combo.currentIndexChanged.connect(self._on_template)
        layout.addWidget(self.template_combo)

        self.svg = MoleculeSvgWidget(460, 280)
        layout.addWidget(self.svg)

        edit_row = QHBoxLayout()
        edit_row.addWidget(QLabel("SMILES:"))
        self.smiles_edit = QLineEdit(initial)
        self.smiles_edit.textChanged.connect(self._on_text)
        edit_row.addWidget(self.smiles_edit, 1)
        layout.addLayout(edit_row)

        self.status = QLabel("")
        layout.addWidget(self.status)

        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)

        self._on_text(initial)

    def _on_template(self) -> None:
        smiles = self.template_combo.currentData()
        if smiles:
            self.smiles_edit.setText(smiles)

    def _on_text(self, text: str) -> None:
        mol = mol_utils.parse_smiles(text)
        ok = mol is not None
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(ok)
        if ok:
            self._smiles = mol_utils.canonical_smiles(mol)
            self.svg.set_mol(mol)
            self.status.setText(f"✓ Valid ({mol.GetNumAtoms()} atoms)")
            self.status.setStyleSheet("color:#3fb950;")
        else:
            self.svg.clear_molecule()
            self.status.setText("✗ Invalid SMILES" if text else "")
            self.status.setStyleSheet("color:#f85149;")

    def smiles(self) -> str:
        return self._smiles


class ProgressDialog(QDialog):
    """A modal progress dialog driven by a worker's progress signal."""

    cancelled = pyqtSignal()

    def __init__(self, title: str, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumWidth(420)
        self.setModal(True)
        layout = QVBoxLayout(self)
        self.label = QLabel("Starting…")
        self.label.setWordWrap(True)
        layout.addWidget(self.label)
        self.bar = QProgressBar()
        self.bar.setRange(0, 100)
        layout.addWidget(self.bar)
        btn_row = QHBoxLayout()
        btn_row.addStretch(1)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self._on_cancel)
        btn_row.addWidget(self.cancel_button)
        layout.addLayout(btn_row)

    def update_progress(self, current: int, total: int, message: str) -> None:
        pct = int((current / total) * 100) if total else 0
        self.bar.setValue(min(100, pct))
        if message:
            self.label.setText(message)

    def _on_cancel(self) -> None:
        self.cancelled.emit()
        self.reject()


class ExpandDialog(QDialog):
    """
    A large, resizable, non-modal window that hosts an expanded copy of a UI
    section. The hosted widget is built fresh by the caller so it always
    reflects the current application state. Esc or Close dismisses it.
    """

    def __init__(
        self, title: str, content: QWidget, parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(f"{title} — expanded")
        self.setModal(False)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.resize(960, 720)

        layout = QVBoxLayout(self)
        layout.addWidget(content, 1)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.close)
        buttons.button(QDialogButtonBox.StandardButton.Close).clicked.connect(
            self.close
        )
        layout.addWidget(buttons)

    @staticmethod
    def open_for(
        title: str, factory: Callable[[], Optional[QWidget]], parent: QWidget
    ) -> Optional["ExpandDialog"]:
        """Build content via ``factory`` and show it in a new ExpandDialog."""
        content = factory()
        if content is None:
            return None
        dialog = ExpandDialog(title, content, parent.window())
        dialog.show()
        dialog.raise_()
        return dialog


class AboutDialog(QDialog):
    """Application about box."""

    def __init__(self, version: str, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("About OpenTargetAI")
        self.setMinimumWidth(480)
        layout = QVBoxLayout(self)

        title = QLabel("OpenTargetAI")
        title.setStyleSheet("color:#58a6ff;font-size:22px;font-weight:bold;")
        layout.addWidget(title)
        layout.addWidget(
            QLabel(f"Open Source Ligand-Based Target Prediction Platform — v{version}")
        )

        body = QLabel(
            "OpenTargetAI predicts candidate protein targets for a query molecule "
            "using ligand-based target fishing over public bioactivity data "
            "(ChEMBL, BindingDB). It computes molecular-fingerprint similarity, "
            "aggregates bioactivity evidence from similar compounds, and reports "
            "transparent, explainable confidence scores.\n\n"
            "This tool is for research and educational use only. It is NOT a "
            "reproduction of SwissTargetPrediction and its predictions must not "
            "be used for clinical or regulatory decisions.\n\n"
            "Built with RDKit, PyQt6, PyQtGraph, Plotly, NetworkX, scikit-learn."
        )
        body.setWordWrap(True)
        body.setStyleSheet("color:#c9d1d9;")
        layout.addWidget(body)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.reject)
        buttons.accepted.connect(self.accept)
        buttons.button(QDialogButtonBox.StandardButton.Close).clicked.connect(
            self.accept
        )
        layout.addWidget(buttons)
