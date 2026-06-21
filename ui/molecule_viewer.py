"""
Center panel: interactive molecule viewer for OpenTargetAI.

Shows the 2D depiction of the current query molecule together with its
physicochemical descriptors and drug-likeness indicators.
"""

from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)
from rdkit import Chem

from similarity import mol_utils
from ui.dialogs import ExpandDialog
from ui.widgets import (
    BadgeRow,
    MoleculeSvgWidget,
    PropertyList,
    SectionTitle,
    expand_button,
)


class MoleculeViewer(QWidget):
    """Displays the 2D structure and computed properties of a molecule."""

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._mol: Optional[Chem.Mol] = None
        self._build_ui()

    def _build_ui(self) -> None:
        outer = QVBoxLayout(self)
        outer.setContentsMargins(10, 10, 10, 10)
        outer.addWidget(SectionTitle("Molecule Viewer"))

        # 2D structure
        struct_group = QGroupBox("2D Structure")
        sg = QVBoxLayout(struct_group)
        header = QHBoxLayout()
        header.addStretch(1)
        self.expand_btn = expand_button("Expand the 2D structure")
        self.expand_btn.clicked.connect(self._expand_structure)
        header.addWidget(self.expand_btn)
        sg.addLayout(header)
        self.svg = MoleculeSvgWidget(480, 340)
        sg.addWidget(self.svg)
        self.formula_label = QLabel("")
        self.formula_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.formula_label.setStyleSheet(
            "color:#58a6ff;font-size:15px;font-weight:bold;"
        )
        sg.addWidget(self.formula_label)
        outer.addWidget(struct_group, 3)

        # Drug-likeness badges
        dl_group = QGroupBox("Drug-likeness")
        dl = QVBoxLayout(dl_group)
        self.badges = BadgeRow()
        dl.addWidget(self.badges)
        self.alerts_label = QLabel("")
        self.alerts_label.setWordWrap(True)
        self.alerts_label.setStyleSheet("color:#d29922;font-size:11px;")
        dl.addWidget(self.alerts_label)
        outer.addWidget(dl_group)

        # Descriptor table (scrollable)
        prop_group = QGroupBox("Physicochemical Properties")
        pg = QVBoxLayout(prop_group)
        self.properties = PropertyList()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.properties)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        pg.addWidget(scroll)
        outer.addWidget(prop_group, 2)

    # ── Public API ────────────────────────────────────────────────────────

    def set_molecule(self, smiles: str) -> None:
        """Update the viewer for a new query molecule."""
        mol = mol_utils.parse_smiles(smiles)
        self._mol = mol
        if mol is None:
            self.clear()
            return

        self.svg.set_mol(mol)
        profile = mol_utils.compute_profile(mol)
        if profile is None:
            return

        self.formula_label.setText(profile.formula)
        self.properties.set_rows(profile.as_rows())

        badges: list[tuple[str, str]] = []
        for name, (passed, _detail) in profile.rules.items():
            badges.append((name, "ok" if passed else "bad"))
        badges.append((f"QED {profile.qed:.2f}", "info"))
        self.badges.set_badges(badges)

        if profile.alerts:
            self.alerts_label.setText(
                "⚠ Structural alerts (PAINS): " + ", ".join(profile.alerts[:4])
            )
        else:
            self.alerts_label.setText("No PAINS structural alerts detected.")
            self.alerts_label.setStyleSheet("color:#3fb950;font-size:11px;")

    def clear(self) -> None:
        """Reset to the empty state."""
        self._mol = None
        self.svg.clear_molecule()
        self.formula_label.setText("")
        self.properties.clear()
        self.badges.set_badges([])
        self.alerts_label.setText("")

    def _expand_structure(self) -> None:
        """Open the current 2D structure in a large resizable window."""
        if self._mol is None:
            return

        def factory() -> QWidget:
            big = MoleculeSvgWidget(880, 640)
            big.set_mol(self._mol)
            return big

        ExpandDialog.open_for("2D Structure", factory, self)
