"""
Reusable UI widgets for OpenTargetAI.

Small, theme-aware building blocks shared across the panels: an SVG molecule
viewer, a key/value property list, drug-likeness badges and a stat card.
"""

from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import QByteArray, Qt
from PyQt6.QtSvgWidgets import QSvgWidget
from PyQt6.QtWidgets import (
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QSizePolicy,
    QToolButton,
    QVBoxLayout,
    QWidget,
)
from rdkit import Chem

from similarity import mol_utils


def expand_button(tooltip: str = "Expand to a larger window") -> QToolButton:
    """A small flat ⤢ button used in section headers to open an expanded view."""
    btn = QToolButton()
    btn.setText("⤢")
    btn.setToolTip(tooltip)
    btn.setCursor(Qt.CursorShape.PointingHandCursor)
    btn.setStyleSheet(
        "QToolButton{color:#8b949e;border:1px solid #30363d;border-radius:4px;"
        "padding:1px 6px;font-size:13px;font-weight:bold;background:#21262d;}"
        "QToolButton:hover{color:#58a6ff;border-color:#58a6ff;}"
    )
    return btn


class MoleculeSvgWidget(QSvgWidget):
    """A QSvgWidget that renders an RDKit molecule using the dark theme."""

    def __init__(
        self, width: int = 440, height: int = 320, parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(parent)
        self._w = width
        self._h = height
        self.setMinimumSize(240, 180)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.clear_molecule()

    def set_mol(
        self,
        mol: Optional[Chem.Mol],
        highlight_atoms: Optional[list[int]] = None,
        legend: str = "",
    ) -> None:
        """Render an RDKit molecule (optionally with highlighted atoms)."""
        w = max(self.width(), 200)
        h = max(self.height(), 160)
        svg = mol_utils.render_svg(
            mol, w, h, highlight_atoms=highlight_atoms, dark=True, legend=legend
        )
        self.load(QByteArray(svg.encode("utf-8")))

    def set_smiles(self, smiles: str, **kwargs) -> None:
        """Render a molecule from a SMILES string."""
        self.set_mol(mol_utils.parse_smiles(smiles), **kwargs)

    def clear_molecule(self) -> None:
        """Show the empty-state placeholder."""
        self.set_mol(None)


class SectionTitle(QLabel):
    """A styled section heading."""

    def __init__(self, text: str, parent: Optional[QWidget] = None) -> None:
        super().__init__(text, parent)
        self.setObjectName("sectionTitle")


class Badge(QLabel):
    """A small coloured pill used for pass/fail and confidence indicators."""

    COLORS = {
        "ok": ("#0f3d23", "#3fb950"),
        "warn": ("#3d2e0f", "#d29922"),
        "bad": ("#3d1418", "#f85149"),
        "info": ("#0d2b4d", "#58a6ff"),
        "muted": ("#21262d", "#8b949e"),
    }

    def __init__(
        self, text: str, kind: str = "info", parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(text, parent)
        self.set_kind(kind)
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)

    def set_kind(self, kind: str) -> None:
        bg, fg = self.COLORS.get(kind, self.COLORS["info"])
        self.setStyleSheet(
            f"background-color:{bg}; color:{fg}; border-radius:9px;"
            f"padding:2px 10px; font-weight:bold; font-size:11px;"
        )


class PropertyList(QWidget):
    """A label / value property display — one field per row.

    Uses a form layout so every field occupies its own row, the value wraps
    within its cell, and long values fall onto a new line under the label when
    the panel is narrow (rather than overflowing into the next column).
    """

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._form = QFormLayout(self)
        self._form.setContentsMargins(4, 4, 4, 4)
        self._form.setVerticalSpacing(6)
        self._form.setHorizontalSpacing(12)
        self._form.setLabelAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self._form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        self._form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        )

    def set_rows(self, rows: list[tuple[str, str]]) -> None:
        """Replace the displayed rows."""
        self.clear()
        for label, value in rows:
            lbl = QLabel(f"{label}")
            lbl.setStyleSheet("color:#8b949e;")
            val = QLabel(str(value))
            val.setStyleSheet("color:#c9d1d9; font-weight:bold;")
            val.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            val.setWordWrap(True)
            self._form.addRow(lbl, val)

    def clear(self) -> None:
        while self._form.count():
            item = self._form.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()


class BadgeRow(QWidget):
    """A horizontal flow of badges (e.g. drug-likeness rules)."""

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._layout = QHBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._layout.setSpacing(6)
        self._layout.addStretch(1)

    def set_badges(self, badges: list[tuple[str, str]]) -> None:
        """Set badges as a list of (text, kind) tuples."""
        while self._layout.count():
            item = self._layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        for text, kind in badges:
            self._layout.addWidget(Badge(text, kind))
        self._layout.addStretch(1)


class StatCard(QFrame):
    """A compact card showing a single headline statistic."""

    def __init__(
        self, title: str, value: str = "—", parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(parent)
        self.setObjectName("statCard")
        self.setStyleSheet(
            "#statCard{background:#161b22;border:1px solid #30363d;"
            "border-radius:8px;}"
        )
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 8, 12, 8)
        self._value = QLabel(value)
        self._value.setObjectName("statValue")
        self._value.setStyleSheet("color:#58a6ff;font-size:20px;font-weight:bold;")
        self._title = QLabel(title)
        self._title.setStyleSheet("color:#8b949e;font-size:11px;")
        layout.addWidget(self._value)
        layout.addWidget(self._title)

    def set_value(self, value: str) -> None:
        self._value.setText(str(value))
