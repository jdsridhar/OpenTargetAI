"""
Prediction results table for OpenTargetAI.

A sortable table of predicted targets. Selecting a row emits ``target_selected``
with the underlying :class:`TargetPrediction` so the details, explainer, network
and structures panels can update.
"""

from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QHeaderView,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ui.widgets import SectionTitle

_CONF_COLORS = {
    "High": "#3fb950",
    "Medium": "#d29922",
    "Low": "#8b949e",
}


class _NumericItem(QTableWidgetItem):
    """A table item that sorts by an associated numeric sort-key."""

    def __init__(self, text: str, sort_value: float) -> None:
        super().__init__(text)
        self._sort_value = sort_value
        self.setFlags(self.flags() & ~Qt.ItemFlag.ItemIsEditable)

    def __lt__(self, other: object) -> bool:
        if isinstance(other, _NumericItem):
            return self._sort_value < other._sort_value
        return super().__lt__(other)  # type: ignore[arg-type]


class ResultsView(QWidget):
    """Sortable prediction results table."""

    target_selected = pyqtSignal(object)  # TargetPrediction

    HEADERS = [
        "Rank", "Target Name", "Gene", "UniProt", "Confidence", "Score",
        "Class", "Ligands", "Evidence", "Mean Sim", "Mean pChEMBL",
    ]

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._predictions: list = []
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.addWidget(SectionTitle("Predicted Targets"))

        self.table = QTableWidget(0, len(self.HEADERS))
        self.table.setHorizontalHeaderLabels(self.HEADERS)
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.table.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection
        )
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(6, QHeaderView.ResizeMode.Stretch)
        self.table.itemSelectionChanged.connect(self._on_selection_changed)
        layout.addWidget(self.table)

    # ── Public API ────────────────────────────────────────────────────────

    def set_results(self, predictions: list) -> None:
        """Populate the table from a list of TargetPrediction objects."""
        self._predictions = predictions
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(predictions))

        for row, pred in enumerate(predictions):
            mean_pchembl = (
                f"{pred.mean_pchembl:.2f}" if pred.mean_pchembl else "N/A"
            )
            cells = [
                _NumericItem(str(row + 1), row + 1),
                self._text_item(pred.target_name),
                self._text_item(pred.gene_symbol),
                self._text_item(pred.uniprot),
                _NumericItem(pred.confidence_label, pred.confidence_score),
                _NumericItem(f"{pred.confidence_score:.3f}", pred.confidence_score),
                self._text_item(pred.target_class),
                _NumericItem(str(pred.n_supporting_ligands), pred.n_supporting_ligands),
                _NumericItem(str(pred.evidence_count), pred.evidence_count),
                _NumericItem(f"{pred.mean_similarity:.3f}", pred.mean_similarity),
                _NumericItem(
                    mean_pchembl, pred.mean_pchembl if pred.mean_pchembl else -1
                ),
            ]
            # Store the prediction object on the first cell.
            cells[0].setData(Qt.ItemDataRole.UserRole, pred)
            # Colour the confidence cell.
            color = QColor(_CONF_COLORS.get(pred.confidence_label, "#8b949e"))
            cells[4].setForeground(color)
            font = cells[4].font()
            font.setBold(True)
            cells[4].setFont(font)
            for col, item in enumerate(cells):
                self.table.setItem(row, col, item)

        self.table.setSortingEnabled(True)
        self.table.resizeColumnsToContents()
        if predictions:
            self.table.selectRow(0)

    def clear(self) -> None:
        self._predictions = []
        self.table.setRowCount(0)

    def _text_item(self, text: str) -> QTableWidgetItem:
        item = QTableWidgetItem(text or "")
        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
        return item

    def _on_selection_changed(self) -> None:
        items = self.table.selectedItems()
        if not items:
            return
        row = items[0].row()
        first = self.table.item(row, 0)
        if first is None:
            return
        pred = first.data(Qt.ItemDataRole.UserRole)
        if pred is not None:
            self.target_selected.emit(pred)
