"""
Advanced visualisation tab for OpenTargetAI.

Assembles the in-application distribution charts (similarity, confidence,
activity and target-class) and offers an "Open Interactive Dashboard" action
that exports a standalone Plotly HTML report and opens it in the browser.
"""

from __future__ import annotations

import os
import tempfile
import webbrowser
from typing import Optional

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from visualization import charts
from ui.widgets import SectionTitle


class ChartsView(QWidget):
    """Distribution charts and interactive export for a prediction result."""

    status_message = pyqtSignal(str)

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._result = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)

        header = QHBoxLayout()
        header.addWidget(SectionTitle("Advanced Visualisation"))
        header.addStretch(1)
        self.plotly_button = QPushButton("Open Interactive Dashboard ↗")
        self.plotly_button.clicked.connect(self._open_dashboard)
        self.plotly_button.setEnabled(False)
        header.addWidget(self.plotly_button)
        layout.addLayout(header)

        grid = QGridLayout()

        self.sim_hist = charts.HistogramWidget("Neighbour Similarity Distribution")
        self.conf_hist = charts.HistogramWidget("Target Confidence Distribution")
        self.act_hist = charts.HistogramWidget("Supporting Activity (pChEMBL)")
        self.class_pie = charts.PieChartWidget("Target Class Distribution")

        for w, (r, c) in (
            (self._wrap(self.sim_hist), (0, 0)),
            (self._wrap(self.conf_hist), (0, 1)),
            (self._wrap(self.act_hist), (1, 0)),
            (self._wrap(self.class_pie), (1, 1)),
        ):
            grid.addWidget(w, r, c)
        layout.addLayout(grid, 1)

    @staticmethod
    def _wrap(widget: QWidget) -> QWidget:
        box = QGroupBox()
        inner = QVBoxLayout(box)
        inner.setContentsMargins(4, 4, 4, 4)
        inner.addWidget(widget)
        return box

    # ── Public API ────────────────────────────────────────────────────────

    def set_result(self, result) -> None:
        """Refresh all charts from a prediction result."""
        self._result = result
        self.plotly_button.setEnabled(result is not None and bool(result.predictions))
        if result is None:
            return
        self.sim_hist.set_values(
            charts.similarity_values(result),
            bins=15,
            color="#3fb950",
            x_label="Tanimoto / similarity",
            x_range=(0.0, 1.0),
        )
        self.conf_hist.set_values(
            charts.confidence_values(result),
            bins=12,
            color="#58a6ff",
            x_label="Confidence score",
            x_range=(0.0, 1.0),
        )
        self.act_hist.set_values(
            charts.activity_values(result),
            bins=14,
            color="#d29922",
            x_label="pChEMBL",
        )
        self.class_pie.set_data(charts.target_class_counts(result))

    def clear(self) -> None:
        self._result = None
        self.plotly_button.setEnabled(False)
        for hist in (self.sim_hist, self.conf_hist, self.act_hist):
            hist.set_values([])
        self.class_pie.set_data({})

    def _open_dashboard(self) -> None:
        if self._result is None:
            return
        try:
            path = os.path.join(
                tempfile.gettempdir(), "opentargetai_dashboard.html"
            )
            charts.export_plotly_html(self._result, path)
            webbrowser.open(f"file://{os.path.abspath(path)}")
            self.status_message.emit(f"Interactive dashboard opened: {path}")
        except Exception as exc:  # noqa: BLE001
            self.status_message.emit(f"Could not open dashboard: {exc}")
