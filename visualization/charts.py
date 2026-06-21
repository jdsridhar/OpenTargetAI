"""
Advanced visualisation widgets and exports for OpenTargetAI.

Provides lightweight, dependency-free PyQtGraph chart widgets for in-application
display (similarity / confidence / activity distributions and a target-class pie
chart) plus a Plotly exporter that writes a standalone interactive HTML report
(opened in the user's browser, since the bundled Qt build has no WebEngine).
"""

from __future__ import annotations

import logging
from collections import Counter
from typing import Optional

import numpy as np
import pyqtgraph as pg
from PyQt6.QtCore import QRectF, Qt
from PyQt6.QtGui import QColor, QFont, QPainter, QPen
from PyQt6.QtWidgets import QVBoxLayout, QWidget

logger = logging.getLogger(__name__)

_BG = "#0d1117"
_PIE_PALETTE = [
    "#58a6ff", "#3fb950", "#d29922", "#f0883e", "#bc8cff",
    "#f85149", "#39c5cf", "#db61a2", "#a371f7", "#56d364",
    "#e3b341", "#ff7b72",
]


# ── Data extraction helpers ────────────────────────────────────────────────


def similarity_values(result) -> list[float]:
    """All neighbour similarity scores from a prediction's search result."""
    sr = getattr(result, "search_result", None)
    if not sr:
        return []
    return [c.similarity for c in sr.similar_compounds]


def confidence_values(result) -> list[float]:
    """Confidence scores for every predicted target."""
    return [p.confidence_score for p in result.predictions]


def activity_values(result) -> list[float]:
    """Pooled supporting pChEMBL values across all predictions."""
    vals: list[float] = []
    for pred in result.predictions:
        for comp in pred.supporting_compounds:
            v = comp.get("pchembl_value")
            if v:
                vals.append(v)
    return vals


def target_class_counts(result) -> dict[str, int]:
    """Count predicted targets by target class."""
    counter: Counter[str] = Counter()
    for pred in result.predictions:
        counter[pred.target_class or "Unclassified"] += 1
    return dict(counter)


# ── PyQtGraph histogram ────────────────────────────────────────────────────


class HistogramWidget(QWidget):
    """A themed histogram backed by a PyQtGraph bar chart."""

    def __init__(self, title: str = "", parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.plot = pg.PlotWidget(background=_BG)
        self.plot.setMenuEnabled(False)
        self.plot.showGrid(x=False, y=True, alpha=0.2)
        self._title = title
        if title:
            self.plot.setTitle(title, color="#c9d1d9", size="10pt")
        layout.addWidget(self.plot)

    def set_values(
        self,
        values: list[float],
        bins: int = 12,
        color: str = "#58a6ff",
        x_label: str = "",
        y_label: str = "Count",
        x_range: Optional[tuple[float, float]] = None,
    ) -> None:
        """Render a histogram of ``values``."""
        self.plot.clear()
        if not values:
            self._show_placeholder()
            return
        arr = np.asarray(values, dtype=float)
        rng = x_range or (float(arr.min()), float(arr.max() + 1e-9))
        counts, edges = np.histogram(arr, bins=bins, range=rng)
        width = (edges[1] - edges[0]) * 0.9
        bar = pg.BarGraphItem(
            x=edges[:-1] + width / 2,
            height=counts,
            width=width,
            brush=pg.mkBrush(color),
            pen=pg.mkPen("#0d1117"),
        )
        self.plot.addItem(bar)
        self.plot.setLabel("bottom", x_label, color="#8b949e")
        self.plot.setLabel("left", y_label, color="#8b949e")

    def _show_placeholder(self) -> None:
        text = pg.TextItem("No data available", color="#484f58", anchor=(0.5, 0.5))
        self.plot.addItem(text)
        text.setPos(0.5, 0.5)


# ── Custom pie chart ───────────────────────────────────────────────────────


class PieChartWidget(QWidget):
    """A QPainter-based pie chart with a legend (target-class distribution)."""

    def __init__(self, title: str = "", parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._data: dict[str, int] = {}
        self._title = title
        self.setMinimumHeight(240)

    def set_data(self, data: dict[str, int]) -> None:
        """Set the {label: count} mapping and repaint."""
        self._data = {k: v for k, v in data.items() if v > 0}
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt override
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor(_BG))

        if self._title:
            painter.setPen(QColor("#c9d1d9"))
            painter.setFont(QFont("Segoe UI", 10, QFont.Weight.Bold))
            painter.drawText(
                QRectF(0, 4, self.width(), 22),
                Qt.AlignmentFlag.AlignCenter,
                self._title,
            )

        if not self._data:
            painter.setPen(QColor("#484f58"))
            painter.setFont(QFont("Segoe UI", 10))
            painter.drawText(
                self.rect(), Qt.AlignmentFlag.AlignCenter, "No data available"
            )
            painter.end()
            return

        total = sum(self._data.values())
        top = 28
        diameter = min(self.width() * 0.5, self.height() - top - 16)
        pie_rect = QRectF(16, top, diameter, diameter)

        start_angle = 90 * 16  # start at top, Qt angles are in 1/16th degrees
        legend_x = pie_rect.right() + 24
        legend_y = top + 4
        painter.setFont(QFont("Segoe UI", 9))

        for i, (label, count) in enumerate(
            sorted(self._data.items(), key=lambda kv: kv[1], reverse=True)
        ):
            span = -int(360 * 16 * count / total)
            color = QColor(_PIE_PALETTE[i % len(_PIE_PALETTE)])
            painter.setBrush(color)
            painter.setPen(QPen(QColor(_BG), 2))
            painter.drawPie(pie_rect, start_angle, span)
            start_angle += span

            # Legend entry
            painter.setBrush(color)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.drawRect(QRectF(legend_x, legend_y, 12, 12))
            painter.setPen(QColor("#c9d1d9"))
            pct = 100.0 * count / total
            painter.drawText(
                QRectF(legend_x + 18, legend_y - 2, self.width() - legend_x - 24, 16),
                Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                f"{label}  ({count}, {pct:.0f}%)",
            )
            legend_y += 20

        painter.end()


# ── Plotly interactive export ──────────────────────────────────────────────


def export_plotly_html(result, output_path: str) -> str:
    """
    Write a standalone interactive HTML dashboard (Plotly) summarising the
    prediction. Returns the output path. Raises ``RuntimeError`` if Plotly is
    unavailable.
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("plotly is required for interactive HTML export") from exc

    preds = result.predictions
    sims = similarity_values(result)
    confs = confidence_values(result)
    acts = activity_values(result)
    classes = target_class_counts(result)

    fig = make_subplots(
        rows=2,
        cols=2,
        specs=[
            [{"type": "bar"}, {"type": "domain"}],
            [{"type": "bar"}, {"type": "histogram"}],
        ],
        subplot_titles=(
            "Top predicted targets (confidence)",
            "Target class distribution",
            "Neighbour similarity distribution",
            "Supporting activity (pChEMBL)",
        ),
    )

    top = preds[:15]
    if top:
        fig.add_trace(
            go.Bar(
                x=[p.confidence_score for p in top][::-1],
                y=[p.gene_symbol or p.target_name[:18] for p in top][::-1],
                orientation="h",
                marker_color="#58a6ff",
                name="Confidence",
            ),
            row=1,
            col=1,
        )
    if classes:
        fig.add_trace(
            go.Pie(
                labels=list(classes.keys()),
                values=list(classes.values()),
                hole=0.35,
            ),
            row=1,
            col=2,
        )
    if sims:
        fig.add_trace(
            go.Histogram(x=sims, marker_color="#3fb950", nbinsx=15, name="Similarity"),
            row=2,
            col=1,
        )
    if acts:
        fig.add_trace(
            go.Histogram(x=acts, marker_color="#d29922", nbinsx=15, name="pChEMBL"),
            row=2,
            col=2,
        )

    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="#0d1117",
        plot_bgcolor="#0d1117",
        showlegend=False,
        title_text=f"OpenTargetAI — Interactive Prediction Dashboard<br>"
        f"<sub>Query: {result.query_smiles}</sub>",
        height=820,
    )
    fig.write_html(output_path, include_plotlyjs="cdn", full_html=True)
    logger.info("Interactive Plotly dashboard written to %s", output_path)
    return output_path
