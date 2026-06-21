"""
Interactive target-prediction network for OpenTargetAI.

Renders a force-directed graph centred on the query molecule, connected to the
predicted targets (edge weight = prediction confidence) and to the supporting
similar ligands (edge weight = chemical similarity). Built on PyQtGraph so zoom,
pan and node selection come for free.
"""

from __future__ import annotations

import logging
from typing import Optional

import networkx as nx
import numpy as np
import pyqtgraph as pg
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

logger = logging.getLogger(__name__)

# Colour palette aligned with the application theme.
_BG = "#0d1117"
_QUERY_COLOR = "#f0883e"      # orange
_TARGET_COLORS = {
    "High": "#3fb950",        # green
    "Medium": "#d29922",      # amber
    "Low": "#8b949e",         # grey
}
_LIGAND_COLOR = "#58a6ff"     # blue


class NetworkView(QWidget):
    """A pan/zoom-able prediction network widget."""

    node_clicked = pyqtSignal(dict)  # emits the metadata dict of the clicked node

    def __init__(
        self, parent: Optional[QWidget] = None, show_expand: bool = True
    ) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        header = QHBoxLayout()
        self._info = QLabel(
            "Run a prediction to populate the molecule–target network."
        )
        self._info.setObjectName("sectionTitle")
        header.addWidget(self._info, 1)
        if show_expand:
            self._expand_btn = QToolButton()
            self._expand_btn.setText("⤢")
            self._expand_btn.setToolTip("Open the network in a larger window")
            self._expand_btn.setCursor(Qt.CursorShape.PointingHandCursor)
            self._expand_btn.setStyleSheet(
                "QToolButton{color:#8b949e;border:1px solid #30363d;"
                "border-radius:4px;padding:1px 6px;font-size:13px;"
                "font-weight:bold;background:#21262d;}"
                "QToolButton:hover{color:#58a6ff;border-color:#58a6ff;}"
            )
            self._expand_btn.clicked.connect(self._expand)
            header.addWidget(self._expand_btn)
        layout.addLayout(header)

        pg.setConfigOptions(antialias=True)
        self.plot = pg.PlotWidget(background=_BG)
        self.plot.hideAxis("bottom")
        self.plot.hideAxis("left")
        self.plot.setAspectLocked(True)
        self.plot.setMenuEnabled(False)
        layout.addWidget(self.plot, 1)

        self._scatter: Optional[pg.ScatterPlotItem] = None
        self._labels: list[pg.TextItem] = []
        self._last_result = None

    # ── Public API ────────────────────────────────────────────────────────

    def clear(self) -> None:
        """Remove all graph items from the view."""
        self.plot.clear()
        self._labels.clear()
        self._scatter = None

    def set_prediction(
        self,
        result,
        max_targets: int = 12,
        max_ligands: int = 18,
    ) -> None:
        """Build and display the network for a :class:`PredictionResult`."""
        self.clear()
        self._last_result = result
        if result is None or not result.predictions:
            self._info.setText("No predictions to visualise.")
            return

        graph, node_meta = self._build_graph(result, max_targets, max_ligands)
        if graph.number_of_nodes() <= 1:
            self._info.setText("Not enough data to build a network.")
            return

        pos = self._layout(graph)
        self._draw_edges(graph, pos)
        self._draw_nodes(graph, pos, node_meta)

        n_t = sum(1 for d in node_meta.values() if d["type"] == "target")
        n_l = sum(1 for d in node_meta.values() if d["type"] == "ligand")
        self._info.setText(
            f"Query → {n_t} predicted targets, {n_l} supporting ligands "
            f"(scroll to zoom, drag to pan, click a node for details)."
        )
        self.plot.autoRange()

    # ── Graph construction ────────────────────────────────────────────────

    def _build_graph(
        self, result, max_targets: int, max_ligands: int
    ) -> tuple[nx.Graph, dict]:
        """Assemble a NetworkX graph and a per-node metadata dictionary."""
        graph = nx.Graph()
        meta: dict[str, dict] = {}

        query_id = "QUERY"
        graph.add_node(query_id)
        meta[query_id] = {
            "type": "query",
            "label": "Query",
            "smiles": result.query_smiles,
            "color": _QUERY_COLOR,
            "size": 28,
        }

        seen_ligands: dict[int, str] = {}
        for pred in result.predictions[:max_targets]:
            tnode = f"T{pred.target_id}"
            graph.add_node(tnode)
            graph.add_edge(
                query_id, tnode, weight=max(0.05, pred.confidence_score)
            )
            meta[tnode] = {
                "type": "target",
                "label": pred.gene_symbol or pred.target_name[:14],
                "color": _TARGET_COLORS.get(pred.confidence_label, "#8b949e"),
                "size": 14 + 14 * pred.confidence_score,
                "confidence": pred.confidence_score,
                "target_id": pred.target_id,
                "target_name": pred.target_name,
            }

            # Link the strongest supporting ligands to this target.
            for comp in pred.supporting_compounds[:4]:
                cid = comp["compound_id"]
                if len(seen_ligands) >= max_ligands and cid not in seen_ligands:
                    continue
                lnode = seen_ligands.get(cid)
                if lnode is None:
                    lnode = f"L{cid}"
                    seen_ligands[cid] = lnode
                    graph.add_node(lnode)
                    meta[lnode] = {
                        "type": "ligand",
                        "label": comp.get("chembl_id") or f"cpd{cid}",
                        "color": _LIGAND_COLOR,
                        "size": 9,
                        "smiles": comp.get("smiles", ""),
                        "compound_id": cid,
                    }
                graph.add_edge(
                    lnode, tnode, weight=max(0.05, comp.get("similarity", 0.1))
                )

        return graph, meta

    @staticmethod
    def _layout(graph: nx.Graph) -> dict:
        """Compute a deterministic force-directed layout."""
        try:
            return nx.spring_layout(
                graph, weight="weight", seed=42, k=1.4, iterations=80
            )
        except Exception:  # noqa: BLE001
            return nx.circular_layout(graph)

    # ── Drawing ───────────────────────────────────────────────────────────

    def _draw_edges(self, graph: nx.Graph, pos: dict) -> None:
        """Draw each edge as a line whose width encodes its weight."""
        for u, v, data in graph.edges(data=True):
            x = [pos[u][0], pos[v][0]]
            y = [pos[u][1], pos[v][1]]
            weight = data.get("weight", 0.2)
            color = QColor("#30363d")
            color.setAlphaF(0.35 + 0.55 * min(1.0, weight))
            pen = pg.mkPen(color=color, width=0.8 + 3.2 * min(1.0, weight))
            self.plot.addItem(pg.PlotCurveItem(x=x, y=y, pen=pen))

    def _draw_nodes(self, graph: nx.Graph, pos: dict, meta: dict) -> None:
        """Draw nodes as a clickable scatter plot with text labels."""
        spots = []
        for node in graph.nodes():
            info = meta[node]
            x, y = pos[node]
            spots.append(
                {
                    "pos": (x, y),
                    "size": info["size"],
                    "brush": pg.mkBrush(info["color"]),
                    "pen": pg.mkPen("#0d1117", width=1.5),
                    "symbol": "o",
                    "data": info,
                }
            )
            label = pg.TextItem(
                info["label"], color="#c9d1d9", anchor=(0.5, 1.6)
            )
            label.setPos(x, y)
            self.plot.addItem(label)
            self._labels.append(label)

        self._scatter = pg.ScatterPlotItem(spots=spots, pxMode=True)
        self._scatter.sigClicked.connect(self._on_clicked)
        self.plot.addItem(self._scatter)

    def _on_clicked(self, _scatter, points) -> None:
        """Relay a node click as a ``node_clicked`` signal."""
        if not len(points):
            return
        data = points[0].data()
        if isinstance(data, dict):
            self.node_clicked.emit(data)

    def _expand(self) -> None:
        """Open the current network in a large, resizable, non-modal window."""
        if self._last_result is None or not self._last_result.predictions:
            return
        dialog = QDialog(self.window())
        dialog.setWindowTitle("Target network — expanded")
        dialog.setModal(False)
        dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        dialog.resize(1000, 760)
        layout = QVBoxLayout(dialog)
        view = NetworkView(show_expand=False)
        view.set_prediction(self._last_result)
        layout.addWidget(view, 1)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(dialog.close)
        buttons.button(QDialogButtonBox.StandardButton.Close).clicked.connect(
            dialog.close
        )
        layout.addWidget(buttons)
        dialog.show()
        dialog.raise_()
