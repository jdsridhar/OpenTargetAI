"""
Target details panel for OpenTargetAI.

Shows the full annotation of a selected predicted target (name, gene, UniProt,
class, organism, protein function, disease associations) together with its
supporting-ligand count and bioactivity statistics, and a link out to UniProt.
"""

from __future__ import annotations

import webbrowser
from typing import Optional

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from ui.widgets import Badge, PropertyList, SectionTitle, StatCard


class TargetDetailsPanel(QWidget):
    """Detailed annotation view for a single predicted target."""

    def __init__(self, db_manager, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.db = db_manager
        self._uniprot: Optional[str] = None
        self._build_ui()

    def _build_ui(self) -> None:
        outer = QVBoxLayout(self)
        outer.setContentsMargins(10, 10, 10, 10)

        header = QHBoxLayout()
        header.addWidget(SectionTitle("Target Details"))
        header.addStretch(1)
        self.confidence_badge = Badge("—", "muted")
        header.addWidget(self.confidence_badge)
        outer.addLayout(header)

        self.title_label = QLabel("Select a target to view details.")
        self.title_label.setWordWrap(True)
        self.title_label.setStyleSheet(
            "color:#58a6ff;font-size:15px;font-weight:bold;"
        )
        outer.addWidget(self.title_label)

        # Stat cards
        cards = QHBoxLayout()
        self.card_ligands = StatCard("Supporting ligands")
        self.card_evidence = StatCard("Evidence records")
        self.card_pchembl = StatCard("Mean pChEMBL")
        cards.addWidget(self.card_ligands)
        cards.addWidget(self.card_evidence)
        cards.addWidget(self.card_pchembl)
        outer.addLayout(cards)

        # Scrollable body
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        body = QWidget()
        body_layout = QVBoxLayout(body)

        ann_group = QGroupBox("Annotation")
        ann_layout = QVBoxLayout(ann_group)
        self.properties = PropertyList()
        ann_layout.addWidget(self.properties)
        body_layout.addWidget(ann_group)

        func_group = QGroupBox("Protein Function")
        func_layout = QVBoxLayout(func_group)
        self.function_label = QLabel("")
        self.function_label.setWordWrap(True)
        self.function_label.setStyleSheet("color:#c9d1d9;")
        func_layout.addWidget(self.function_label)
        body_layout.addWidget(func_group)

        disease_group = QGroupBox("Disease Associations")
        disease_layout = QVBoxLayout(disease_group)
        self.disease_label = QLabel("")
        self.disease_label.setWordWrap(True)
        self.disease_label.setStyleSheet("color:#c9d1d9;")
        disease_layout.addWidget(self.disease_label)
        body_layout.addWidget(disease_group)

        activity_group = QGroupBox("Bioactivity Statistics (pChEMBL)")
        activity_layout = QVBoxLayout(activity_group)
        self.activity_props = PropertyList()
        activity_layout.addWidget(self.activity_props)
        body_layout.addWidget(activity_group)

        body_layout.addStretch(1)
        scroll.setWidget(body)
        outer.addWidget(scroll, 1)

        # External link
        link_row = QHBoxLayout()
        self.uniprot_button = QPushButton("Open in UniProt")
        self.uniprot_button.clicked.connect(self._open_uniprot)
        self.uniprot_button.setEnabled(False)
        link_row.addWidget(self.uniprot_button)
        link_row.addStretch(1)
        outer.addLayout(link_row)

    # ── Public API ────────────────────────────────────────────────────────

    def set_target(self, pred) -> None:
        """Populate from a TargetPrediction, enriching with the DB target row."""
        self._uniprot = pred.uniprot
        self.title_label.setText(pred.target_name)

        kind = {"High": "ok", "Medium": "warn", "Low": "muted"}.get(
            pred.confidence_label, "muted"
        )
        self.confidence_badge.setText(
            f"{pred.confidence_label}  ({pred.confidence_score:.2f})"
        )
        self.confidence_badge.set_kind(kind)

        self.card_ligands.set_value(str(pred.n_supporting_ligands))
        self.card_evidence.set_value(str(pred.evidence_count))
        self.card_pchembl.set_value(
            f"{pred.mean_pchembl:.2f}" if pred.mean_pchembl else "N/A"
        )

        # Enrich from the database (description, diseases, family).
        target = self.db.get_target_by_id(pred.target_id) or {}

        self.properties.set_rows(
            [
                ("Gene symbol", pred.gene_symbol or "—"),
                ("UniProt", pred.uniprot or "—"),
                ("Organism", pred.organism or "—"),
                ("Target class", pred.target_class or "—"),
                ("Protein family", target.get("protein_family") or "—"),
                ("Target score", f"{pred.target_score:.4f}"),
                ("Enrichment", f"{pred.enrichment_score:.4f}"),
                ("Activity types", ", ".join(pred.activity_types) or "—"),
            ]
        )

        self.function_label.setText(
            target.get("description") or "No functional annotation available."
        )
        diseases = target.get("disease_associations")
        self.disease_label.setText(
            diseases if diseases else "No disease associations recorded."
        )

        stats = pred.activity_stats()
        if stats:
            self.activity_props.set_rows(
                [
                    ("Measurements", str(stats["n"])),
                    ("Mean", f"{stats['mean']:.2f}"),
                    ("Median", f"{stats['median']:.2f}"),
                    ("Min", f"{stats['min']:.2f}"),
                    ("Max", f"{stats['max']:.2f}"),
                    ("Std. dev.", f"{stats['std']:.2f}"),
                ]
            )
        else:
            self.activity_props.set_rows([("Bioactivity data", "Not available")])

        self.uniprot_button.setEnabled(
            bool(self._uniprot)
            and not self._uniprot.startswith(("NOUNIPROT", "BDB_"))
        )

    def clear(self) -> None:
        self.title_label.setText("Select a target to view details.")
        self.properties.clear()
        self.activity_props.clear()
        self.function_label.setText("")
        self.disease_label.setText("")
        self.confidence_badge.setText("—")
        self.confidence_badge.set_kind("muted")
        for card in (self.card_ligands, self.card_evidence, self.card_pchembl):
            card.set_value("—")
        self.uniprot_button.setEnabled(False)

    def _open_uniprot(self) -> None:
        if self._uniprot:
            webbrowser.open(f"https://www.uniprot.org/uniprotkb/{self._uniprot}/entry")
