"""
Similarity explainer for OpenTargetAI.

Makes every prediction interpretable. For a selected target it shows:

* the confidence broken down into its four evidence components;
* the supporting ligands (with similarity and bioactivity);
* a side-by-side depiction of the query and a chosen supporting ligand with the
  shared Morgan substructure highlighted on the query;
* the common Bemis–Murcko scaffold / maximum common substructure of the
  supporting evidence.
"""

from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from similarity import mol_utils
from similarity.fingerprint_engine import FingerprintEngine
from ui.dialogs import ExpandDialog
from ui.widgets import MoleculeSvgWidget, SectionTitle, expand_button


class _ComponentStat(QFrame):
    """A compact cell showing one confidence component as a percentage."""

    def __init__(self, label: str, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setStyleSheet(
            "QFrame{background:#161b22;border:1px solid #30363d;border-radius:6px;}"
        )
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 4, 10, 4)
        layout.setSpacing(0)
        self.value = QLabel("—")
        self.value.setStyleSheet(
            "color:#58a6ff;font-weight:bold;font-size:16px;border:0;"
        )
        name = QLabel(label)
        name.setStyleSheet("color:#8b949e;font-size:11px;border:0;")
        layout.addWidget(self.value)
        layout.addWidget(name)

    def set_value(self, value: Optional[float]) -> None:
        self.value.setText(f"{value * 100:.0f}%" if value is not None else "—")


class ExplainerView(QWidget):
    """Explainable-prediction panel for one target."""

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.fp_engine = FingerprintEngine()
        self._query_mol = None
        self._pred = None
        self._current_highlight: list[int] = []
        self._current_ref_mol = None
        self._current_ligand_caption = "Supporting ligand"
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)
        layout.addWidget(SectionTitle("Why was this target predicted?"))

        self.headline = QLabel("Select a predicted target to see the evidence.")
        self.headline.setWordWrap(True)
        self.headline.setStyleSheet("color:#c9d1d9;")
        layout.addWidget(self.headline)

        # Confidence breakdown — one compact row of numbers (multiple columns).
        conf_group = QGroupBox("Confidence breakdown")
        cg = QHBoxLayout(conf_group)
        cg.setContentsMargins(8, 4, 8, 4)
        cg.setSpacing(8)
        self.stat_similarity = _ComponentStat("Similarity")
        self.stat_bioactivity = _ComponentStat("Bioactivity")
        self.stat_occurrence = _ComponentStat("Occurrence")
        self.stat_consistency = _ComponentStat("Consistency")
        for s in (
            self.stat_similarity,
            self.stat_bioactivity,
            self.stat_occurrence,
            self.stat_consistency,
        ):
            cg.addWidget(s)
        conf_group.setMaximumHeight(72)
        layout.addWidget(conf_group)

        # ── Resizable area: structures (top) over supporting ligands (bottom) ──
        body = QSplitter(Qt.Orientation.Vertical)
        body.setHandleWidth(8)
        body.setChildrenCollapsible(False)

        # Structures: query | ligand in their own horizontal splitter.
        struct_group = QGroupBox("Shared substructure (query vs. selected ligand)")
        sg = QVBoxLayout(struct_group)
        sg.setContentsMargins(6, 6, 6, 6)
        struct_header = QHBoxLayout()
        struct_header.addStretch(1)
        self.expand_struct_btn = expand_button("Expand the structure comparison")
        self.expand_struct_btn.clicked.connect(self._expand_structures)
        struct_header.addWidget(self.expand_struct_btn)
        sg.addLayout(struct_header)
        struct_split = QSplitter(Qt.Orientation.Horizontal)
        struct_split.setHandleWidth(8)
        struct_split.setChildrenCollapsible(False)

        q_panel = QWidget()
        q_box = QVBoxLayout(q_panel)
        q_box.setContentsMargins(0, 0, 0, 0)
        q_title = QLabel("Query (matched atoms highlighted)")
        q_title.setStyleSheet("color:#8b949e;")
        q_box.addWidget(q_title)
        self.query_svg = MoleculeSvgWidget(280, 200)
        self.query_svg.setMinimumSize(180, 150)
        q_box.addWidget(self.query_svg, 1)

        l_panel = QWidget()
        l_box = QVBoxLayout(l_panel)
        l_box.setContentsMargins(0, 0, 0, 0)
        self.ligand_caption = QLabel("Supporting ligand")
        self.ligand_caption.setStyleSheet("color:#8b949e;")
        l_box.addWidget(self.ligand_caption)
        self.ligand_svg = MoleculeSvgWidget(280, 200)
        self.ligand_svg.setMinimumSize(180, 150)
        l_box.addWidget(self.ligand_svg, 1)

        struct_split.addWidget(q_panel)
        struct_split.addWidget(l_panel)
        struct_split.setSizes([400, 400])
        sg.addWidget(struct_split, 1)

        self.scaffold_label = QLabel("")
        self.scaffold_label.setWordWrap(True)
        self.scaffold_label.setStyleSheet(
            "color:#8b949e;font-family:Consolas,monospace;font-size:11px;"
        )
        sg.addWidget(self.scaffold_label)
        body.addWidget(struct_group)

        # Supporting ligands table
        sup_group = QGroupBox("Supporting ligands (click a row to compare)")
        sup_layout = QVBoxLayout(sup_group)
        sup_layout.setContentsMargins(6, 6, 6, 6)
        sup_header = QHBoxLayout()
        sup_header.addStretch(1)
        self.expand_ligands_btn = expand_button("Expand the supporting-ligands table")
        self.expand_ligands_btn.clicked.connect(self._expand_ligands)
        sup_header.addWidget(self.expand_ligands_btn)
        sup_layout.addLayout(sup_header)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(
            ["ChEMBL/ID", "Similarity", "Activity", "Value", "pChEMBL"]
        )
        self.table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.table.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection
        )
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self._on_row_selected)
        sup_layout.addWidget(self.table)
        body.addWidget(sup_group)

        body.setStretchFactor(0, 3)
        body.setStretchFactor(1, 2)
        body.setSizes([300, 200])
        layout.addWidget(body, 1)

    # ── Public API ────────────────────────────────────────────────────────

    def set_target(self, query_smiles: str, pred) -> None:
        """Show the evidence behind a single prediction."""
        self._query_mol = mol_utils.parse_smiles(query_smiles)
        self._pred = pred

        self.headline.setText(
            f"<b>{pred.gene_symbol or pred.target_name}</b> is predicted from "
            f"<b>{pred.n_supporting_ligands}</b> similar ligand(s) "
            f"(mean similarity {pred.mean_similarity:.2f}, "
            f"max {pred.max_similarity:.2f})."
        )

        self.stat_similarity.set_value(pred.similarity_weight)
        self.stat_bioactivity.set_value(pred.bioactivity_weight)
        self.stat_occurrence.set_value(pred.occurrence_frequency)
        self.stat_consistency.set_value(pred.activity_consistency)

        self._populate_table(pred)
        self._update_scaffold(pred)

        if pred.supporting_compounds:
            self.table.selectRow(0)
        else:
            self.query_svg.set_mol(self._query_mol)
            self.ligand_svg.clear_molecule()

    def clear(self) -> None:
        self._pred = None
        self._query_mol = None
        self._current_highlight = []
        self._current_ref_mol = None
        self._current_ligand_caption = "Supporting ligand"
        self.headline.setText("Select a predicted target to see the evidence.")
        self.table.setRowCount(0)
        self.query_svg.clear_molecule()
        self.ligand_svg.clear_molecule()
        self.scaffold_label.setText("")
        for s in (
            self.stat_similarity,
            self.stat_bioactivity,
            self.stat_occurrence,
            self.stat_consistency,
        ):
            s.set_value(None)

    # ── Internals ─────────────────────────────────────────────────────────

    def _populate_table(self, pred) -> None:
        comps = pred.supporting_compounds
        self.table.setRowCount(len(comps))
        for row, comp in enumerate(comps):
            value = comp.get("activity_value")
            value_str = (
                f"{value:g} {comp.get('units', '')}".strip()
                if value is not None
                else "—"
            )
            pchembl = comp.get("pchembl_value")
            cells = [
                comp.get("chembl_id") or f"cpd{comp['compound_id']}",
                f"{comp.get('similarity', 0):.3f}",
                comp.get("activity_type") or "—",
                value_str,
                f"{pchembl:.2f}" if pchembl else "—",
            ]
            for col, text in enumerate(cells):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                if col == 0:
                    item.setData(Qt.ItemDataRole.UserRole, comp)
                self.table.setItem(row, col, item)
        self.table.resizeColumnsToContents()

    def _update_scaffold(self, pred) -> None:
        """Show the common scaffold of the query and its supporting ligands."""
        if self._query_mol is None:
            self.scaffold_label.setText("")
            return
        scaffold = mol_utils.murcko_scaffold_smiles(self._query_mol)
        mols = [self._query_mol]
        for comp in pred.supporting_compounds[:6]:
            m = mol_utils.parse_smiles(comp.get("smiles", ""))
            if m is not None:
                mols.append(m)
        mcs = mol_utils.maximum_common_substructure(mols)
        mcs_smarts = mcs.smartsString if mcs and mcs.numAtoms > 0 else "—"
        self.scaffold_label.setText(
            f"Query Murcko scaffold:  {scaffold or '—'}\n"
            f"Maximum common substructure (query + supporters):  {mcs_smarts}"
        )

    def _on_row_selected(self) -> None:
        items = self.table.selectedItems()
        if not items or self._query_mol is None:
            return
        row = items[0].row()
        first = self.table.item(row, 0)
        comp = first.data(Qt.ItemDataRole.UserRole) if first else None
        if not comp:
            return
        ref_mol = mol_utils.parse_smiles(comp.get("smiles", ""))
        # Highlight the atoms shared (by Morgan environment) with this ligand.
        highlight = (
            self.fp_engine.shared_morgan_atoms(self._query_mol, ref_mol)
            if ref_mol is not None
            else []
        )
        self._current_highlight = highlight
        self._current_ref_mol = ref_mol
        self._current_ligand_caption = (
            f"{comp.get('chembl_id') or 'compound'} "
            f"(similarity {comp.get('similarity', 0):.2f})"
        )
        self.query_svg.set_mol(self._query_mol, highlight_atoms=highlight)
        self.ligand_svg.set_mol(ref_mol)
        self.ligand_caption.setText(f"Supporting ligand: {self._current_ligand_caption}")

    def _expand_structures(self) -> None:
        """Open the query / supporting-ligand comparison in a large window."""
        if self._query_mol is None:
            return

        def factory() -> QWidget:
            container = QWidget()
            hbox = QHBoxLayout(container)
            q_panel = QVBoxLayout()
            q_title = QLabel("Query (matched atoms highlighted)")
            q_title.setStyleSheet("color:#8b949e;")
            q_svg = MoleculeSvgWidget(440, 560)
            q_svg.set_mol(self._query_mol, highlight_atoms=self._current_highlight)
            q_panel.addWidget(q_title)
            q_panel.addWidget(q_svg, 1)

            l_panel = QVBoxLayout()
            l_title = QLabel(f"Supporting ligand: {self._current_ligand_caption}")
            l_title.setStyleSheet("color:#8b949e;")
            l_svg = MoleculeSvgWidget(440, 560)
            l_svg.set_mol(self._current_ref_mol)
            l_panel.addWidget(l_title)
            l_panel.addWidget(l_svg, 1)

            hbox.addLayout(q_panel)
            hbox.addLayout(l_panel)
            return container

        ExpandDialog.open_for("Structure comparison", factory, self)

    def _expand_ligands(self) -> None:
        """Open the full supporting-ligands table in a large window."""
        if not self._pred or not self._pred.supporting_compounds:
            return
        comps = self._pred.supporting_compounds

        def factory() -> QWidget:
            table = QTableWidget(len(comps), 6)
            table.setHorizontalHeaderLabels(
                ["ChEMBL/ID", "Similarity", "Activity", "Value", "pChEMBL", "SMILES"]
            )
            table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
            table.setAlternatingRowColors(True)
            table.verticalHeader().setVisible(False)
            for row, comp in enumerate(comps):
                value = comp.get("activity_value")
                value_str = (
                    f"{value:g} {comp.get('units', '')}".strip()
                    if value is not None
                    else "—"
                )
                pchembl = comp.get("pchembl_value")
                cells = [
                    comp.get("chembl_id") or f"cpd{comp['compound_id']}",
                    f"{comp.get('similarity', 0):.3f}",
                    comp.get("activity_type") or "—",
                    value_str,
                    f"{pchembl:.2f}" if pchembl else "—",
                    comp.get("smiles", ""),
                ]
                for col, text in enumerate(cells):
                    item = QTableWidgetItem(text)
                    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                    table.setItem(row, col, item)
            table.resizeColumnsToContents()
            table.horizontalHeader().setStretchLastSection(True)
            return table

        ExpandDialog.open_for("Supporting ligands", factory, self)
