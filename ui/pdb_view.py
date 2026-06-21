"""
Structures & Docking tab for OpenTargetAI.

For the selected predicted target it retrieves experimental 3D structures
(PDBe / RCSB), lets the user download a coordinate file, and assembles a
ready-to-dock AutoDock Vina package (protein + 3D ligand + grid box + config).
All network and file IO runs in background workers to keep the UI responsive.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from integration.docking_prep import DockingPreparation
from integration.pdb_client import PDBClient, PDBLookupError
from ui.widgets import SectionTitle
from workers.workers import GenericWorker


class StructuresView(QWidget):
    """PDB structure browser and docking-package builder."""

    status_message = pyqtSignal(str)

    def __init__(self, project_root: str, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.project_root = project_root
        self._pred = None
        self._query_smiles = ""
        self._structures: list = []
        self._worker: Optional[GenericWorker] = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.addWidget(SectionTitle("Structures & Docking Preparation"))

        self.target_label = QLabel("Select a predicted target to fetch structures.")
        self.target_label.setStyleSheet("color:#58a6ff;font-weight:bold;")
        layout.addWidget(self.target_label)

        # Action buttons
        actions = QHBoxLayout()
        self.fetch_button = QPushButton("Retrieve PDB Structures")
        self.fetch_button.clicked.connect(self._fetch_structures)
        self.fetch_button.setEnabled(False)
        self.download_button = QPushButton("Download PDB")
        self.download_button.clicked.connect(self._download_selected)
        self.download_button.setEnabled(False)
        self.dock_button = QPushButton("Prepare for Docking")
        self.dock_button.clicked.connect(self._prepare_docking)
        self.dock_button.setEnabled(False)
        actions.addWidget(self.fetch_button)
        actions.addWidget(self.download_button)
        actions.addWidget(self.dock_button)
        actions.addStretch(1)
        layout.addLayout(actions)

        # Structures table
        struct_group = QGroupBox("Available Structures")
        sg = QVBoxLayout(struct_group)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(
            ["PDB ID", "Resolution", "Method", "Chain", "Coverage"]
        )
        self.table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.table.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection
        )
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self._on_structure_selected)
        sg.addWidget(self.table)
        layout.addWidget(struct_group, 1)

        # Log
        log_group = QGroupBox("Activity Log")
        lg = QVBoxLayout(log_group)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(140)
        lg.addWidget(self.log)
        layout.addWidget(log_group)

    # ── Public API ────────────────────────────────────────────────────────

    def set_target(self, query_smiles: str, pred) -> None:
        """Set the active target / query for structure and docking actions."""
        self._pred = pred
        self._query_smiles = query_smiles
        self._structures = []
        self.table.setRowCount(0)
        self.download_button.setEnabled(False)
        has_uniprot = bool(pred.uniprot) and not pred.uniprot.startswith(
            ("NOUNIPROT", "BDB_")
        )
        self.target_label.setText(
            f"{pred.target_name}  ({pred.uniprot or 'no UniProt'})"
        )
        self.fetch_button.setEnabled(has_uniprot)
        self.dock_button.setEnabled(bool(query_smiles))
        if not has_uniprot:
            self._append_log(
                "This target has no resolvable UniProt accession; structure "
                "retrieval is unavailable. Docking prep can still build a "
                "ligand-only package."
            )

    def clear(self) -> None:
        self._pred = None
        self._query_smiles = ""
        self.table.setRowCount(0)
        self.target_label.setText("Select a predicted target to fetch structures.")
        for b in (self.fetch_button, self.download_button, self.dock_button):
            b.setEnabled(False)

    # ── Structure retrieval ───────────────────────────────────────────────

    def _fetch_structures(self) -> None:
        if not self._pred:
            return
        uniprot = self._pred.uniprot
        self._append_log(f"Retrieving structures for {uniprot} …")
        self.fetch_button.setEnabled(False)

        def task():
            return PDBClient().get_structures_for_uniprot(uniprot, max_results=30)

        self._run_worker(task, self._on_structures_ready, self._on_worker_error)

    def _on_structures_ready(self, structures: list) -> None:
        self.fetch_button.setEnabled(True)
        self._structures = structures or []
        self.table.setRowCount(len(self._structures))
        for row, s in enumerate(self._structures):
            cov = f"{s.coverage:.0%}" if s.coverage is not None else "—"
            cells = [
                s.pdb_id,
                s.resolution_str,
                s.experimental_method,
                s.chain_id or "—",
                cov,
            ]
            for col, text in enumerate(cells):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(row, col, item)
        self.table.resizeColumnsToContents()
        if self._structures:
            self._append_log(f"Found {len(self._structures)} structure(s).")
            self.table.selectRow(0)
        else:
            self._append_log("No experimental structures found for this target.")

    def _on_structure_selected(self) -> None:
        self.download_button.setEnabled(bool(self.table.selectedItems()))

    def _selected_pdb_id(self) -> Optional[str]:
        items = self.table.selectedItems()
        if not items:
            return None
        return self.table.item(items[0].row(), 0).text()

    def _download_selected(self) -> None:
        pdb_id = self._selected_pdb_id()
        if not pdb_id:
            return
        dest = os.path.join(self.project_root, "data", "structures")
        self._append_log(f"Downloading {pdb_id} …")
        self.download_button.setEnabled(False)

        def task():
            return PDBClient().download_structure(pdb_id, dest, "pdb")

        def done(path: str):
            self.download_button.setEnabled(True)
            self._append_log(f"Downloaded to {path}")
            self.status_message.emit(f"Downloaded {pdb_id}")

        self._run_worker(task, done, self._on_worker_error)

    # ── Docking preparation ───────────────────────────────────────────────

    def _prepare_docking(self) -> None:
        if not self._pred or not self._query_smiles:
            return
        base = QFileDialog.getExistingDirectory(
            self, "Choose output directory for docking package",
            os.path.join(self.project_root, "data"),
        )
        if not base:
            return
        pred = self._pred
        gene = pred.gene_symbol or f"target{pred.target_id}"
        out_dir = os.path.join(base, f"docking_{gene}")

        # Use a previously downloaded receptor if one matches the selection.
        receptor = None
        pdb_id = self._selected_pdb_id()
        if pdb_id:
            candidate = Path(self.project_root) / "data" / "structures" / f"{pdb_id.lower()}.pdb"
            if candidate.exists():
                receptor = str(candidate)

        self._append_log(f"Preparing docking package in {out_dir} …")
        self.dock_button.setEnabled(False)
        smiles = self._query_smiles

        def task():
            dp = DockingPreparation()
            pkg = dp.prepare(
                out_dir,
                ligand_smiles=smiles,
                target_name=pred.target_name,
                uniprot=pred.uniprot,
                receptor_pdb_path=receptor,
            )
            dp.export_zip(pkg)
            return pkg

        def done(pkg):
            self.dock_button.setEnabled(True)
            self._append_log(
                f"Docking package ready: {pkg.output_dir}\n"
                f"  Grid centre {tuple(round(c, 1) for c in pkg.grid_box.center)}, "
                f"size {tuple(round(s, 1) for s in pkg.grid_box.size)}\n"
                f"  {pkg.grid_box.source}\n"
                f"  ZIP: {pkg.zip_file}"
            )
            self.status_message.emit("Docking package prepared")

        self._run_worker(task, done, self._on_worker_error)

    # ── Worker plumbing ───────────────────────────────────────────────────

    def _run_worker(self, task, on_done, on_error) -> None:
        worker = GenericWorker(task)
        worker.finished.connect(on_done)
        worker.error.connect(on_error)
        worker.finished.connect(lambda *_: self._clear_worker())
        worker.error.connect(lambda *_: self._clear_worker())
        self._worker = worker
        worker.start()

    def _clear_worker(self) -> None:
        self._worker = None

    def _on_worker_error(self, message: str) -> None:
        for b in (self.fetch_button, self.download_button, self.dock_button):
            b.setEnabled(True)
        # Re-disable fetch/download if there is no valid target/selection.
        self._append_log(f"Error: {message}")
        self.status_message.emit("Operation failed")

    def _append_log(self, text: str) -> None:
        self.log.appendPlainText(text)
