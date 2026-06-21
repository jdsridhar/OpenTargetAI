"""
Main application window for OpenTargetAI.

Assembles the prediction workspace (input → viewer → controls), the results,
explainer, network, charts, structures/docking and validation tabs, and wires
their signals to the search / prediction / import / export back-ends. Long
operations run in background QThread workers so the UI stays responsive.
"""

from __future__ import annotations

import logging
import os
from typing import Optional

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction
from PyQt6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QSplitter,
    QStatusBar,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from database.db_manager import DatabaseManager
from prediction.prediction_engine import TargetPredictionEngine
from similarity.search_engine import SimilaritySearchEngine
from ui.charts_view import ChartsView
from ui.dialogs import AboutDialog, ProgressDialog
from ui.explainer_view import ExplainerView
from ui.molecule_input import MoleculeInputPanel
from ui.molecule_viewer import MoleculeViewer
from ui.pdb_view import StructuresView
from ui.prediction_controls import PredictionControls
from ui.results_view import ResultsView
from ui.target_details import TargetDetailsPanel
from ui.validation_view import ValidationView
from visualization.network_view import NetworkView
from workers.workers import ImportWorker, PredictionWorker

logger = logging.getLogger(__name__)

APP_VERSION = "1.0.0"


class MainWindow(QMainWindow):
    """The top-level OpenTargetAI window."""

    def __init__(self, db_manager: DatabaseManager) -> None:
        super().__init__()
        self.db = db_manager
        self.project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        # Shared back-end engines (the search engine caches fingerprints).
        self.search_engine = SimilaritySearchEngine(db_manager)
        self.prediction_engine = TargetPredictionEngine()

        # Runtime state
        self._query_smiles: str = ""
        self._result = None
        self._prediction_worker: Optional[PredictionWorker] = None
        self._import_worker: Optional[ImportWorker] = None
        self._progress_dialog: Optional[ProgressDialog] = None

        self.setWindowTitle("OpenTargetAI — Ligand-Based Target Prediction")
        self.resize(1480, 920)

        self._build_menu()
        self._build_central()
        self._build_status_bar()
        self._connect_signals()
        self._refresh_stats()
        # NOTE: the "load demo data?" prompt is intentionally NOT shown here.
        # It must run after the window is visible and the splash screen has
        # closed (see app.py / offer_demo_if_empty), otherwise the modal dialog
        # opens behind the always-on-top splash and the app looks frozen.

    # ── Construction ──────────────────────────────────────────────────────

    def _build_central(self) -> None:
        self.tabs = QTabWidget()
        self.predict_tab = self._build_predict_tab()
        self.results_tab = self._build_results_tab()
        self.tabs.addTab(self.predict_tab, "Predict")
        self.tabs.addTab(self.results_tab, "Results")

        self.network_view = NetworkView()
        self.tabs.addTab(self.network_view, "Network")

        self.charts_view = ChartsView()
        self.tabs.addTab(self.charts_view, "Charts")

        self.structures_view = StructuresView(self.project_root)
        self.tabs.addTab(self.structures_view, "Structures & Docking")

        self.validation_view = ValidationView(self.db, self.project_root)
        self.tabs.addTab(self.validation_view, "Validation")

        self.setCentralWidget(self.tabs)

    def _build_predict_tab(self) -> QWidget:
        self.input_panel = MoleculeInputPanel()
        self.viewer = MoleculeViewer()
        self.controls = PredictionControls()

        # Keep each panel usable while still draggable.
        self.input_panel.setMinimumWidth(260)
        self.viewer.setMinimumWidth(360)
        self.controls.setMinimumWidth(260)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self.input_panel)
        splitter.addWidget(self.viewer)
        splitter.addWidget(self.controls)
        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 4)
        splitter.setStretchFactor(2, 2)
        self._tune_splitter(splitter, [340, 700, 340])

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(splitter)
        return container

    def _build_results_tab(self) -> QWidget:
        self.results_view = ResultsView()
        self.target_details = TargetDetailsPanel(self.db)
        self.explainer = ExplainerView()

        self.target_details.setMinimumWidth(320)
        self.explainer.setMinimumWidth(380)
        self.results_view.setMinimumHeight(160)

        bottom = QSplitter(Qt.Orientation.Horizontal)
        bottom.addWidget(self.target_details)
        bottom.addWidget(self.explainer)
        bottom.setStretchFactor(0, 2)
        bottom.setStretchFactor(1, 3)
        self._tune_splitter(bottom, [460, 720])

        vertical = QSplitter(Qt.Orientation.Vertical)
        vertical.addWidget(self.results_view)
        vertical.addWidget(bottom)
        vertical.setStretchFactor(0, 2)
        vertical.setStretchFactor(1, 3)
        self._tune_splitter(vertical, [320, 560])
        return vertical

    @staticmethod
    def _tune_splitter(splitter: QSplitter, sizes: list[int]) -> None:
        """Make a splitter easy to drag and prevent panels collapsing to zero."""
        splitter.setHandleWidth(8)
        splitter.setChildrenCollapsible(False)
        splitter.setOpaqueResize(True)
        splitter.setSizes(sizes)

    def _build_menu(self) -> None:
        menubar = self.menuBar()

        file_menu = menubar.addMenu("&File")
        self.export_actions = {}
        for fmt, label in (
            ("csv", "Export Results as CSV…"),
            ("excel", "Export Results as Excel…"),
            ("html", "Export Results as HTML…"),
            ("pdf", "Export Results as PDF…"),
        ):
            action = QAction(label, self)
            action.triggered.connect(lambda _checked, f=fmt: self._export(f))
            action.setEnabled(False)
            file_menu.addAction(action)
            self.export_actions[fmt] = action
        file_menu.addSeparator()
        exit_action = QAction("E&xit", self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        db_menu = menubar.addMenu("&Database")
        demo_action = QAction("Load Demo Dataset", self)
        demo_action.triggered.connect(lambda: self._start_import("demo"))
        db_menu.addAction(demo_action)
        chembl_action = QAction("Import ChEMBL TSV…", self)
        chembl_action.triggered.connect(lambda: self._import_file("chembl"))
        db_menu.addAction(chembl_action)
        bindingdb_action = QAction("Import BindingDB TSV…", self)
        bindingdb_action.triggered.connect(lambda: self._import_file("bindingdb"))
        db_menu.addAction(bindingdb_action)
        db_menu.addSeparator()
        stats_action = QAction("Database Statistics", self)
        stats_action.triggered.connect(self._show_stats)
        db_menu.addAction(stats_action)
        clear_cache_action = QAction("Clear Prediction Cache", self)
        clear_cache_action.triggered.connect(self._clear_cache)
        db_menu.addAction(clear_cache_action)

        run_menu = menubar.addMenu("&Run")
        run_action = QAction("Run Prediction", self)
        run_action.setShortcut("F5")
        run_action.triggered.connect(self._run_prediction)
        run_menu.addAction(run_action)
        validate_action = QAction("Open Validation", self)
        validate_action.triggered.connect(
            lambda: self.tabs.setCurrentWidget(self.validation_view)
        )
        run_menu.addAction(validate_action)

        help_menu = menubar.addMenu("&Help")
        about_action = QAction("About OpenTargetAI", self)
        about_action.triggered.connect(self._show_about)
        help_menu.addAction(about_action)

    def _build_status_bar(self) -> None:
        self.setStatusBar(QStatusBar())
        self.stats_label = QLabel("")
        self.statusBar().addPermanentWidget(self.stats_label)
        self.statusBar().showMessage("Ready.")

    def _connect_signals(self) -> None:
        self.input_panel.molecule_changed.connect(self._on_molecule_changed)
        self.input_panel.status_message.connect(self.statusBar().showMessage)
        self.controls.run_requested.connect(lambda _params: self._run_prediction())
        self.controls.cancel_requested.connect(self._cancel_prediction)
        self.results_view.target_selected.connect(self._on_target_selected)
        self.network_view.node_clicked.connect(self._on_node_clicked)
        self.charts_view.status_message.connect(self.statusBar().showMessage)
        self.structures_view.status_message.connect(self.statusBar().showMessage)
        self.validation_view.status_message.connect(self.statusBar().showMessage)

    # ── Molecule handling ─────────────────────────────────────────────────

    def _on_molecule_changed(self, smiles: str) -> None:
        self._query_smiles = smiles
        self.viewer.set_molecule(smiles)
        self.statusBar().showMessage("Molecule loaded — ready to predict.")

    # ── Prediction ────────────────────────────────────────────────────────

    def _run_prediction(self) -> None:
        if self._prediction_worker is not None:
            return
        if not self._query_smiles:
            self.statusBar().showMessage("Load a valid molecule first.")
            QMessageBox.information(
                self, "No molecule", "Please load a valid query molecule first."
            )
            return

        stats = self.db.get_database_stats()
        if stats["compounds"] == 0:
            if self._confirm_demo():
                self._start_import("demo")
            return

        params = self.controls.parameters()
        self.controls.set_running(True)
        self.statusBar().showMessage("Running prediction…")

        self._prediction_worker = PredictionWorker(
            self.search_engine,
            self.prediction_engine,
            self._query_smiles,
            params["fp_type"],
            params["metric"],
            params["threshold"],
            params["top_n"],
            params["organism"],
            min_ligands=params["min_ligands"],
        )
        self._prediction_worker.progress.connect(self.controls.update_progress)
        self._prediction_worker.finished.connect(self._on_prediction_finished)
        self._prediction_worker.error.connect(self._on_prediction_error)
        self._prediction_worker.start()

    def _cancel_prediction(self) -> None:
        if self._prediction_worker is not None:
            self._prediction_worker.cancel()
            self.statusBar().showMessage("Cancelling…")

    def _on_prediction_finished(self, result) -> None:
        self.controls.set_running(False)
        self._prediction_worker = None
        self._result = result

        self.results_view.set_results(result.predictions)
        self.network_view.set_prediction(result)
        self.charts_view.set_result(result)

        for action in self.export_actions.values():
            action.setEnabled(bool(result.predictions))

        if not result.predictions:
            self.target_details.clear()
            self.explainer.clear()
            self.structures_view.clear()
            self.statusBar().showMessage(
                "No targets predicted — try lowering the threshold or loading "
                "more data."
            )
            QMessageBox.information(
                self,
                "No predictions",
                "No targets were predicted for this molecule with the current "
                "settings.\n\nTry lowering the similarity threshold, increasing "
                "Top-N, or importing more bioactivity data.",
            )
            return

        n = len(result.predictions)
        self.statusBar().showMessage(
            f"Prediction complete: {n} target(s) from "
            f"{result.n_similar_compounds} similar compound(s)."
        )
        self.tabs.setCurrentWidget(self.results_tab)
        self.results_view.table.setFocus()

    def _on_prediction_error(self, message: str) -> None:
        self.controls.set_running(False)
        self._prediction_worker = None
        logger.error("Prediction error: %s", message)
        self.statusBar().showMessage("Prediction failed.")
        QMessageBox.critical(
            self, "Prediction error", f"The prediction failed:\n\n{message}"
        )

    # ── Target selection ──────────────────────────────────────────────────

    def _on_target_selected(self, pred) -> None:
        self.target_details.set_target(pred)
        self.explainer.set_target(self._query_smiles, pred)
        self.structures_view.set_target(self._query_smiles, pred)

    def _on_node_clicked(self, node: dict) -> None:
        if node.get("type") == "target" and self._result:
            tid = node.get("target_id")
            for pred in self._result.predictions:
                if pred.target_id == tid:
                    self._on_target_selected(pred)
                    self.tabs.setCurrentWidget(self.results_tab)
                    break
        elif node.get("type") == "query":
            self.statusBar().showMessage(f"Query molecule: {node.get('smiles', '')}")

    # ── Import ────────────────────────────────────────────────────────────

    def _import_file(self, source: str) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            f"Open {source.upper()} export",
            "",
            "Tabular data (*.tsv *.csv *.txt *.gz);;All files (*)",
        )
        if path:
            self._start_import(source, path)

    def _start_import(self, source: str, file_path: Optional[str] = None) -> None:
        if self._import_worker is not None:
            QMessageBox.information(
                self, "Import in progress", "An import is already running."
            )
            return
        from database.importer import DataImporter

        importer = DataImporter(self.db)
        self._import_worker = ImportWorker(importer, source, file_path)

        title = {
            "demo": "Loading demo dataset",
            "chembl": "Importing ChEMBL data",
            "bindingdb": "Importing BindingDB data",
        }.get(source, "Importing data")
        self._progress_dialog = ProgressDialog(title, self)
        self._import_worker.progress.connect(self._progress_dialog.update_progress)
        self._import_worker.finished.connect(self._on_import_finished)
        self._import_worker.error.connect(self._on_import_error)
        self._progress_dialog.cancelled.connect(self._import_worker.requestInterruption)

        self._import_worker.start()
        self._progress_dialog.show()

    def _on_import_finished(self, count: int) -> None:
        if self._progress_dialog:
            self._progress_dialog.close()
            self._progress_dialog = None
        self._import_worker = None
        # New data invalidates the cached fingerprints.
        self.search_engine.invalidate_cache()
        self._refresh_stats()
        self.statusBar().showMessage(f"Import complete: {count} activity records.")
        QMessageBox.information(
            self, "Import complete", f"Imported {count} activity records."
        )

    def _on_import_error(self, message: str) -> None:
        if self._progress_dialog:
            self._progress_dialog.close()
            self._progress_dialog = None
        self._import_worker = None
        logger.error("Import error: %s", message)
        QMessageBox.critical(self, "Import error", f"Import failed:\n\n{message}")

    # ── Export ────────────────────────────────────────────────────────────

    def _export(self, fmt: str) -> None:
        if not self._result or not self._result.predictions:
            QMessageBox.information(self, "Nothing to export", "Run a prediction first.")
            return
        from reports.exporter import ResultExporter

        filters = {
            "csv": ("CSV files (*.csv)", ".csv"),
            "excel": ("Excel files (*.xlsx)", ".xlsx"),
            "html": ("HTML files (*.html)", ".html"),
            "pdf": ("PDF files (*.pdf)", ".pdf"),
        }
        flt, ext = filters[fmt]
        default = os.path.join(self.project_root, f"opentargetai_report{ext}")
        path, _ = QFileDialog.getSaveFileName(self, "Export results", default, flt)
        if not path:
            return
        try:
            exporter = ResultExporter(self._result, self._query_smiles)
            method = {
                "csv": exporter.export_csv,
                "excel": exporter.export_excel,
                "html": exporter.export_html,
                "pdf": exporter.export_pdf,
            }[fmt]
            method(path)
            self.statusBar().showMessage(f"Exported to {path}")
            QMessageBox.information(self, "Export complete", f"Saved to:\n{path}")
        except Exception as exc:  # noqa: BLE001
            logger.exception("Export failed")
            QMessageBox.critical(self, "Export error", f"Export failed:\n\n{exc}")

    # ── Database utilities ────────────────────────────────────────────────

    def _refresh_stats(self) -> None:
        stats = self.db.get_database_stats()
        self.stats_label.setText(
            f"Compounds: {stats['compounds']:,}  |  "
            f"Targets: {stats['targets']:,}  |  "
            f"Activities: {stats['activities']:,}"
        )

    def _show_stats(self) -> None:
        stats = self.db.get_database_stats()
        organisms = "\n".join(
            f"  • {o['organism'] or 'Unknown'}: {o['cnt']:,}"
            for o in stats["organisms"]
        )
        QMessageBox.information(
            self,
            "Database Statistics",
            f"Compounds:   {stats['compounds']:,}\n"
            f"Targets:     {stats['targets']:,}\n"
            f"Activities:  {stats['activities']:,}\n\n"
            f"Top organisms:\n{organisms or '  (none)'}",
        )

    def _clear_cache(self) -> None:
        self.db.clear_prediction_cache()
        self.search_engine.invalidate_cache()
        self.statusBar().showMessage("Prediction cache cleared.")

    # ── Misc ──────────────────────────────────────────────────────────────

    def offer_demo_if_empty(self) -> None:
        """If the database is empty, offer to load the demo dataset.

        Must be called *after* the main window is shown and the splash screen
        has closed (the dialog is modal and would otherwise be hidden behind the
        always-on-top splash).
        """
        if self.db.get_database_stats()["compounds"] == 0:
            if self._confirm_demo():
                self._start_import("demo")

    def _confirm_demo(self) -> bool:
        reply = QMessageBox.question(
            self,
            "Load demo data?",
            "The database is empty. Load the curated demo dataset (37 drugs, "
            "11 targets) so you can try the application immediately?\n\n"
            "You can import full ChEMBL/BindingDB exports later from the "
            "Database menu.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        return reply == QMessageBox.StandardButton.Yes

    def _show_about(self) -> None:
        AboutDialog(APP_VERSION, self).exec()

    def closeEvent(self, event) -> None:  # noqa: N802 - Qt override
        if self._prediction_worker is not None:
            self._prediction_worker.cancel()
            self._prediction_worker.wait(2000)
        event.accept()
