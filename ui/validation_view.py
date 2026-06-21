"""
Validation tab for OpenTargetAI.

Configures and runs the leave-one-out benchmark, then displays the headline
metrics (Top-k accuracy, MRR, Precision/Recall, ROC-AUC) and a saveable text
report. The benchmark runs in a background worker.
"""

from __future__ import annotations

import os
from typing import Optional

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from benchmark.validator import Validator, format_report
from ui.prediction_controls import FINGERPRINTS, METRICS
from ui.widgets import SectionTitle, StatCard
from workers.workers import ValidationWorker


class ValidationView(QWidget):
    """Benchmark configuration, execution and reporting."""

    status_message = pyqtSignal(str)

    def __init__(self, db_manager, project_root: str, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.db = db_manager
        self.project_root = project_root
        self._worker: Optional[ValidationWorker] = None
        self._last_metrics: Optional[dict] = None
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.addWidget(SectionTitle("Validation / Benchmark"))

        intro = QLabel(
            "Leave-one-out cross-validation: each compound is masked from the "
            "neighbour set, its targets are predicted, and the ranking is "
            "compared with its known targets."
        )
        intro.setWordWrap(True)
        intro.setStyleSheet("color:#8b949e;")
        layout.addWidget(intro)

        # Config
        config_group = QGroupBox("Configuration")
        form = QFormLayout(config_group)
        self.fp_combo = QComboBox()
        for label, key in FINGERPRINTS:
            self.fp_combo.addItem(label, key)
        form.addRow("Fingerprint:", self.fp_combo)
        self.metric_combo = QComboBox()
        for label, key in METRICS:
            self.metric_combo.addItem(label, key)
        form.addRow("Metric:", self.metric_combo)
        self.threshold_spin = QSpinBox()
        self.threshold_spin.setRange(5, 95)
        self.threshold_spin.setValue(30)
        self.threshold_spin.setSuffix(" %")
        form.addRow("Threshold:", self.threshold_spin)
        self.max_queries_spin = QSpinBox()
        self.max_queries_spin.setRange(10, 5000)
        self.max_queries_spin.setValue(100)
        form.addRow("Max queries:", self.max_queries_spin)
        layout.addWidget(config_group)

        run_row = QHBoxLayout()
        self.run_button = QPushButton("▶  Run Benchmark")
        self.run_button.setObjectName("runButton")
        self.run_button.clicked.connect(self._run)
        self.save_button = QPushButton("Save Report")
        self.save_button.clicked.connect(self._save_report)
        self.save_button.setEnabled(False)
        run_row.addWidget(self.run_button)
        run_row.addWidget(self.save_button)
        run_row.addStretch(1)
        layout.addLayout(run_row)

        self.progress = QProgressBar()
        self.progress.setVisible(False)
        layout.addWidget(self.progress)

        # Metric cards
        cards = QGridLayout()
        self.card_top1 = StatCard("Top-1 accuracy")
        self.card_top5 = StatCard("Top-5 accuracy")
        self.card_top10 = StatCard("Top-10 accuracy")
        self.card_mrr = StatCard("Mean Reciprocal Rank")
        self.card_roc = StatCard("ROC-AUC")
        self.card_eval = StatCard("Queries evaluated")
        for w, (r, c) in (
            (self.card_top1, (0, 0)),
            (self.card_top5, (0, 1)),
            (self.card_top10, (0, 2)),
            (self.card_mrr, (1, 0)),
            (self.card_roc, (1, 1)),
            (self.card_eval, (1, 2)),
        ):
            cards.addWidget(w, r, c)
        layout.addLayout(cards)

        # Report text
        report_group = QGroupBox("Report")
        rg = QVBoxLayout(report_group)
        self.report = QPlainTextEdit()
        self.report.setReadOnly(True)
        self.report.setStyleSheet("font-family:Consolas,monospace;")
        rg.addWidget(self.report)
        layout.addWidget(report_group, 1)

    # ── Run ───────────────────────────────────────────────────────────────

    def _run(self) -> None:
        if self._worker is not None:
            return
        validator = Validator(
            self.db,
            fp_type=self.fp_combo.currentData(),
            metric=self.metric_combo.currentData(),
            threshold=self.threshold_spin.value() / 100.0,
            max_queries=self.max_queries_spin.value(),
        )
        self.run_button.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 100)
        self.report.setPlainText("Running benchmark …")

        self._worker = ValidationWorker(validator)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_finished)
        self._worker.error.connect(self._on_error)
        self._worker.start()
        self.status_message.emit("Benchmark running…")

    def _on_progress(self, current: int, total: int, message: str) -> None:
        pct = int((current / total) * 100) if total else 0
        self.progress.setValue(min(100, pct))

    def _on_finished(self, metrics: dict) -> None:
        self._worker = None
        self.run_button.setEnabled(True)
        self.progress.setVisible(False)
        self._last_metrics = metrics
        self.report.setPlainText(format_report(metrics))
        self.status_message.emit("Benchmark complete")

        if metrics.get("error"):
            for card in (
                self.card_top1, self.card_top5, self.card_top10,
                self.card_mrr, self.card_roc, self.card_eval,
            ):
                card.set_value("—")
            self.save_button.setEnabled(False)
            return

        self.card_top1.set_value(f"{metrics['top1_accuracy'] * 100:.0f}%")
        self.card_top5.set_value(f"{metrics['top5_accuracy'] * 100:.0f}%")
        self.card_top10.set_value(f"{metrics['top10_accuracy'] * 100:.0f}%")
        self.card_mrr.set_value(f"{metrics['mean_reciprocal_rank']:.3f}")
        roc = metrics.get("roc_auc")
        self.card_roc.set_value(f"{roc:.3f}" if roc is not None else "N/A")
        self.card_eval.set_value(str(metrics["evaluated"]))
        self.save_button.setEnabled(True)

    def _on_error(self, message: str) -> None:
        self._worker = None
        self.run_button.setEnabled(True)
        self.progress.setVisible(False)
        self.report.setPlainText(f"Benchmark failed:\n{message}")
        self.status_message.emit("Benchmark failed")

    def _save_report(self) -> None:
        if not self._last_metrics:
            return
        default = os.path.join(self.project_root, "validation_report.txt")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save validation report", default, "Text files (*.txt)"
        )
        if not path:
            return
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(format_report(self._last_metrics))
        self.status_message.emit(f"Report saved to {path}")
