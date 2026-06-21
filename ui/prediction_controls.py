"""
Right panel: prediction controls for OpenTargetAI.

Exposes the species filter, fingerprint type, similarity metric, threshold,
Top-N neighbour count and minimum-ligand support, plus the Run Prediction
action and a progress indicator.
"""

from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QPushButton,
    QSlider,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

# (display label, internal fingerprint key)
FINGERPRINTS = [
    ("Morgan ECFP4", "ecfp4"),
    ("Morgan ECFP6", "ecfp6"),
    ("MACCS Keys", "maccs"),
    ("Atom Pair", "ap"),
    ("Topological", "topo"),
]
METRICS = [("Tanimoto", "tanimoto"), ("Dice", "dice"), ("Cosine", "cosine")]
SPECIES = ["Homo sapiens", "Mus musculus", "Rattus norvegicus", "All"]
SPECIES_LABELS = {
    "Homo sapiens": "Human",
    "Mus musculus": "Mouse",
    "Rattus norvegicus": "Rat",
    "All": "All species",
}


class PredictionControls(QWidget):
    """Similarity / prediction parameter controls."""

    run_requested = pyqtSignal(dict)
    cancel_requested = pyqtSignal()

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self._running = False
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)

        title = QLabel("Prediction Controls")
        title.setObjectName("sectionTitle")
        layout.addWidget(title)

        # Species
        species_group = QGroupBox("Species")
        sg = QVBoxLayout(species_group)
        self.species_combo = QComboBox()
        for org in SPECIES:
            self.species_combo.addItem(SPECIES_LABELS[org], org)
        sg.addWidget(self.species_combo)
        layout.addWidget(species_group)

        # Similarity settings
        sim_group = QGroupBox("Similarity Settings")
        form = QFormLayout(sim_group)

        self.fp_combo = QComboBox()
        for label, key in FINGERPRINTS:
            self.fp_combo.addItem(label, key)
        form.addRow("Fingerprint:", self.fp_combo)

        self.metric_combo = QComboBox()
        for label, key in METRICS:
            self.metric_combo.addItem(label, key)
        form.addRow("Metric:", self.metric_combo)

        # Threshold slider + readout
        thr_row = QHBoxLayout()
        self.threshold_slider = QSlider(Qt.Orientation.Horizontal)
        self.threshold_slider.setRange(5, 95)
        self.threshold_slider.setValue(40)
        self.threshold_value = QLabel("0.40")
        self.threshold_value.setMinimumWidth(36)
        self.threshold_slider.valueChanged.connect(
            lambda v: self.threshold_value.setText(f"{v / 100:.2f}")
        )
        thr_row.addWidget(self.threshold_slider)
        thr_row.addWidget(self.threshold_value)
        form.addRow("Threshold:", thr_row)

        self.top_n_spin = QSpinBox()
        self.top_n_spin.setRange(5, 1000)
        self.top_n_spin.setValue(50)
        form.addRow("Top N similar:", self.top_n_spin)

        self.min_ligands_spin = QSpinBox()
        self.min_ligands_spin.setRange(1, 20)
        self.min_ligands_spin.setValue(1)
        form.addRow("Min. ligands:", self.min_ligands_spin)

        layout.addWidget(sim_group)

        # Run button + progress
        self.run_button = QPushButton("▶  Run Prediction")
        self.run_button.setObjectName("runButton")
        self.run_button.clicked.connect(self._on_run_clicked)
        layout.addWidget(self.run_button)

        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.progress.setVisible(False)
        layout.addWidget(self.progress)

        self.progress_label = QLabel("")
        self.progress_label.setWordWrap(True)
        self.progress_label.setStyleSheet("color:#8b949e;font-size:11px;")
        layout.addWidget(self.progress_label)

        layout.addStretch(1)

    # ── Public API ────────────────────────────────────────────────────────

    def parameters(self) -> dict:
        """Return the current control values as a parameter dict."""
        return {
            "organism": self.species_combo.currentData(),
            "fp_type": self.fp_combo.currentData(),
            "metric": self.metric_combo.currentData(),
            "threshold": self.threshold_slider.value() / 100.0,
            "top_n": self.top_n_spin.value(),
            "min_ligands": self.min_ligands_spin.value(),
        }

    def set_running(self, running: bool) -> None:
        """Toggle between idle and running states."""
        self._running = running
        self.progress.setVisible(running)
        if running:
            self.run_button.setText("■  Cancel")
            self.progress.setValue(0)
        else:
            self.run_button.setText("▶  Run Prediction")
            self.progress_label.setText("")

    def update_progress(self, current: int, total: int, message: str) -> None:
        pct = int((current / total) * 100) if total else 0
        self.progress.setValue(min(100, pct))
        self.progress_label.setText(message)

    # ── Handlers ──────────────────────────────────────────────────────────

    def _on_run_clicked(self) -> None:
        if self._running:
            self.cancel_requested.emit()
        else:
            self.run_requested.emit(self.parameters())
