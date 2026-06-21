"""
Background Workers for OpenTargetAI.
QThread-based workers for non-blocking UI operations.
"""

import logging
import traceback
from typing import Any, Callable, Optional

from PyQt6.QtCore import QThread, pyqtSignal

logger = logging.getLogger(__name__)


class PredictionWorker(QThread):
    """Worker thread for running target prediction."""

    progress = pyqtSignal(int, int, str)       # current, total, message
    finished = pyqtSignal(object)               # PredictionResult
    error = pyqtSignal(str)                     # error message

    def __init__(
        self,
        search_engine,
        prediction_engine,
        query_smiles: str,
        fp_type: str,
        metric: str,
        threshold: float,
        top_n: int,
        organism: str,
        min_ligands: int = 1,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.search_engine = search_engine
        self.prediction_engine = prediction_engine
        self.query_smiles = query_smiles
        self.fp_type = fp_type
        self.metric = metric
        self.threshold = threshold
        self.top_n = top_n
        self.organism = organism
        self.min_ligands = min_ligands
        self._cancelled = False

    def cancel(self) -> None:
        """Request cancellation of the current operation."""
        self._cancelled = True

    def run(self) -> None:
        """Execute search and prediction in background thread."""
        try:
            self.progress.emit(0, 100, "Starting similarity search...")

            search_result = self.search_engine.search(
                query_smiles=self.query_smiles,
                fp_type=self.fp_type,
                metric=self.metric,
                threshold=self.threshold,
                top_n=self.top_n,
                organism=self.organism,
                progress_callback=self._emit_progress,
            )

            if self._cancelled:
                return

            self.progress.emit(80, 100, "Aggregating target evidence...")

            prediction_result = self.prediction_engine.predict(
                search_result, min_ligands=self.min_ligands
            )

            if self._cancelled:
                return

            self.progress.emit(100, 100, "Complete")
            self.finished.emit(prediction_result)

        except Exception as e:
            logger.exception("PredictionWorker error")
            self.error.emit(f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    def _emit_progress(self, current: int, total: int, message: str) -> None:
        if not self._cancelled:
            scaled = int((current / max(1, total)) * 70)
            self.progress.emit(scaled, 100, message)


class ImportWorker(QThread):
    """Worker thread for database import operations."""

    progress = pyqtSignal(int, int, str)
    finished = pyqtSignal(int)           # records imported
    error = pyqtSignal(str)

    def __init__(
        self,
        importer,
        import_type: str,    # 'chembl', 'bindingdb', 'demo'
        file_path: Optional[str] = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.importer = importer
        self.import_type = import_type
        self.file_path = file_path

    def run(self) -> None:
        """Run import in background thread."""
        try:
            def callback(cur, tot, msg):
                self.progress.emit(cur, tot, msg)

            self.importer.progress_callback = callback

            if self.import_type == "chembl":
                count = self.importer.import_chembl_activities_tsv(self.file_path)
            elif self.import_type == "bindingdb":
                count = self.importer.import_bindingdb_tsv(self.file_path)
            elif self.import_type == "demo":
                count = self.importer.import_demo_data()
            else:
                raise ValueError(f"Unknown import type: {self.import_type}")

            self.finished.emit(count)

        except Exception as e:
            logger.exception("ImportWorker error")
            self.error.emit(f"{type(e).__name__}: {e}")


class ValidationWorker(QThread):
    """Worker thread for benchmark validation."""

    progress = pyqtSignal(int, int, str)
    finished = pyqtSignal(dict)          # metrics dict
    error = pyqtSignal(str)

    def __init__(
        self,
        validator,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.validator = validator

    def run(self) -> None:
        try:
            def callback(cur, tot, msg):
                self.progress.emit(cur, tot, msg)

            metrics = self.validator.run_benchmark(progress_callback=callback)
            self.finished.emit(metrics)

        except Exception as e:
            logger.exception("ValidationWorker error")
            self.error.emit(f"{type(e).__name__}: {e}")


class GenericWorker(QThread):
    """General-purpose worker for arbitrary callable tasks."""

    progress = pyqtSignal(int, int, str)
    finished = pyqtSignal(object)
    error = pyqtSignal(str)

    def __init__(self, fn: Callable, *args, **kwargs) -> None:
        super().__init__()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs

    def run(self) -> None:
        try:
            result = self.fn(*self.args, **self.kwargs)
            self.finished.emit(result)
        except Exception as e:
            logger.exception("GenericWorker error")
            self.error.emit(f"{type(e).__name__}: {e}")
