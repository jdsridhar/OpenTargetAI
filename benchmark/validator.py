"""
Benchmark / Validation module for OpenTargetAI.

Measures prediction performance with leave-one-out (LOO) cross-validation: for
each query compound the compound itself is masked from the neighbour set, the
pipeline predicts targets, and the ranked predictions are compared against the
compound's experimentally known targets.

Reported metrics
----------------
* Top-1 / Top-5 / Top-10 accuracy (a known target appears within the top-k)
* Mean Reciprocal Rank (MRR)
* Precision@10 and Recall@10 (macro-averaged)
* ROC-AUC over the pooled (label, confidence) pairs across all queries
"""

from __future__ import annotations

import logging
import random
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

from database.db_manager import DatabaseManager
from prediction.prediction_engine import TargetPredictionEngine
from similarity.fingerprint_engine import FingerprintEngine, normalize_fp_type
from similarity.search_engine import SimilaritySearchEngine

logger = logging.getLogger(__name__)

ProgressCallback = Callable[[int, int, str], None]


@dataclass
class BenchmarkResult:
    """Container for benchmark metrics and configuration."""

    total_queries: int = 0
    evaluated: int = 0
    top1_accuracy: float = 0.0
    top5_accuracy: float = 0.0
    top10_accuracy: float = 0.0
    mean_reciprocal_rank: float = 0.0
    mean_precision_at_10: float = 0.0
    mean_recall_at_10: float = 0.0
    roc_auc: Optional[float] = None
    fp_type: str = "ecfp4"
    metric: str = "tanimoto"
    threshold: float = 0.4
    organism: str = "Homo sapiens"
    error: Optional[str] = None

    def as_dict(self) -> dict:
        """Return a plain dict (used by workers / the UI / report writers)."""
        return {
            "total_queries": self.total_queries,
            "evaluated": self.evaluated,
            "top1_accuracy": self.top1_accuracy,
            "top5_accuracy": self.top5_accuracy,
            "top10_accuracy": self.top10_accuracy,
            "mean_reciprocal_rank": self.mean_reciprocal_rank,
            "mean_precision_at_10": self.mean_precision_at_10,
            "mean_recall_at_10": self.mean_recall_at_10,
            "roc_auc": self.roc_auc,
            "fp_type": self.fp_type,
            "metric": self.metric,
            "threshold": self.threshold,
            "organism": self.organism,
            "error": self.error,
        }


class Validator:
    """Benchmarks the prediction system using leave-one-out validation."""

    def __init__(
        self,
        db_manager: DatabaseManager,
        fp_type: str = "ecfp4",
        metric: str = "tanimoto",
        threshold: float = 0.4,
        top_n: int = 20,
        organism: str = "Homo sapiens",
        max_queries: int = 100,
        seed: int = 42,
    ) -> None:
        self.db = db_manager
        self.fp_type = normalize_fp_type(fp_type)
        self.metric = metric
        self.threshold = threshold
        self.top_n = top_n
        self.organism = organism
        self.max_queries = max_queries
        self.seed = seed

        self.search_engine = SimilaritySearchEngine(db_manager)
        self.pred_engine = TargetPredictionEngine()
        self.fp_engine = FingerprintEngine()

    def run_benchmark(
        self, progress_callback: Optional[ProgressCallback] = None
    ) -> dict:
        """Execute the LOO benchmark and return a metrics dict."""
        cb = progress_callback or (lambda c, t, m: None)
        result = BenchmarkResult(
            fp_type=self.fp_type,
            metric=self.metric,
            threshold=self.threshold,
            organism=self.organism,
        )

        with self.db.get_connection() as conn:
            rows = conn.execute(
                """
                SELECT DISTINCT c.compound_id, c.smiles
                FROM compounds c
                JOIN activities a ON c.compound_id = a.compound_id
                WHERE c.smiles IS NOT NULL
                """
            ).fetchall()

        if not rows:
            result.error = "No compounds with known targets in database."
            return result.as_dict()

        queries = list(rows)
        rng = random.Random(self.seed)
        rng.shuffle(queries)
        queries = queries[: self.max_queries]
        total = len(queries)
        result.total_queries = total

        if total == 0:
            result.error = "No suitable compounds found for benchmarking."
            return result.as_dict()

        top1 = top5 = top10 = 0
        reciprocal_ranks: list[float] = []
        precisions: list[float] = []
        recalls: list[float] = []
        # Pooled (label, score) pairs for a global ROC-AUC.
        roc_labels: list[int] = []
        roc_scores: list[float] = []

        for i, row in enumerate(queries):
            compound_id = row["compound_id"]
            smiles = row["smiles"]
            cb(i, total, f"Benchmarking compound {i + 1}/{total}")

            try:
                true_targets = set(self._get_true_targets(compound_id))
                if not true_targets:
                    continue

                # Leave-one-out: mask the query compound from the neighbour set.
                search_result = self.search_engine.search(
                    query_smiles=smiles,
                    fp_type=self.fp_type,
                    metric=self.metric,
                    threshold=self.threshold,
                    top_n=self.top_n,
                    organism=self.organism if self.organism != "All" else None,
                    exclude_ids={compound_id},
                )

                pred_result = self.pred_engine.predict(search_result)
                predictions = pred_result.predictions
                predicted_ids = [p.target_id for p in predictions]
                if not predicted_ids:
                    continue

                rank = self._first_hit_rank(predicted_ids, true_targets)
                if rank:
                    top1 += rank <= 1
                    top5 += rank <= 5
                    top10 += rank <= 10
                    reciprocal_ranks.append(1.0 / rank)
                else:
                    reciprocal_ranks.append(0.0)

                top10_pred = set(predicted_ids[:10])
                tp = len(top10_pred & true_targets)
                precisions.append(tp / len(top10_pred) if top10_pred else 0.0)
                recalls.append(tp / len(true_targets) if true_targets else 0.0)

                for pred in predictions:
                    roc_labels.append(1 if pred.target_id in true_targets else 0)
                    roc_scores.append(pred.confidence_score)

            except Exception as exc:  # noqa: BLE001
                logger.debug("Benchmark error on compound %s: %s", compound_id, exc)
                continue

        cb(total, total, "Computing metrics...")

        evaluated = len(reciprocal_ranks)
        if evaluated == 0:
            result.error = "No predictions could be evaluated."
            return result.as_dict()

        result.evaluated = evaluated
        result.top1_accuracy = top1 / evaluated
        result.top5_accuracy = top5 / evaluated
        result.top10_accuracy = top10 / evaluated
        result.mean_reciprocal_rank = float(np.mean(reciprocal_ranks))
        result.mean_precision_at_10 = float(np.mean(precisions)) if precisions else 0.0
        result.mean_recall_at_10 = float(np.mean(recalls)) if recalls else 0.0
        result.roc_auc = self._roc_auc(roc_labels, roc_scores)

        logger.info("Benchmark complete: %s", result.as_dict())
        return result.as_dict()

    def _get_true_targets(self, compound_id: int) -> list[int]:
        """Return the known target IDs for a compound."""
        with self.db.get_connection() as conn:
            rows = conn.execute(
                "SELECT DISTINCT target_id FROM activities WHERE compound_id=?",
                (compound_id,),
            ).fetchall()
        return [r["target_id"] for r in rows]

    @staticmethod
    def _first_hit_rank(predicted_ids: list[int], true_set: set[int]) -> int:
        """Return the 1-indexed rank of the first correct prediction (0 if none)."""
        for i, tid in enumerate(predicted_ids, 1):
            if tid in true_set:
                return i
        return 0

    @staticmethod
    def _roc_auc(labels: list[int], scores: list[float]) -> Optional[float]:
        """Compute ROC-AUC over pooled (label, score) pairs, if well defined."""
        if not labels or len(set(labels)) < 2:
            return None
        try:
            from sklearn.metrics import roc_auc_score

            return float(roc_auc_score(labels, scores))
        except Exception as exc:  # noqa: BLE001
            logger.debug("ROC-AUC computation failed: %s", exc)
            return None


def format_report(metrics: dict) -> str:
    """Render a benchmark metrics dict as a human-readable text report."""
    if metrics.get("error"):
        return f"Benchmark could not be run: {metrics['error']}"

    def pct(x: Optional[float]) -> str:
        return f"{x * 100:.1f}%" if x is not None else "N/A"

    def num(x: Optional[float]) -> str:
        return f"{x:.3f}" if x is not None else "N/A"

    lines = [
        "OpenTargetAI — Validation Report",
        "=" * 40,
        f"Fingerprint:          {metrics['fp_type']}",
        f"Similarity metric:    {metrics['metric']}",
        f"Threshold:            {metrics['threshold']}",
        f"Organism:             {metrics['organism']}",
        f"Queries sampled:      {metrics['total_queries']}",
        f"Queries evaluated:    {metrics['evaluated']}",
        "-" * 40,
        f"Top-1 accuracy:       {pct(metrics['top1_accuracy'])}",
        f"Top-5 accuracy:       {pct(metrics['top5_accuracy'])}",
        f"Top-10 accuracy:      {pct(metrics['top10_accuracy'])}",
        f"Mean Reciprocal Rank: {num(metrics['mean_reciprocal_rank'])}",
        f"Precision@10:         {num(metrics['mean_precision_at_10'])}",
        f"Recall@10:            {num(metrics['mean_recall_at_10'])}",
        f"ROC-AUC:              {num(metrics.get('roc_auc'))}",
        "=" * 40,
    ]
    return "\n".join(lines)
