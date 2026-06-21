"""
Confidence model for OpenTargetAI target predictions.

Encapsulates the logic that turns the aggregated evidence for a candidate
target into an interpretable confidence score (0–1) and a categorical label
(High / Medium / Low). Keeping this separate from the scoring engine makes the
heuristics transparent, testable, and easy to tune.

The confidence reflects four orthogonal lines of evidence:

* **Similarity strength** – how chemically close the supporting ligands are to
  the query (penalising weak matches via a squared-similarity weight).
* **Bioactivity potency** – how potent those ligands are at the target
  (derived from pChEMBL / pActivity values).
* **Occurrence frequency** – what fraction of the retrieved neighbours hit the
  target (recurrence is corroborating evidence).
* **Activity consistency** – how reproducible the reported potencies are
  (tight distributions are more trustworthy than scattered ones).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from statistics import pstdev
from typing import Optional


@dataclass
class ConfidenceAssessment:
    """Result of evaluating the confidence of a single target prediction."""

    confidence_score: float          # composite 0–1 value
    confidence_label: str            # "High" / "Medium" / "Low"

    # Individual evidence components (each 0–1) for explainability.
    similarity_weight: float
    bioactivity_weight: float
    occurrence_frequency: float
    activity_consistency: float

    # Convenience statistics surfaced to the UI / reports.
    mean_similarity: float = 0.0
    max_similarity: float = 0.0
    mean_pchembl: Optional[float] = None
    components: dict[str, float] = field(default_factory=dict)


class ConfidenceModel:
    """
    Computes confidence assessments from per-target evidence.

    The composite score is a weighted blend of the four evidence components,
    deliberately dominated by similarity strength (the primary signal in
    ligand-based target fishing) with potency, recurrence and consistency as
    supporting modifiers.
    """

    HIGH_THRESHOLD = 0.65
    MEDIUM_THRESHOLD = 0.40

    # Component weights (must sum to 1.0).
    W_SIMILARITY = 0.45
    W_BIOACTIVITY = 0.20
    W_OCCURRENCE = 0.20
    W_CONSISTENCY = 0.15

    def evaluate(
        self,
        similarity_scores: list[float],
        pchembl_values: list[float],
        n_ligands: int,
        total_neighbours: int,
    ) -> ConfidenceAssessment:
        """
        Evaluate confidence for one target.

        Args:
            similarity_scores: similarity of every supporting ligand to the query.
            pchembl_values: available pChEMBL/pActivity values for the target.
            n_ligands: number of distinct supporting ligands.
            total_neighbours: total neighbours retrieved by the search.
        """
        sims = [s for s in similarity_scores if s is not None]
        if not sims:
            sims = [0.0]

        mean_sim = sum(sims) / len(sims)
        max_sim = max(sims)

        # Similarity weight: mean of squared similarities rewards close matches
        # and suppresses contributions from weak neighbours.
        similarity_weight = sum(s * s for s in sims) / len(sims)

        # Bioactivity weight from mean pChEMBL (≈5 weak → ≈9 very potent).
        if pchembl_values:
            mean_pchembl = sum(pchembl_values) / len(pchembl_values)
            bioactivity_weight = _clamp((mean_pchembl - 4.0) / 5.0)
        else:
            mean_pchembl = None
            bioactivity_weight = 0.5  # neutral when no potency data

        # Occurrence frequency, gently saturated so a couple of hits already
        # count as meaningful corroboration.
        raw_freq = n_ligands / max(1, total_neighbours)
        occurrence_frequency = _clamp(raw_freq * 5.0)

        # Activity consistency: tight pChEMBL distributions score high.
        activity_consistency = self._consistency(pchembl_values)

        confidence_score = _clamp(
            self.W_SIMILARITY * similarity_weight
            + self.W_BIOACTIVITY * bioactivity_weight
            + self.W_OCCURRENCE * occurrence_frequency
            + self.W_CONSISTENCY * activity_consistency
        )

        return ConfidenceAssessment(
            confidence_score=confidence_score,
            confidence_label=self.classify(confidence_score),
            similarity_weight=similarity_weight,
            bioactivity_weight=bioactivity_weight,
            occurrence_frequency=occurrence_frequency,
            activity_consistency=activity_consistency,
            mean_similarity=mean_sim,
            max_similarity=max_sim,
            mean_pchembl=mean_pchembl,
            components={
                "similarity": similarity_weight,
                "bioactivity": bioactivity_weight,
                "occurrence": occurrence_frequency,
                "consistency": activity_consistency,
            },
        )

    @staticmethod
    def _consistency(pchembl_values: list[float]) -> float:
        """
        Map the spread of pChEMBL values to a 0–1 consistency score. A single
        (or no) measurement is treated as neutral (0.5); tighter distributions
        approach 1.0, scattered ones approach 0.0.
        """
        vals = [v for v in pchembl_values if v is not None]
        if len(vals) < 2:
            return 0.5
        spread = pstdev(vals)
        # 1 log unit of scatter ≈ moderate disagreement.
        return _clamp(1.0 - spread / 2.0)

    def classify(self, score: float) -> str:
        """Map a composite score to a categorical confidence label."""
        if score >= self.HIGH_THRESHOLD:
            return "High"
        if score >= self.MEDIUM_THRESHOLD:
            return "Medium"
        return "Low"

    @staticmethod
    def enrichment(target_score: float, n_ligands: int) -> float:
        """
        A simple enrichment score that rewards recurrent evidence on top of the
        weighted target score (log-scaled to avoid runaway values).
        """
        return target_score * math.log1p(n_ligands)


def _clamp(value: float, low: float = 0.0, high: float = 1.0) -> float:
    """Clamp a value into the inclusive ``[low, high]`` range."""
    return max(low, min(high, value))
