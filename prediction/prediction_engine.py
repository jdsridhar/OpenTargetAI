"""
Target Prediction Engine for OpenTargetAI.

Aggregates evidence from chemically similar compounds to predict likely protein
targets for a query molecule (ligand-based target fishing).

Algorithm
---------
1. Retrieve similar compounds from the similarity search.
2. Collect all known targets of those similar compounds.
3. For each target, aggregate similarity-weighted bioactivity evidence.
4. Score each target and derive an interpretable confidence (see
   :class:`prediction.confidence.ConfidenceModel`).
5. Rank, normalise and classify the predictions.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from prediction.confidence import ConfidenceModel
from similarity.search_engine import SearchResult, SimilarCompound

logger = logging.getLogger(__name__)


@dataclass
class TargetPrediction:
    """Represents a predicted target with aggregated, explainable evidence."""

    target_id: int
    target_name: str
    gene_symbol: str
    uniprot: str
    organism: str
    target_class: str

    # Scoring
    confidence_score: float = 0.0        # 0–1 overall confidence
    target_score: float = 0.0            # weighted evidence score
    enrichment_score: float = 0.0        # enrichment vs. random expectation
    normalized_rank: float = 0.0         # 0–1, lower = better rank

    # Evidence summary
    n_supporting_ligands: int = 0
    mean_similarity: float = 0.0
    max_similarity: float = 0.0
    mean_pchembl: Optional[float] = None
    activity_types: list[str] = field(default_factory=list)
    supporting_compounds: list[dict] = field(default_factory=list)

    # Confidence components (for the explainer panel)
    similarity_weight: float = 0.0
    bioactivity_weight: float = 0.0
    occurrence_frequency: float = 0.0
    activity_consistency: float = 0.0

    # Classification
    confidence_label: str = "Low"        # High / Medium / Low

    @property
    def evidence_count(self) -> int:
        """Total number of supporting activity records."""
        return len(self.supporting_compounds)

    def activity_stats(self) -> dict[str, float]:
        """Summary statistics of the supporting bioactivities (pChEMBL units)."""
        vals = [
            c["pchembl_value"]
            for c in self.supporting_compounds
            if c.get("pchembl_value")
        ]
        if not vals:
            return {}
        arr = np.array(vals, dtype=float)
        return {
            "n": int(arr.size),
            "mean": float(arr.mean()),
            "median": float(np.median(arr)),
            "min": float(arr.min()),
            "max": float(arr.max()),
            "std": float(arr.std()),
        }


@dataclass
class PredictionResult:
    """Complete prediction output for a single query molecule."""

    query_smiles: str
    predictions: list[TargetPrediction] = field(default_factory=list)
    search_result: Optional[SearchResult] = None
    n_similar_compounds: int = 0
    n_targets_considered: int = 0


class TargetPredictionEngine:
    """Core prediction algorithm implementing weighted ligand-based target fishing."""

    def __init__(self, confidence_model: Optional[ConfidenceModel] = None) -> None:
        self.confidence_model = confidence_model or ConfidenceModel()

    def predict(
        self,
        search_result: SearchResult,
        min_ligands: int = 1,
    ) -> PredictionResult:
        """
        Generate target predictions from a similarity search result.

        Args:
            search_result: output from :meth:`SimilaritySearchEngine.search`.
            min_ligands: minimum number of distinct supporting ligands required
                for a target to be reported.
        """
        similar = search_result.similar_compounds
        if not similar:
            logger.warning("No similar compounds – cannot predict targets")
            return PredictionResult(
                query_smiles=search_result.query_smiles,
                search_result=search_result,
            )

        # Step 1: aggregate all target evidence across the neighbours.
        target_evidence: dict[int, dict] = {}
        for compound in similar:
            self._aggregate_compound_evidence(compound, target_evidence)

        if not target_evidence:
            logger.info("No targets found among similar compounds")
            return PredictionResult(
                query_smiles=search_result.query_smiles,
                search_result=search_result,
                n_similar_compounds=len(similar),
            )

        # Step 2: score every candidate target.
        predictions = self._score_targets(target_evidence, len(similar))

        # Step 3: filter by minimum supporting ligands.
        predictions = [
            p for p in predictions if p.n_supporting_ligands >= min_ligands
        ]

        # Step 4: rank and normalise.
        predictions = self._rank_predictions(predictions)

        logger.info(
            "Predicted %d targets from %d similar compounds",
            len(predictions),
            len(similar),
        )

        return PredictionResult(
            query_smiles=search_result.query_smiles,
            predictions=predictions,
            search_result=search_result,
            n_similar_compounds=len(similar),
            n_targets_considered=len(target_evidence),
        )

    def _aggregate_compound_evidence(
        self,
        compound: SimilarCompound,
        target_evidence: dict[int, dict],
    ) -> None:
        """Accumulate evidence for each target from one similar compound."""
        for act in compound.targets:
            tid = act["target_id"]
            ev = target_evidence.get(tid)
            if ev is None:
                ev = target_evidence[tid] = {
                    "target_id": tid,
                    "target_name": act.get("target_name", ""),
                    "gene_symbol": act.get("gene_symbol", ""),
                    "uniprot": act.get("uniprot", ""),
                    "organism": act.get("organism", ""),
                    "target_class": act.get("target_class", ""),
                    "similarity_scores": [],
                    "pchembl_values": [],
                    "activity_types": set(),
                    "supporting_compounds": [],
                }
            ev["similarity_scores"].append(compound.similarity)
            pchembl = act.get("pchembl_value")
            if pchembl and pchembl > 0:
                ev["pchembl_values"].append(pchembl)
            atype = act.get("activity_type", "")
            if atype:
                ev["activity_types"].add(atype)
            ev["supporting_compounds"].append(
                {
                    "compound_id": compound.compound_id,
                    "smiles": compound.smiles,
                    "chembl_id": compound.chembl_id,
                    "similarity": compound.similarity,
                    "activity_type": atype,
                    "activity_value": act.get("activity_value"),
                    "units": act.get("units", ""),
                    "pchembl_value": pchembl,
                }
            )

    def _score_targets(
        self,
        target_evidence: dict[int, dict],
        total_similar: int,
    ) -> list[TargetPrediction]:
        """Compute weighted scores and confidence for every candidate target."""
        predictions: list[TargetPrediction] = []

        for tid, ev in target_evidence.items():
            sims = ev["similarity_scores"]
            pchembl_vals = ev["pchembl_values"]
            n_ligands = len({c["compound_id"] for c in ev["supporting_compounds"]})

            assessment = self.confidence_model.evaluate(
                similarity_scores=sims,
                pchembl_values=pchembl_vals,
                n_ligands=n_ligands,
                total_neighbours=total_similar,
            )

            # Weighted target score = similarity × bioactivity × occurrence.
            target_score = (
                assessment.similarity_weight
                * assessment.bioactivity_weight
                * assessment.occurrence_frequency
            )
            enrichment_score = self.confidence_model.enrichment(
                target_score, n_ligands
            )

            predictions.append(
                TargetPrediction(
                    target_id=tid,
                    target_name=ev["target_name"],
                    gene_symbol=ev["gene_symbol"],
                    uniprot=ev["uniprot"],
                    organism=ev["organism"],
                    target_class=ev["target_class"],
                    confidence_score=assessment.confidence_score,
                    target_score=target_score,
                    enrichment_score=enrichment_score,
                    n_supporting_ligands=n_ligands,
                    mean_similarity=assessment.mean_similarity,
                    max_similarity=assessment.max_similarity,
                    mean_pchembl=assessment.mean_pchembl,
                    activity_types=sorted(ev["activity_types"]),
                    supporting_compounds=sorted(
                        ev["supporting_compounds"],
                        key=lambda c: c["similarity"],
                        reverse=True,
                    ),
                    similarity_weight=assessment.similarity_weight,
                    bioactivity_weight=assessment.bioactivity_weight,
                    occurrence_frequency=assessment.occurrence_frequency,
                    activity_consistency=assessment.activity_consistency,
                    confidence_label=assessment.confidence_label,
                )
            )

        return predictions

    @staticmethod
    def _rank_predictions(
        predictions: list[TargetPrediction],
    ) -> list[TargetPrediction]:
        """Sort predictions by confidence and assign normalised ranks."""
        if not predictions:
            return []
        predictions.sort(
            key=lambda p: (p.confidence_score, p.enrichment_score), reverse=True
        )
        n = len(predictions)
        for i, pred in enumerate(predictions):
            pred.normalized_rank = i / n
        return predictions
