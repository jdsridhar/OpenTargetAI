"""
Similarity Search Engine for OpenTargetAI.
Searches the local database for compounds similar to a query molecule.
"""

import logging
from dataclasses import dataclass, field
from typing import Optional

from rdkit import Chem, DataStructs

from database.db_manager import DatabaseManager
from similarity.fingerprint_engine import FingerprintEngine

logger = logging.getLogger(__name__)


@dataclass
class SimilarCompound:
    """Represents a database compound with its similarity score."""
    compound_id: int
    smiles: str
    chembl_id: Optional[str]
    similarity: float
    targets: list[dict] = field(default_factory=list)


@dataclass
class SearchResult:
    """Complete search result set."""
    query_smiles: str
    fp_type: str
    metric: str
    threshold: float
    similar_compounds: list[SimilarCompound] = field(default_factory=list)
    total_searched: int = 0
    search_time_ms: float = 0.0


class SimilaritySearchEngine:
    """
    Performs fingerprint-based similarity searches against the local database.
    Supports Tanimoto, Dice, and Cosine metrics with configurable thresholds.
    """

    def __init__(self, db_manager: DatabaseManager) -> None:
        self.db = db_manager
        self.fp_engine = FingerprintEngine()
        self._fp_cache: dict[str, list[tuple[int, DataStructs.ExplicitBitVect]]] = {}

    def search(
        self,
        query_smiles: str,
        fp_type: str = "ecfp4",
        metric: str = "tanimoto",
        threshold: float = 0.4,
        top_n: int = 50,
        organism: Optional[str] = None,
        progress_callback=None,
        exclude_ids: Optional[set[int]] = None,
    ) -> SearchResult:
        """
        Main search entry point.

        Args:
            exclude_ids: compound IDs to omit from the neighbour set. Used by the
                benchmark module to perform leave-one-out validation (so a query
                compound cannot retrieve itself).

        Returns a :class:`SearchResult` with ranked similar compounds.
        """
        import time
        t0 = time.time()

        mol = Chem.MolFromSmiles(query_smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {query_smiles}")

        query_fp = self.fp_engine.generate_fp_object(mol, fp_type)
        if query_fp is None:
            raise ValueError(f"Could not generate {fp_type} fingerprint for query")

        # Load or cache database fingerprints
        db_entries = self._load_fingerprints(fp_type, progress_callback)
        total = len(db_entries)

        if not db_entries:
            logger.warning("No fingerprints found in database. Import data first.")
            return SearchResult(
                query_smiles=query_smiles,
                fp_type=fp_type,
                metric=metric,
                threshold=threshold,
                total_searched=0,
                search_time_ms=(time.time() - t0) * 1000,
            )

        # Bulk similarity computation
        compound_ids = [e[0] for e in db_entries]
        db_fps = [e[1] for e in db_entries]

        if progress_callback:
            progress_callback(0, total, "Computing similarities...")

        similarities = self.fp_engine.bulk_similarity(query_fp, db_fps, metric)

        if progress_callback:
            progress_callback(total, total, "Filtering results...")

        # Filter and sort (optionally masking specific compounds for LOO validation)
        exclude_ids = exclude_ids or set()
        hits = [
            (compound_ids[i], similarities[i])
            for i in range(len(similarities))
            if similarities[i] >= threshold and compound_ids[i] not in exclude_ids
        ]
        hits.sort(key=lambda x: x[1], reverse=True)
        hits = hits[:top_n]

        # Fetch compound details and targets
        similar_compounds = self._enrich_hits(hits, organism)

        elapsed = (time.time() - t0) * 1000
        logger.info(
            f"Similarity search: {len(similar_compounds)} hits "
            f"from {total} compounds in {elapsed:.0f}ms"
        )

        return SearchResult(
            query_smiles=query_smiles,
            fp_type=fp_type,
            metric=metric,
            threshold=threshold,
            similar_compounds=similar_compounds,
            total_searched=total,
            search_time_ms=elapsed,
        )

    def _load_fingerprints(
        self,
        fp_type: str,
        progress_callback=None,
    ) -> list[tuple[int, DataStructs.ExplicitBitVect]]:
        """
        Load all fingerprints of the given type from the database.
        Results are cached in memory for repeated searches.
        """
        if fp_type in self._fp_cache:
            return self._fp_cache[fp_type]

        raw = self.db.get_all_fingerprints(fp_type)
        total = len(raw)

        if progress_callback:
            progress_callback(0, total, f"Loading {total} fingerprints...")

        entries = []
        for i, row in enumerate(raw):
            fp = self.fp_engine.bytes_to_fp(row["fingerprint"])
            if fp is not None:
                entries.append((row["compound_id"], fp))
            if progress_callback and i % 10000 == 0:
                progress_callback(i, total, f"Loading fingerprints: {i}/{total}")

        self._fp_cache[fp_type] = entries
        logger.info(f"Loaded {len(entries)} valid {fp_type} fingerprints")
        return entries

    def _enrich_hits(
        self,
        hits: list[tuple[int, float]],
        organism: Optional[str],
    ) -> list[SimilarCompound]:
        """Fetch compound SMILES and associated targets for each hit."""
        if not hits:
            return []

        compound_ids = [h[0] for h in hits]
        sim_map = {h[0]: h[1] for h in hits}

        # Get all target activities for these compounds
        activities = self.db.get_compound_targets(compound_ids)

        # Filter by organism if specified
        if organism and organism != "All":
            activities = [
                a for a in activities
                if organism.lower() in a.get("organism", "").lower()
            ]

        # Group activities by compound_id
        targets_by_compound: dict[int, list[dict]] = {}
        for act in activities:
            cid = act["compound_id"]
            if cid not in targets_by_compound:
                targets_by_compound[cid] = []
            targets_by_compound[cid].append(act)

        # Build result objects
        results = []
        for cid, sim in hits:
            compound = self.db.get_compound_by_id(cid)
            if compound is None:
                continue
            results.append(SimilarCompound(
                compound_id=cid,
                smiles=compound.get("smiles", ""),
                chembl_id=compound.get("chembl_id"),
                similarity=sim,
                targets=targets_by_compound.get(cid, []),
            ))

        return results

    def invalidate_cache(self) -> None:
        """Clear the in-memory fingerprint cache."""
        self._fp_cache.clear()
        logger.info("Fingerprint cache cleared")
