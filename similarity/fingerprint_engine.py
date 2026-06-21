"""
Fingerprint Engine for OpenTargetAI.

Generates and manages molecular fingerprints using RDKit's modern
``rdFingerprintGenerator`` API. Fingerprints are serialised with the compact
native ``ExplicitBitVect.ToBinary`` representation (≈40-260 bytes each) rather
than pickle, which keeps the database small for >2M compounds and avoids the
security concerns of unpickling arbitrary blobs.

Supported fingerprint types
----------------------------
* ``ecfp4`` – Morgan, radius 2, 2048 bits (ECFP4 analogue)
* ``ecfp6`` – Morgan, radius 3, 2048 bits (ECFP6 analogue)
* ``maccs`` – MACCS structural keys, 167 bits
* ``ap``    – Atom-pair fingerprint, folded to 2048 bits
* ``topo``  – RDKit topological (path) fingerprint, 2048 bits
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator

logger = logging.getLogger(__name__)

# Public identifiers used throughout the application and the database columns.
FP_TYPES = ["ecfp4", "ecfp6", "maccs", "ap", "topo"]

FP_DISPLAY_NAMES = {
    "ecfp4": "Morgan ECFP4 (r=2)",
    "ecfp6": "Morgan ECFP6 (r=3)",
    "maccs": "MACCS Keys (167)",
    "ap": "Atom Pair",
    "topo": "Topological (RDKit)",
}

# Map a few human-facing labels back onto the internal keys so the UI can pass
# either form safely.
FP_ALIASES = {
    "morgan ecfp4": "ecfp4",
    "morgan ecfp6": "ecfp6",
    "ecfp4": "ecfp4",
    "ecfp6": "ecfp6",
    "maccs": "maccs",
    "atom pair": "ap",
    "atompair": "ap",
    "ap": "ap",
    "topological": "topo",
    "topo": "topo",
}


def normalize_fp_type(fp_type: str) -> str:
    """Normalise a user/UI supplied fingerprint label to an internal key."""
    return FP_ALIASES.get(fp_type.strip().lower(), fp_type.strip().lower())


class FingerprintEngine:
    """
    Generates Morgan (ECFP), MACCS, atom-pair and topological fingerprints and
    computes pairwise / bulk similarity.

    The RDKit fingerprint *generators* are stateless and reusable, so they are
    created once per engine instance and shared across all calls. The engine is
    therefore cheap to construct and safe to reuse from worker threads (each
    worker should hold its own instance).
    """

    def __init__(self, fp_size: int = 2048) -> None:
        self.fp_size = fp_size
        self._morgan2 = rdFingerprintGenerator.GetMorganGenerator(
            radius=2, fpSize=fp_size
        )
        self._morgan3 = rdFingerprintGenerator.GetMorganGenerator(
            radius=3, fpSize=fp_size
        )
        self._atompair = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=fp_size)
        self._rdkit = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=fp_size)

    # ── Generator dispatch ────────────────────────────────────────────────

    def generate_fp_object(
        self, mol: Chem.Mol, fp_type: str = "ecfp4"
    ) -> Optional[DataStructs.ExplicitBitVect]:
        """
        Generate a fingerprint as an RDKit ``ExplicitBitVect`` object suitable
        for similarity calculation. Returns ``None`` on failure.
        """
        if mol is None:
            return None
        fp_type = normalize_fp_type(fp_type)
        try:
            if fp_type == "ecfp4":
                return self._morgan2.GetFingerprint(mol)
            if fp_type == "ecfp6":
                return self._morgan3.GetFingerprint(mol)
            if fp_type == "maccs":
                return MACCSkeys.GenMACCSKeys(mol)
            if fp_type == "ap":
                return self._atompair.GetFingerprint(mol)
            if fp_type == "topo":
                return self._rdkit.GetFingerprint(mol)
            raise ValueError(f"Unknown fp_type: {fp_type}")
        except Exception as exc:  # noqa: BLE001 - defensive, logged
            logger.debug("FP object generation error (%s): %s", fp_type, exc)
            return None

    def fingerprint_from_smiles(
        self, smiles: str, fp_type: str = "ecfp4"
    ) -> Optional[DataStructs.ExplicitBitVect]:
        """Generate a fingerprint object directly from a SMILES string."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return self.generate_fp_object(mol, fp_type)

    # ── Binary serialisation (single type → bytes) ────────────────────────

    def generate_bytes(self, mol: Chem.Mol, fp_type: str = "ecfp4") -> Optional[bytes]:
        """Generate a single fingerprint and serialise it to bytes for storage."""
        fp = self.generate_fp_object(mol, fp_type)
        return self._fp_to_bytes(fp) if fp is not None else None

    def generate_all_fingerprints(
        self, mol: Chem.Mol
    ) -> dict[str, Optional[bytes]]:
        """Generate every supported fingerprint type for a molecule as bytes."""
        return {fp_type: self.generate_bytes(mol, fp_type) for fp_type in FP_TYPES}

    @staticmethod
    def _fp_to_bytes(fp: DataStructs.ExplicitBitVect) -> bytes:
        """Serialise a fingerprint to its compact native binary form."""
        return fp.ToBinary()

    @staticmethod
    def bytes_to_fp(
        blob: Optional[bytes],
    ) -> Optional[DataStructs.ExplicitBitVect]:
        """
        Deserialise a fingerprint from a stored blob.

        Primary path uses the native ``ExplicitBitVect`` binary constructor.
        A pickle fallback is retained so databases created by older builds of
        OpenTargetAI keep working.
        """
        if blob is None:
            return None
        try:
            return DataStructs.ExplicitBitVect(bytes(blob))
        except Exception:  # noqa: BLE001 - try legacy formats
            try:
                import pickle

                obj = pickle.loads(blob)
                if isinstance(obj, DataStructs.ExplicitBitVect):
                    return obj
                if isinstance(obj, str):
                    return DataStructs.CreateFromBitString(obj)
            except Exception as exc:  # noqa: BLE001
                logger.debug("FP deserialisation error: %s", exc)
        return None

    # ── Batch processing ──────────────────────────────────────────────────

    def batch_generate(
        self, smiles_list: list[str], fp_type: str = "ecfp4"
    ) -> list[Optional[DataStructs.ExplicitBitVect]]:
        """
        Generate fingerprint objects for a batch of SMILES strings.
        Returns a list aligned with the input (``None`` where invalid).
        """
        results: list[Optional[DataStructs.ExplicitBitVect]] = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles) if smiles else None
            results.append(self.generate_fp_object(mol, fp_type) if mol else None)
        return results

    # ── Similarity calculation ────────────────────────────────────────────

    @staticmethod
    def tanimoto(
        fp1: DataStructs.ExplicitBitVect, fp2: DataStructs.ExplicitBitVect
    ) -> float:
        """Tanimoto (Jaccard) similarity."""
        return DataStructs.TanimotoSimilarity(fp1, fp2)

    @staticmethod
    def dice(
        fp1: DataStructs.ExplicitBitVect, fp2: DataStructs.ExplicitBitVect
    ) -> float:
        """Dice similarity."""
        return DataStructs.DiceSimilarity(fp1, fp2)

    @staticmethod
    def cosine(
        fp1: DataStructs.ExplicitBitVect, fp2: DataStructs.ExplicitBitVect
    ) -> float:
        """Cosine similarity computed from the dense bit arrays."""
        arr1 = np.zeros((fp1.GetNumBits(),), dtype=np.float32)
        arr2 = np.zeros((fp2.GetNumBits(),), dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fp1, arr1)
        DataStructs.ConvertToNumpyArray(fp2, arr2)
        denom = float(np.linalg.norm(arr1) * np.linalg.norm(arr2))
        if denom == 0.0:
            return 0.0
        return float(np.dot(arr1, arr2) / denom)

    def compute_similarity(
        self,
        fp1: DataStructs.ExplicitBitVect,
        fp2: DataStructs.ExplicitBitVect,
        metric: str = "tanimoto",
    ) -> float:
        """Dispatch a single pairwise similarity computation."""
        metric = metric.strip().lower()
        if metric == "tanimoto":
            return self.tanimoto(fp1, fp2)
        if metric == "dice":
            return self.dice(fp1, fp2)
        if metric == "cosine":
            return self.cosine(fp1, fp2)
        raise ValueError(f"Unknown similarity metric: {metric}")

    def bulk_similarity(
        self,
        query_fp: DataStructs.ExplicitBitVect,
        db_fps: list[DataStructs.ExplicitBitVect],
        metric: str = "tanimoto",
    ) -> list[float]:
        """
        Compute similarity between a query and many database fingerprints.
        Uses RDKit's optimised C++ bulk routines for Tanimoto/Dice and a
        vectorised NumPy implementation for cosine.
        """
        metric = metric.strip().lower()
        if not db_fps:
            return []
        if metric == "tanimoto":
            return list(DataStructs.BulkTanimotoSimilarity(query_fp, db_fps))
        if metric == "dice":
            return list(DataStructs.BulkDiceSimilarity(query_fp, db_fps))
        if metric == "cosine":
            return self._bulk_cosine(query_fp, db_fps)
        raise ValueError(f"Unknown similarity metric: {metric}")

    @staticmethod
    def _bulk_cosine(
        query_fp: DataStructs.ExplicitBitVect,
        db_fps: list[DataStructs.ExplicitBitVect],
    ) -> list[float]:
        """Vectorised cosine similarity over a list of bit vectors."""
        q = np.zeros((query_fp.GetNumBits(),), dtype=np.float32)
        DataStructs.ConvertToNumpyArray(query_fp, q)
        q_norm = float(np.linalg.norm(q))
        if q_norm == 0.0:
            return [0.0] * len(db_fps)
        out: list[float] = []
        for fp in db_fps:
            v = np.zeros((fp.GetNumBits(),), dtype=np.float32)
            DataStructs.ConvertToNumpyArray(fp, v)
            denom = q_norm * float(np.linalg.norm(v))
            out.append(float(np.dot(q, v) / denom) if denom else 0.0)
        return out

    # ── Explainability helpers ────────────────────────────────────────────

    def morgan_bit_info(
        self, mol: Chem.Mol, radius: int = 2
    ) -> tuple[DataStructs.ExplicitBitVect, dict[int, tuple]]:
        """
        Return a Morgan fingerprint together with its bit-info map mapping each
        on-bit to the ``(atom_idx, radius)`` environments that set it. Used by
        the similarity explainer to highlight matched substructures.
        """
        gen = self._morgan2 if radius == 2 else self._morgan3
        ao = rdFingerprintGenerator.AdditionalOutput()
        ao.AllocateBitInfoMap()
        fp = gen.GetFingerprint(mol, additionalOutput=ao)
        return fp, ao.GetBitInfoMap()

    def shared_morgan_atoms(
        self, query_mol: Chem.Mol, ref_mol: Chem.Mol, radius: int = 2
    ) -> list[int]:
        """
        Identify atoms in ``query_mol`` that participate in Morgan environments
        also present in ``ref_mol`` (i.e. shared on-bits). Returns a sorted list
        of query atom indices to highlight.
        """
        try:
            q_fp, q_info = self.morgan_bit_info(query_mol, radius)
            r_fp, _ = self.morgan_bit_info(ref_mol, radius)
            shared_bits = set(q_fp.GetOnBits()) & set(r_fp.GetOnBits())
            atoms: set[int] = set()
            for bit in shared_bits:
                for atom_idx, rad in q_info.get(bit, ()):  # type: ignore[union-attr]
                    atoms.add(atom_idx)
                    if rad > 0:
                        env = Chem.FindAtomEnvironmentOfRadiusN(
                            query_mol, rad, atom_idx
                        )
                        for bond_id in env:
                            bond = query_mol.GetBondWithIdx(bond_id)
                            atoms.add(bond.GetBeginAtomIdx())
                            atoms.add(bond.GetEndAtomIdx())
            return sorted(atoms)
        except Exception as exc:  # noqa: BLE001
            logger.debug("shared_morgan_atoms error: %s", exc)
            return []
