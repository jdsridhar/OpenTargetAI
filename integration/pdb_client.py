"""
PDB structure retrieval for OpenTargetAI.

Looks up experimentally determined 3D structures for a predicted target via its
UniProt accession, using the public PDBe "best structures" API (a simple GET
endpoint that returns resolution and experimental method), and downloads the
coordinate files from RCSB.

All network access is defensive: failures (offline, rate-limited, unknown
accession) raise :class:`PDBLookupError` with a readable message rather than
crashing the caller. Calls should be made from a background worker thread.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

PDBE_BEST_STRUCTURES = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{acc}"
RCSB_DOWNLOAD = "https://files.rcsb.org/download/{pdb_id}.{fmt}"
DEFAULT_TIMEOUT = 20


class PDBLookupError(RuntimeError):
    """Raised when structure lookup or download fails."""


@dataclass
class PDBStructure:
    """A single experimentally determined structure mapped to a UniProt target."""

    pdb_id: str
    resolution: Optional[float]
    experimental_method: str
    chain_id: str = ""
    coverage: Optional[float] = None
    uniprot: str = ""

    @property
    def resolution_str(self) -> str:
        return f"{self.resolution:.2f} Å" if self.resolution else "N/A"


class PDBClient:
    """Thin client over the PDBe / RCSB public REST APIs."""

    def __init__(self, timeout: int = DEFAULT_TIMEOUT) -> None:
        self.timeout = timeout

    def get_structures_for_uniprot(
        self, uniprot: str, max_results: int = 25
    ) -> list[PDBStructure]:
        """
        Return experimental structures for a UniProt accession, best
        (highest-resolution / coverage) first.
        """
        if not uniprot or uniprot.startswith(("NOUNIPROT", "BDB_")):
            raise PDBLookupError(
                "Target has no resolvable UniProt accession for structure lookup."
            )
        try:
            import requests
        except ImportError as exc:  # pragma: no cover
            raise PDBLookupError("The 'requests' package is required.") from exc

        url = PDBE_BEST_STRUCTURES.format(acc=uniprot.strip())
        try:
            resp = requests.get(url, timeout=self.timeout)
        except Exception as exc:  # noqa: BLE001 - network error
            raise PDBLookupError(f"Network error contacting PDBe: {exc}") from exc

        if resp.status_code == 404:
            return []
        if resp.status_code != 200:
            raise PDBLookupError(
                f"PDBe returned HTTP {resp.status_code} for {uniprot}."
            )

        try:
            data = resp.json()
        except Exception as exc:  # noqa: BLE001
            raise PDBLookupError("Could not parse PDBe response.") from exc

        entries = data.get(uniprot.strip(), [])
        structures: list[PDBStructure] = []
        seen: set[str] = set()
        for entry in entries:
            pdb_id = str(entry.get("pdb_id", "")).upper()
            if not pdb_id or pdb_id in seen:
                continue
            seen.add(pdb_id)
            structures.append(
                PDBStructure(
                    pdb_id=pdb_id,
                    resolution=_safe_float(entry.get("resolution")),
                    experimental_method=entry.get("experimental_method")
                    or "Unknown",
                    chain_id=entry.get("chain_id", ""),
                    coverage=_safe_float(entry.get("coverage")),
                    uniprot=uniprot.strip(),
                )
            )
            if len(structures) >= max_results:
                break
        return structures

    def download_structure(
        self, pdb_id: str, dest_dir: str, fmt: str = "pdb"
    ) -> str:
        """Download a structure file from RCSB into ``dest_dir``; return its path."""
        try:
            import requests
        except ImportError as exc:  # pragma: no cover
            raise PDBLookupError("The 'requests' package is required.") from exc

        pdb_id = pdb_id.strip().lower()
        dest = Path(dest_dir)
        dest.mkdir(parents=True, exist_ok=True)
        out_path = dest / f"{pdb_id}.{fmt}"
        url = RCSB_DOWNLOAD.format(pdb_id=pdb_id, fmt=fmt)
        try:
            resp = requests.get(url, timeout=self.timeout)
        except Exception as exc:  # noqa: BLE001
            raise PDBLookupError(f"Network error downloading {pdb_id}: {exc}") from exc
        if resp.status_code != 200:
            raise PDBLookupError(
                f"Could not download {pdb_id} (HTTP {resp.status_code})."
            )
        out_path.write_bytes(resp.content)
        logger.info("Downloaded %s to %s", pdb_id, out_path)
        return str(out_path)


def _safe_float(value) -> Optional[float]:
    """Best-effort float conversion."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
