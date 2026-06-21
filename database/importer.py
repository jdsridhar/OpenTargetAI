"""
Database Importer for OpenTargetAI.
Handles importing data from ChEMBL and BindingDB exports into local SQLite.
"""

import csv
import gzip
import json
import logging
import sqlite3
import time
from pathlib import Path
from typing import Callable, Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from database.db_manager import DatabaseManager
from similarity.fingerprint_engine import FingerprintEngine

logger = logging.getLogger(__name__)

ProgressCallback = Callable[[int, int, str], None]


class DataImporter:
    """
    Imports ChEMBL and BindingDB data into the local database.
    Supports resume on interrupted imports.
    """

    BATCH_SIZE = 500

    def __init__(
        self,
        db_manager: DatabaseManager,
        progress_callback: Optional[ProgressCallback] = None,
    ) -> None:
        self.db = db_manager
        self.fp_engine = FingerprintEngine()
        self.progress_callback = progress_callback or (lambda cur, tot, msg: None)

    # ── ChEMBL Import ────────────────────────────────────────────────────

    def import_chembl_activities_tsv(self, file_path: str) -> int:
        """
        Import a ChEMBL activities TSV export.
        Expected columns: molregno, canonical_smiles, chembl_id,
            target_chembl_id, uniprot_id, gene_name, target_name,
            organism, standard_type, standard_value, standard_units,
            pchembl_value, assay_type
        Returns total records imported.
        """
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        import_id = self.db.log_import_start("chembl", file_path)

        try:
            open_fn = gzip.open if path.suffix == ".gz" else open
            mode = "rt" if path.suffix == ".gz" else "r"

            with open_fn(path, mode, encoding="utf-8", errors="replace") as f:
                reader = csv.DictReader(f, delimiter="\t")
                rows = list(reader)

            total = len(rows)
            imported = 0
            skipped = 0

            for batch_start in range(0, total, self.BATCH_SIZE):
                batch = rows[batch_start: batch_start + self.BATCH_SIZE]
                for row in batch:
                    try:
                        compound_id = self._process_chembl_compound(row)
                        if compound_id < 0:
                            skipped += 1
                            continue

                        target_id = self._process_chembl_target(row)
                        if target_id < 0:
                            skipped += 1
                            continue

                        self.db.insert_activity({
                            "compound_id": compound_id,
                            "target_id": target_id,
                            "activity_type": row.get("standard_type", ""),
                            "activity_value": self._safe_float(row.get("standard_value")),
                            "units": row.get("standard_units", ""),
                            "pchembl_value": self._safe_float(row.get("pchembl_value")),
                            "assay_type": row.get("assay_type", ""),
                            "source": "ChEMBL",
                            "reference": row.get("target_chembl_id", ""),
                        })
                        imported += 1
                    except Exception as e:
                        logger.debug(f"Row error: {e}")
                        skipped += 1

                progress_pct = batch_start + len(batch)
                self.progress_callback(
                    progress_pct, total,
                    f"ChEMBL: {imported} imported, {skipped} skipped"
                )

            self.db.log_import_complete(import_id, imported)
            logger.info(f"ChEMBL import complete: {imported} records")
            return imported

        except Exception as e:
            self.db.log_import_error(import_id, str(e))
            logger.exception("ChEMBL import failed")
            raise

    def _process_chembl_compound(self, row: dict) -> int:
        """Process a single ChEMBL compound row, return compound_id."""
        smiles = row.get("canonical_smiles", "").strip()
        if not smiles:
            return -1

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return -1

        record = self._build_compound_record(
            mol, smiles, chembl_id=row.get("chembl_id", "") or None
        )
        return self.db.insert_compound(record)

    def _process_chembl_target(self, row: dict) -> int:
        """Process a single ChEMBL target row, return target_id."""
        uniprot = row.get("uniprot_id", "").strip()
        target_name = row.get("target_name", "").strip()
        if not target_name:
            return -1

        return self.db.insert_target({
            "uniprot": uniprot if uniprot else f"NOUNIPROT_{row.get('target_chembl_id', 'UNKNOWN')}",
            "gene_symbol": row.get("gene_name", ""),
            "target_name": target_name,
            "organism": row.get("organism", ""),
            "target_class": row.get("target_type", ""),
            "protein_family": "",
            "description": "",
            "disease_associations": "",
        })

    # ── BindingDB Import ─────────────────────────────────────────────────

    def import_bindingdb_tsv(self, file_path: str) -> int:
        """
        Import a BindingDB TSV export.
        Common columns: Ligand SMILES, Target Name, UniProt (SwissProt) Primary ID
            of Target Chain, Ki (nM), IC50 (nM), Kd (nM), EC50 (nM),
            Target Source Organism According to Curator or DataSource
        Returns total records imported.
        """
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        import_id = self.db.log_import_start("bindingdb", file_path)

        try:
            open_fn = gzip.open if path.suffix == ".gz" else open
            mode = "rt" if path.suffix == ".gz" else "r"

            df = pd.read_csv(path, sep="\t", low_memory=False, on_bad_lines="skip",
                             encoding="utf-8", encoding_errors="replace")
            total = len(df)
            imported = 0
            skipped = 0

            activity_cols = {
                "Ki (nM)": "Ki", "IC50 (nM)": "IC50",
                "Kd (nM)": "Kd", "EC50 (nM)": "EC50",
            }

            for batch_start in range(0, total, self.BATCH_SIZE):
                batch_df = df.iloc[batch_start: batch_start + self.BATCH_SIZE]
                for _, row in batch_df.iterrows():
                    try:
                        smiles = str(row.get("Ligand SMILES", "")).strip()
                        if not smiles or smiles == "nan":
                            skipped += 1
                            continue

                        mol = Chem.MolFromSmiles(smiles)
                        if mol is None:
                            skipped += 1
                            continue

                        record = self._build_compound_record(
                            mol,
                            smiles,
                            bindingdb_id=str(row.get("BindingDB MonomerID", "")) or None,
                        )
                        compound_id = self.db.insert_compound(record)
                        if compound_id < 0:
                            skipped += 1
                            continue

                        uniprot = str(row.get("UniProt (SwissProt) Primary ID of Target Chain", "")).strip()
                        target_name = str(row.get("Target Name", "")).strip()
                        if not target_name or target_name == "nan":
                            skipped += 1
                            continue

                        target_id = self.db.insert_target({
                            "uniprot": uniprot if uniprot and uniprot != "nan" else f"BDB_{hash(target_name) % 999999}",
                            "gene_symbol": str(row.get("UniProt (SwissProt) Entry Name of Target Chain", "")).split("_")[0],
                            "target_name": target_name,
                            "organism": str(row.get("Target Source Organism According to Curator or DataSource", "")),
                            "target_class": "",
                            "protein_family": "",
                            "description": "",
                            "disease_associations": "",
                        })

                        for col, atype in activity_cols.items():
                            val = self._safe_float(row.get(col))
                            if val is not None and val > 0:
                                self.db.insert_activity({
                                    "compound_id": compound_id,
                                    "target_id": target_id,
                                    "activity_type": atype,
                                    "activity_value": val,
                                    "units": "nM",
                                    "pchembl_value": self._nm_to_pchembl(val),
                                    "assay_type": "B",
                                    "source": "BindingDB",
                                    "reference": "",
                                })
                                imported += 1
                                break  # one activity per row is enough

                    except Exception as e:
                        logger.debug(f"BindingDB row error: {e}")
                        skipped += 1

                self.progress_callback(
                    batch_start + len(batch_df), total,
                    f"BindingDB: {imported} imported, {skipped} skipped"
                )

            self.db.log_import_complete(import_id, imported)
            logger.info(f"BindingDB import complete: {imported} records")
            return imported

        except Exception as e:
            self.db.log_import_error(import_id, str(e))
            logger.exception("BindingDB import failed")
            raise

    # ── Demo Data ────────────────────────────────────────────────────────

    def import_demo_data(self) -> int:
        """
        Import the curated demo dataset (see :mod:`database.demo_data`) so the
        application works out of the box without a multi-gigabyte ChEMBL or
        BindingDB download. Returns the number of activity records imported.
        """
        from database import demo_data

        dataset = demo_data.get_demo_dataset()
        logger.info("Importing demo dataset...")
        import_id = self.db.log_import_start("demo", "<built-in>")

        # Targets keyed by gene symbol so activities can reference them by name.
        target_ids: dict[str, int] = {}
        for uniprot, gene, name, organism, tclass, desc, diseases in dataset["targets"]:
            target_ids[gene] = self.db.insert_target(
                {
                    "uniprot": uniprot,
                    "gene_symbol": gene,
                    "target_name": name,
                    "organism": organism,
                    "target_class": tclass,
                    "protein_family": "",
                    "description": desc,
                    "disease_associations": diseases,
                }
            )

        # Compounds keyed by name.
        compound_ids: dict[str, int] = {}
        total = len(dataset["compounds"])
        for i, (name, smiles, chembl_id) in enumerate(dataset["compounds"]):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning("Demo compound %s has invalid SMILES, skipped", name)
                continue
            record = self._build_compound_record(
                mol, smiles, chembl_id=chembl_id, bindingdb_id=None
            )
            compound_ids[name] = self.db.insert_compound(record)
            self.progress_callback(i + 1, total, f"Demo compounds: {i + 1}/{total}")

        imported = 0
        for name, gene, atype, val, units, pchembl in dataset["activities"]:
            c_id = compound_ids.get(name, -1)
            t_id = target_ids.get(gene, -1)
            if c_id > 0 and t_id > 0:
                self.db.insert_activity(
                    {
                        "compound_id": c_id,
                        "target_id": t_id,
                        "activity_type": atype,
                        "activity_value": val,
                        "units": units,
                        "pchembl_value": pchembl,
                        "assay_type": "B",
                        "source": "Demo",
                        "reference": "",
                    }
                )
                imported += 1

        self.db.log_import_complete(import_id, imported)
        logger.info("Demo data imported: %d activity records", imported)
        return imported

    # ── Helpers ───────────────────────────────────────────────────────────

    def _build_compound_record(
        self,
        mol: "Chem.Mol",
        smiles: str,
        chembl_id: Optional[str] = None,
        bindingdb_id: Optional[str] = None,
    ) -> dict:
        """Build a complete compound row (descriptors + all fingerprints)."""
        fps = self.fp_engine.generate_all_fingerprints(mol)
        try:
            inchi = Chem.MolToInchi(mol)
            inchikey = Chem.MolToInchiKey(mol)
        except Exception:  # noqa: BLE001 - some structures break InChI
            inchi = None
            inchikey = None
        return {
            "chembl_id": chembl_id,
            "bindingdb_id": bindingdb_id,
            "smiles": smiles,
            "inchi": inchi,
            "inchikey": inchikey,
            "mol_weight": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "tpsa": rdMolDescriptors.CalcTPSA(mol),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "fingerprint_ecfp4": fps.get("ecfp4"),
            "fingerprint_ecfp6": fps.get("ecfp6"),
            "fingerprint_maccs": fps.get("maccs"),
            "fingerprint_ap": fps.get("ap"),
            "fingerprint_topo": fps.get("topo"),
        }

    @staticmethod
    def _safe_float(val) -> Optional[float]:
        """Safely convert to float."""
        try:
            f = float(val)
            return f if f == f else None  # NaN check
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _nm_to_pchembl(nm_value: float) -> Optional[float]:
        """Convert nM activity value to pChEMBL (pIC50/pKi etc.)."""
        if nm_value and nm_value > 0:
            import math
            return -math.log10(nm_value * 1e-9)
        return None
