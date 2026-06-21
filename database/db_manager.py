"""
Database Manager for OpenTargetAI.
Handles SQLite database initialization, schema creation, and CRUD operations.
"""

import sqlite3
import logging
import pickle
import zlib
from pathlib import Path
from typing import Optional, Any
from contextlib import contextmanager

logger = logging.getLogger(__name__)


class DatabaseManager:
    """Manages the SQLite database for OpenTargetAI."""

    SCHEMA_VERSION = 1

    def __init__(self, db_path: str) -> None:
        self.db_path = db_path
        self._connection: Optional[sqlite3.Connection] = None
        logger.info(f"DatabaseManager initialized: {db_path}")

    @contextmanager
    def get_connection(self):
        """Context manager for database connections."""
        conn = sqlite3.connect(self.db_path, timeout=30)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA synchronous=NORMAL")
        conn.execute("PRAGMA cache_size=-64000")
        conn.execute("PRAGMA foreign_keys=ON")
        try:
            yield conn
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            conn.close()

    def initialize_schema(self) -> None:
        """Create all required tables and indexes."""
        with self.get_connection() as conn:
            self._create_tables(conn)
            self._create_indexes(conn)
            self._insert_schema_version(conn)
        logger.info("Database schema initialized.")

    def _create_tables(self, conn: sqlite3.Connection) -> None:
        """Create all database tables."""
        conn.executescript("""
            CREATE TABLE IF NOT EXISTS schema_version (
                version INTEGER PRIMARY KEY,
                applied_at TEXT DEFAULT CURRENT_TIMESTAMP
            );

            CREATE TABLE IF NOT EXISTS compounds (
                compound_id INTEGER PRIMARY KEY AUTOINCREMENT,
                chembl_id TEXT UNIQUE,
                bindingdb_id TEXT,
                smiles TEXT NOT NULL,
                inchi TEXT,
                inchikey TEXT UNIQUE,
                mol_weight REAL,
                logp REAL,
                tpsa REAL,
                hbd INTEGER,
                hba INTEGER,
                rotatable_bonds INTEGER,
                fingerprint_ecfp4 BLOB,
                fingerprint_ecfp6 BLOB,
                fingerprint_maccs BLOB,
                fingerprint_ap BLOB,
                fingerprint_topo BLOB,
                created_at TEXT DEFAULT CURRENT_TIMESTAMP
            );

            CREATE TABLE IF NOT EXISTS targets (
                target_id INTEGER PRIMARY KEY AUTOINCREMENT,
                uniprot TEXT UNIQUE,
                gene_symbol TEXT,
                target_name TEXT NOT NULL,
                organism TEXT,
                target_class TEXT,
                protein_family TEXT,
                description TEXT,
                disease_associations TEXT,
                created_at TEXT DEFAULT CURRENT_TIMESTAMP
            );

            CREATE TABLE IF NOT EXISTS activities (
                activity_id INTEGER PRIMARY KEY AUTOINCREMENT,
                compound_id INTEGER NOT NULL,
                target_id INTEGER NOT NULL,
                activity_type TEXT,
                activity_value REAL,
                units TEXT,
                pchembl_value REAL,
                assay_type TEXT,
                source TEXT,
                reference TEXT,
                FOREIGN KEY (compound_id) REFERENCES compounds(compound_id),
                FOREIGN KEY (target_id) REFERENCES targets(target_id)
            );

            CREATE TABLE IF NOT EXISTS import_log (
                import_id INTEGER PRIMARY KEY AUTOINCREMENT,
                source TEXT NOT NULL,
                file_path TEXT,
                records_imported INTEGER DEFAULT 0,
                status TEXT DEFAULT 'started',
                started_at TEXT DEFAULT CURRENT_TIMESTAMP,
                completed_at TEXT,
                error_message TEXT
            );

            CREATE TABLE IF NOT EXISTS prediction_cache (
                cache_id INTEGER PRIMARY KEY AUTOINCREMENT,
                query_inchikey TEXT NOT NULL,
                fingerprint_type TEXT NOT NULL,
                similarity_metric TEXT NOT NULL,
                threshold REAL NOT NULL,
                top_n INTEGER NOT NULL,
                organism TEXT NOT NULL,
                results BLOB,
                created_at TEXT DEFAULT CURRENT_TIMESTAMP
            );
        """)

    def _create_indexes(self, conn: sqlite3.Connection) -> None:
        """Create performance indexes."""
        conn.executescript("""
            CREATE INDEX IF NOT EXISTS idx_compounds_inchikey ON compounds(inchikey);
            CREATE INDEX IF NOT EXISTS idx_compounds_chembl_id ON compounds(chembl_id);
            CREATE INDEX IF NOT EXISTS idx_targets_uniprot ON targets(uniprot);
            CREATE INDEX IF NOT EXISTS idx_targets_gene ON targets(gene_symbol);
            CREATE INDEX IF NOT EXISTS idx_targets_organism ON targets(organism);
            CREATE INDEX IF NOT EXISTS idx_activities_compound ON activities(compound_id);
            CREATE INDEX IF NOT EXISTS idx_activities_target ON activities(target_id);
            CREATE INDEX IF NOT EXISTS idx_activities_type ON activities(activity_type);
            CREATE INDEX IF NOT EXISTS idx_prediction_cache ON prediction_cache(query_inchikey, fingerprint_type);
        """)

    def _insert_schema_version(self, conn: sqlite3.Connection) -> None:
        """Insert schema version if not present."""
        conn.execute(
            "INSERT OR IGNORE INTO schema_version (version) VALUES (?)",
            (self.SCHEMA_VERSION,)
        )

    # ── Compound Operations ──────────────────────────────────────────────────

    def insert_compound(self, data: dict) -> int:
        """Insert a compound and return its ID."""
        with self.get_connection() as conn:
            cursor = conn.execute("""
                INSERT OR IGNORE INTO compounds
                    (chembl_id, bindingdb_id, smiles, inchi, inchikey,
                     mol_weight, logp, tpsa, hbd, hba, rotatable_bonds,
                     fingerprint_ecfp4, fingerprint_ecfp6, fingerprint_maccs,
                     fingerprint_ap, fingerprint_topo)
                VALUES
                    (:chembl_id, :bindingdb_id, :smiles, :inchi, :inchikey,
                     :mol_weight, :logp, :tpsa, :hbd, :hba, :rotatable_bonds,
                     :fingerprint_ecfp4, :fingerprint_ecfp6, :fingerprint_maccs,
                     :fingerprint_ap, :fingerprint_topo)
            """, data)
            if cursor.lastrowid:
                return cursor.lastrowid
            row = conn.execute(
                "SELECT compound_id FROM compounds WHERE inchikey=?", (data.get("inchikey"),)
            ).fetchone()
            return row["compound_id"] if row else -1

    def insert_target(self, data: dict) -> int:
        """Insert a target and return its ID."""
        with self.get_connection() as conn:
            cursor = conn.execute("""
                INSERT OR IGNORE INTO targets
                    (uniprot, gene_symbol, target_name, organism, target_class,
                     protein_family, description, disease_associations)
                VALUES
                    (:uniprot, :gene_symbol, :target_name, :organism, :target_class,
                     :protein_family, :description, :disease_associations)
            """, data)
            if cursor.lastrowid:
                return cursor.lastrowid
            row = conn.execute(
                "SELECT target_id FROM targets WHERE uniprot=?", (data.get("uniprot"),)
            ).fetchone()
            return row["target_id"] if row else -1

    def insert_activity(self, data: dict) -> int:
        """Insert an activity record."""
        with self.get_connection() as conn:
            cursor = conn.execute("""
                INSERT INTO activities
                    (compound_id, target_id, activity_type, activity_value,
                     units, pchembl_value, assay_type, source, reference)
                VALUES
                    (:compound_id, :target_id, :activity_type, :activity_value,
                     :units, :pchembl_value, :assay_type, :source, :reference)
            """, data)
            return cursor.lastrowid or -1

    def bulk_insert_compounds(self, records: list[dict]) -> int:
        """Bulk insert compounds, return count inserted."""
        if not records:
            return 0
        keys = records[0].keys()
        placeholders = ", ".join("?" * len(keys))
        cols = ", ".join(keys)
        query = f"INSERT OR IGNORE INTO compounds ({cols}) VALUES ({placeholders})"
        with self.get_connection() as conn:
            cursor = conn.executemany(query, [list(r.values()) for r in records])
            return cursor.rowcount

    def bulk_insert_activities(self, records: list[dict]) -> int:
        """Bulk insert activities, return count inserted."""
        if not records:
            return 0
        with self.get_connection() as conn:
            cursor = conn.executemany("""
                INSERT OR IGNORE INTO activities
                    (compound_id, target_id, activity_type, activity_value,
                     units, pchembl_value, assay_type, source, reference)
                VALUES (?,?,?,?,?,?,?,?,?)
            """, [
                (r.get("compound_id"), r.get("target_id"), r.get("activity_type"),
                 r.get("activity_value"), r.get("units"), r.get("pchembl_value"),
                 r.get("assay_type"), r.get("source"), r.get("reference"))
                for r in records
            ])
            return cursor.rowcount

    # ── Fingerprint Retrieval ─────────────────────────────────────────────

    def get_all_fingerprints(self, fp_type: str = "ecfp4") -> list[dict]:
        """
        Return all compounds with their specified fingerprint.
        fp_type: 'ecfp4', 'ecfp6', 'maccs', 'ap', 'topo'
        """
        col = f"fingerprint_{fp_type}"
        with self.get_connection() as conn:
            rows = conn.execute(
                f"SELECT compound_id, {col} FROM compounds WHERE {col} IS NOT NULL"
            ).fetchall()
        return [{"compound_id": r["compound_id"], "fingerprint": r[col]} for r in rows]

    def get_compound_targets(self, compound_ids: list[int]) -> list[dict]:
        """Get all targets associated with a list of compound IDs."""
        if not compound_ids:
            return []
        placeholders = ",".join("?" * len(compound_ids))
        with self.get_connection() as conn:
            rows = conn.execute(f"""
                SELECT
                    a.compound_id, a.activity_type, a.activity_value,
                    a.pchembl_value, a.units, a.source,
                    t.target_id, t.uniprot, t.gene_symbol, t.target_name,
                    t.organism, t.target_class, t.protein_family
                FROM activities a
                JOIN targets t ON a.target_id = t.target_id
                WHERE a.compound_id IN ({placeholders})
            """, compound_ids).fetchall()
        return [dict(r) for r in rows]

    def get_target_by_id(self, target_id: int) -> Optional[dict]:
        """Get full target details by ID."""
        with self.get_connection() as conn:
            row = conn.execute(
                "SELECT * FROM targets WHERE target_id=?", (target_id,)
            ).fetchone()
        return dict(row) if row else None

    def get_compound_by_id(self, compound_id: int) -> Optional[dict]:
        """Get compound details by ID."""
        with self.get_connection() as conn:
            row = conn.execute(
                "SELECT * FROM compounds WHERE compound_id=?", (compound_id,)
            ).fetchone()
        return dict(row) if row else None

    def search_compounds_by_smiles(self, smiles: str) -> Optional[dict]:
        """Find a compound by its SMILES."""
        with self.get_connection() as conn:
            row = conn.execute(
                "SELECT * FROM compounds WHERE smiles=?", (smiles,)
            ).fetchone()
        return dict(row) if row else None

    # ── Statistics ────────────────────────────────────────────────────────

    def get_database_stats(self) -> dict:
        """Return statistics about the database contents."""
        with self.get_connection() as conn:
            n_compounds = conn.execute("SELECT COUNT(*) FROM compounds").fetchone()[0]
            n_targets = conn.execute("SELECT COUNT(*) FROM targets").fetchone()[0]
            n_activities = conn.execute("SELECT COUNT(*) FROM activities").fetchone()[0]
            organisms = conn.execute(
                "SELECT organism, COUNT(*) as cnt FROM targets GROUP BY organism ORDER BY cnt DESC LIMIT 10"
            ).fetchall()
        return {
            "compounds": n_compounds,
            "targets": n_targets,
            "activities": n_activities,
            "organisms": [dict(r) for r in organisms],
        }

    # ── Import Log ────────────────────────────────────────────────────────

    def log_import_start(self, source: str, file_path: str) -> int:
        """Start an import log entry, return import_id."""
        with self.get_connection() as conn:
            cursor = conn.execute(
                "INSERT INTO import_log (source, file_path, status) VALUES (?,?,?)",
                (source, file_path, "started")
            )
            return cursor.lastrowid or -1

    def log_import_complete(self, import_id: int, records: int) -> None:
        """Mark import as complete."""
        with self.get_connection() as conn:
            conn.execute("""
                UPDATE import_log
                SET status='completed', records_imported=?, completed_at=CURRENT_TIMESTAMP
                WHERE import_id=?
            """, (records, import_id))

    def log_import_error(self, import_id: int, error: str) -> None:
        """Mark import as failed."""
        with self.get_connection() as conn:
            conn.execute("""
                UPDATE import_log
                SET status='failed', error_message=?, completed_at=CURRENT_TIMESTAMP
                WHERE import_id=?
            """, (error, import_id))

    # ── Prediction Cache ──────────────────────────────────────────────────

    def get_cached_prediction(
        self, inchikey: str, fp_type: str, metric: str,
        threshold: float, top_n: int, organism: str
    ) -> Optional[Any]:
        """Return cached prediction results or None."""
        with self.get_connection() as conn:
            row = conn.execute("""
                SELECT results FROM prediction_cache
                WHERE query_inchikey=? AND fingerprint_type=?
                  AND similarity_metric=? AND threshold=?
                  AND top_n=? AND organism=?
                ORDER BY created_at DESC LIMIT 1
            """, (inchikey, fp_type, metric, threshold, top_n, organism)).fetchone()
        if row and row["results"]:
            return pickle.loads(zlib.decompress(row["results"]))
        return None

    def save_prediction_cache(
        self, inchikey: str, fp_type: str, metric: str,
        threshold: float, top_n: int, organism: str, results: Any
    ) -> None:
        """Save prediction results to cache."""
        blob = zlib.compress(pickle.dumps(results))
        with self.get_connection() as conn:
            conn.execute("""
                INSERT INTO prediction_cache
                    (query_inchikey, fingerprint_type, similarity_metric,
                     threshold, top_n, organism, results)
                VALUES (?,?,?,?,?,?,?)
            """, (inchikey, fp_type, metric, threshold, top_n, organism, blob))

    def clear_prediction_cache(self) -> None:
        """Clear all cached predictions."""
        with self.get_connection() as conn:
            conn.execute("DELETE FROM prediction_cache")
