-- ============================================================================
-- OpenTargetAI — SQLite database schema
-- ----------------------------------------------------------------------------
-- This DDL mirrors the schema created programmatically by
-- database/db_manager.py (DatabaseManager.initialize_schema). It is provided
-- as documentation and for manual / external database creation. The
-- application creates and migrates the database automatically on first launch,
-- so running this file by hand is optional.
--
-- Apply with:   sqlite3 data/opentargetai.db < schema.sql
-- ============================================================================

PRAGMA journal_mode = WAL;
PRAGMA foreign_keys = ON;

-- ── Schema version ──────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS schema_version (
    version    INTEGER PRIMARY KEY,
    applied_at TEXT DEFAULT CURRENT_TIMESTAMP
);

-- ── Compounds ───────────────────────────────────────────────────────────────
-- One row per unique compound. Fingerprints are stored as compact native
-- RDKit ExplicitBitVect binary blobs (ToBinary()).
CREATE TABLE IF NOT EXISTS compounds (
    compound_id        INTEGER PRIMARY KEY AUTOINCREMENT,
    chembl_id          TEXT UNIQUE,
    bindingdb_id       TEXT,
    smiles             TEXT NOT NULL,
    inchi              TEXT,
    inchikey           TEXT UNIQUE,
    mol_weight         REAL,
    logp               REAL,
    tpsa               REAL,
    hbd                INTEGER,
    hba                INTEGER,
    rotatable_bonds    INTEGER,
    fingerprint_ecfp4  BLOB,   -- Morgan radius 2, 2048 bits
    fingerprint_ecfp6  BLOB,   -- Morgan radius 3, 2048 bits
    fingerprint_maccs  BLOB,   -- MACCS keys, 167 bits
    fingerprint_ap     BLOB,   -- Atom pair, folded 2048 bits
    fingerprint_topo   BLOB,   -- RDKit topological, 2048 bits
    created_at         TEXT DEFAULT CURRENT_TIMESTAMP
);

-- ── Targets ─────────────────────────────────────────────────────────────────
-- One row per protein target (keyed by UniProt accession where available).
CREATE TABLE IF NOT EXISTS targets (
    target_id            INTEGER PRIMARY KEY AUTOINCREMENT,
    uniprot              TEXT UNIQUE,
    gene_symbol          TEXT,
    target_name          TEXT NOT NULL,
    organism             TEXT,
    target_class         TEXT,
    protein_family       TEXT,
    description          TEXT,
    disease_associations TEXT,
    created_at           TEXT DEFAULT CURRENT_TIMESTAMP
);

-- ── Activities ──────────────────────────────────────────────────────────────
-- Measured bioactivities linking compounds to targets.
CREATE TABLE IF NOT EXISTS activities (
    activity_id     INTEGER PRIMARY KEY AUTOINCREMENT,
    compound_id     INTEGER NOT NULL,
    target_id       INTEGER NOT NULL,
    activity_type   TEXT,           -- IC50, Ki, Kd, EC50, …
    activity_value  REAL,
    units           TEXT,
    pchembl_value   REAL,           -- -log10(value in M)
    assay_type      TEXT,
    source          TEXT,           -- ChEMBL / BindingDB / Demo
    reference       TEXT,
    FOREIGN KEY (compound_id) REFERENCES compounds (compound_id),
    FOREIGN KEY (target_id)   REFERENCES targets (target_id)
);

-- ── Import log ──────────────────────────────────────────────────────────────
-- Tracks import runs (used for resume / auditing).
CREATE TABLE IF NOT EXISTS import_log (
    import_id         INTEGER PRIMARY KEY AUTOINCREMENT,
    source            TEXT NOT NULL,
    file_path         TEXT,
    records_imported  INTEGER DEFAULT 0,
    status            TEXT DEFAULT 'started',  -- started / completed / failed
    started_at        TEXT DEFAULT CURRENT_TIMESTAMP,
    completed_at      TEXT,
    error_message     TEXT
);

-- ── Prediction cache ────────────────────────────────────────────────────────
-- Caches prediction results keyed by query + parameters (compressed blob).
CREATE TABLE IF NOT EXISTS prediction_cache (
    cache_id          INTEGER PRIMARY KEY AUTOINCREMENT,
    query_inchikey    TEXT NOT NULL,
    fingerprint_type  TEXT NOT NULL,
    similarity_metric TEXT NOT NULL,
    threshold         REAL NOT NULL,
    top_n             INTEGER NOT NULL,
    organism          TEXT NOT NULL,
    results           BLOB,
    created_at        TEXT DEFAULT CURRENT_TIMESTAMP
);

-- ── Indexes ─────────────────────────────────────────────────────────────────
CREATE INDEX IF NOT EXISTS idx_compounds_inchikey  ON compounds (inchikey);
CREATE INDEX IF NOT EXISTS idx_compounds_chembl_id ON compounds (chembl_id);
CREATE INDEX IF NOT EXISTS idx_targets_uniprot     ON targets (uniprot);
CREATE INDEX IF NOT EXISTS idx_targets_gene        ON targets (gene_symbol);
CREATE INDEX IF NOT EXISTS idx_targets_organism    ON targets (organism);
CREATE INDEX IF NOT EXISTS idx_activities_compound ON activities (compound_id);
CREATE INDEX IF NOT EXISTS idx_activities_target   ON activities (target_id);
CREATE INDEX IF NOT EXISTS idx_activities_type     ON activities (activity_type);
CREATE INDEX IF NOT EXISTS idx_prediction_cache
    ON prediction_cache (query_inchikey, fingerprint_type);

INSERT OR IGNORE INTO schema_version (version) VALUES (1);
