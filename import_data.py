#!/usr/bin/env python3
"""
Standalone command-line database importer for OpenTargetAI.

Builds / populates the local SQLite database used by the desktop application
without launching the GUI. Supports the curated demo dataset and full ChEMBL /
BindingDB tab-separated exports, generating fingerprints and indexes on the way.

Examples
--------
    # Create the DB and load the built-in demo dataset
    python import_data.py --demo

    # Import a ChEMBL activities TSV export into a specific database
    python import_data.py --chembl chembl_activities.tsv --db data/opentargetai.db

    # Import a BindingDB TSV export
    python import_data.py --bindingdb BindingDB_All.tsv

    # Show database statistics
    python import_data.py --stats
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from database.db_manager import DatabaseManager  # noqa: E402
from database.importer import DataImporter  # noqa: E402

DEFAULT_DB = PROJECT_ROOT / "data" / "opentargetai.db"


def _setup_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )


class _ConsoleProgress:
    """A minimal throttled console progress bar."""

    def __init__(self) -> None:
        self._last = 0.0

    def __call__(self, current: int, total: int, message: str) -> None:
        now = time.time()
        if now - self._last < 0.1 and current < total:
            return
        self._last = now
        total = max(1, total)
        pct = min(100, int(current / total * 100))
        bar = "#" * (pct // 2) + "-" * (50 - pct // 2)
        try:
            sys.stdout.write(f"\r[{bar}] {pct:3d}%  {message[:60]:<60}")
            sys.stdout.flush()
            if current >= total:
                sys.stdout.write("\n")
        except (UnicodeEncodeError, OSError):
            # Never let progress reporting abort an import.
            pass


def _print_stats(db: DatabaseManager) -> None:
    stats = db.get_database_stats()
    print("\nDatabase statistics")
    print("-------------------")
    print(f"  Compounds:   {stats['compounds']:,}")
    print(f"  Targets:     {stats['targets']:,}")
    print(f"  Activities:  {stats['activities']:,}")
    if stats["organisms"]:
        print("  Top organisms:")
        for org in stats["organisms"]:
            print(f"    - {org['organism'] or 'Unknown'}: {org['cnt']:,}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="OpenTargetAI standalone database importer",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--db",
        default=str(DEFAULT_DB),
        help=f"Path to the SQLite database (default: {DEFAULT_DB})",
    )
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--demo", action="store_true", help="Import the curated demo dataset"
    )
    source.add_argument("--chembl", metavar="TSV", help="Import a ChEMBL TSV export")
    source.add_argument(
        "--bindingdb", metavar="TSV", help="Import a BindingDB TSV export"
    )
    parser.add_argument(
        "--stats", action="store_true", help="Print database statistics and exit"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose (debug) logging"
    )
    args = parser.parse_args(argv)

    # Prefer UTF-8 console output where the platform allows it.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")  # type: ignore[attr-defined]
        except Exception:  # noqa: BLE001
            pass

    _setup_logging(args.verbose)

    db_path = Path(args.db)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    db = DatabaseManager(str(db_path))
    db.initialize_schema()
    print(f"Database: {db_path}")

    if args.stats and not (args.demo or args.chembl or args.bindingdb):
        _print_stats(db)
        return 0

    importer = DataImporter(db, progress_callback=_ConsoleProgress())

    try:
        if args.demo:
            print("Importing demo dataset…")
            count = importer.import_demo_data()
        elif args.chembl:
            print(f"Importing ChEMBL export: {args.chembl}")
            count = importer.import_chembl_activities_tsv(args.chembl)
        elif args.bindingdb:
            print(f"Importing BindingDB export: {args.bindingdb}")
            count = importer.import_bindingdb_tsv(args.bindingdb)
        else:
            parser.print_help()
            return 1
    except FileNotFoundError as exc:
        print(f"\nError: file not found — {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # noqa: BLE001
        print(f"\nImport failed: {exc}", file=sys.stderr)
        return 3

    print(f"\nDone. {count:,} activity records imported.")
    _print_stats(db)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
