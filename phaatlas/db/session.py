from __future__ import annotations

import os
from contextlib import contextmanager
from pathlib import Path

from sqlalchemy import create_engine, event
from sqlalchemy.orm import Session, sessionmaker

from phaatlas.config_loader import REPO_ROOT
from phaatlas.db.models import Base
from phaatlas.db.views import PROTEIN_MASTER_EXPORT_VIEW_SQL

DEFAULT_DB_PATH = REPO_ROOT / "PHA_reference" / "pha_reference.sqlite"


def get_db_path() -> Path:
    raw = os.environ.get("PHA_REFERENCE_DATABASE_PATH")
    return Path(raw) if raw else DEFAULT_DB_PATH


def get_engine(db_path: Path | None = None):
    path = db_path or get_db_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    engine = create_engine(f"sqlite:///{path}", future=True)

    # SQLite's default foreign-key enforcement is OFF per-connection --
    # without this, an orphaned family_assignment/source_evidence/phenotype
    # row (e.g. from a bug in ingest code) would insert silently instead of
    # failing fast.
    @event.listens_for(engine, "connect")
    def _set_sqlite_pragma(dbapi_connection, connection_record):
        cursor = dbapi_connection.cursor()
        cursor.execute("PRAGMA foreign_keys=ON")
        cursor.close()

    return engine


# Additive-only, idempotent column migrations. `Base.metadata.create_all`
# only creates tables that don't exist yet -- it never adds a column to an
# already-existing table, so a database created before a model gained a new
# column (e.g. cluster95_id, added after protein/family_assignment were
# already populated) would silently stay on the old schema forever. No
# alembic here (Phase 1 deliberately kept this lightweight); this list is
# the whole migration story, and it's safe to re-run unconditionally --
# each entry is skipped once the column already exists.
_COLUMN_MIGRATIONS: list[tuple[str, str, str]] = [
    # (table, column, SQL type+constraints for ADD COLUMN)
    ("family_assignment", "cluster95_id", "VARCHAR"),
    ("family_assignment", "cluster95_representative", "VARCHAR"),
    ("family_assignment", "cluster95_size", "INTEGER"),
]


def _run_column_migrations(conn) -> None:
    for table, column, coltype in _COLUMN_MIGRATIONS:
        existing = {row[1] for row in conn.exec_driver_sql(f"PRAGMA table_info({table})").fetchall()}
        if column not in existing:
            conn.exec_driver_sql(f"ALTER TABLE {table} ADD COLUMN {column} {coltype}")


def init_db(db_path: Path | None = None) -> None:
    engine = get_engine(db_path)
    Base.metadata.create_all(engine)
    with engine.begin() as conn:
        _run_column_migrations(conn)
        conn.exec_driver_sql("DROP VIEW IF EXISTS protein_master_export")
        conn.exec_driver_sql(PROTEIN_MASTER_EXPORT_VIEW_SQL)


def get_sessionmaker(db_path: Path | None = None) -> sessionmaker:
    # Deliberately NOT cached globally: tests (and any future multi-database
    # use, e.g. a scratch DB for a dry run) point this at different paths
    # within the same process, and a cached sessionmaker would silently
    # keep binding to whichever db_path was passed first.
    engine = get_engine(db_path)
    return sessionmaker(bind=engine, expire_on_commit=False, future=True)


@contextmanager
def session_scope(db_path: Path | None = None):
    Session_ = get_sessionmaker(db_path)
    session: Session = Session_()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()
