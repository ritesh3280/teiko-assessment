"""Opening the SQLite database built by ``load_data.py``."""

from __future__ import annotations

import sqlite3
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path


def connect(db_path: str | Path) -> sqlite3.Connection:
    """Open an existing database read-only.

    Plain ``sqlite3.connect`` silently creates an empty database when the file
    is missing, which later surfaces as a confusing "no such table" error, so a
    missing file is reported up front instead.
    """
    path = Path(db_path).resolve()
    if not path.is_file():
        raise FileNotFoundError(
            f"Database not found at {path}. Build it first with `make pipeline` (or `python load_data.py`)."
        )
    return sqlite3.connect(f"{path.as_uri()}?mode=ro", uri=True)


@contextmanager
def open_db(db_path: str | Path) -> Iterator[sqlite3.Connection]:
    """``with open_db(path) as conn:`` - a read-only connection that is always closed."""
    conn = connect(db_path)
    try:
        yield conn
    finally:
        conn.close()
