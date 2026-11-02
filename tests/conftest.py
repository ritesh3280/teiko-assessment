"""Shared fixtures: a database built from the real CSV in a temporary directory.

Tests never touch the repository's own ``cell_counts.db``.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd
import pytest

import load_data
from analysis.db import connect


@pytest.fixture(scope="session")
def db_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    path = tmp_path_factory.mktemp("db") / "cell_counts.db"
    load_data.build_database(load_data.CSV_PATH, path)
    return path


@pytest.fixture()
def conn(db_path: Path):
    connection: sqlite3.Connection = connect(db_path)  # read-only, like the pipeline and dashboard
    yield connection
    connection.close()


@pytest.fixture(scope="session")
def csv_frame() -> pd.DataFrame:
    """The raw CSV read with pandas: an independent reference for the SQL results."""
    return pd.read_csv(load_data.CSV_PATH, keep_default_na=False)
