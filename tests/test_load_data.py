"""Part 1 - schema and loader."""

from __future__ import annotations

import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import load_data

REPO_ROOT = Path(__file__).resolve().parents[1]
HEADER = ",".join(load_data.REQUIRED_COLUMNS)
ROW = "prj1,sbj1,melanoma,60,F,miraclib,yes,s1,PBMC,0,10,20,30,40,50"


def write_csv(path: Path, rows: list[str]) -> Path:
    path.write_text("\n".join([HEADER, *rows]) + "\n")
    return path


def test_loads_every_row(conn, csv_frame):
    assert load_data.table_counts(conn) == {
        "projects": 3,
        "subjects": 3_500,
        "samples": 10_500,
        "cell_populations": 5,
        "cell_counts": 52_500,
    }
    assert len(csv_frame) == 10_500


def test_database_round_trips_the_csv(conn, csv_frame):
    """Re-assembling the wide CSV from the normalised tables gives back every original value."""
    metadata = pd.read_sql_query("SELECT * FROM sample_metadata", conn)
    counts = pd.read_sql_query("SELECT sample, population, count FROM cell_frequencies", conn).pivot(
        index="sample", columns="population", values="count"
    )
    rebuilt = metadata.merge(counts.reset_index(), on="sample")
    rebuilt["response"] = rebuilt["response"].fillna("")
    rebuilt = rebuilt[list(csv_frame.columns)].sort_values("sample").reset_index(drop=True)
    original = csv_frame.sort_values("sample").reset_index(drop=True)
    pd.testing.assert_frame_equal(rebuilt, original, check_dtype=False)


def test_blank_response_is_stored_as_null(conn):
    untreated = conn.execute("SELECT DISTINCT condition, treatment FROM subjects WHERE response IS NULL").fetchall()
    assert untreated == [("healthy", "none")]
    assert conn.execute("SELECT COUNT(*) FROM subjects WHERE response = ''").fetchone() == (0,)


def test_rebuilding_starts_from_scratch(tmp_path):
    db = tmp_path / "x.db"
    csv_path = write_csv(tmp_path / "a.csv", [ROW])
    load_data.build_database(csv_path, db)
    load_data.build_database(csv_path, db)
    with sqlite3.connect(db) as c:
        assert load_data.table_counts(c)["samples"] == 1
    assert not (tmp_path / "x.db.tmp").exists()


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        ([ROW, ROW], "duplicate sample id"),
        ([ROW, ROW.replace(",60,", ",61,").replace(",s1,", ",s2,")], "conflicting attributes"),
        ([ROW.replace(",0,10,", ",0,ten,")], "must be an integer"),
        ([ROW.replace("melanoma", "")], "'condition' is empty"),
    ],
)
def test_invalid_rows_are_rejected(tmp_path, rows, message):
    with pytest.raises(load_data.DataValidationError, match=message):
        load_data.build_database(write_csv(tmp_path / "bad.csv", rows), tmp_path / "x.db")


def test_missing_column_is_rejected(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("project,subject\nprj1,sbj1\n")
    with pytest.raises(load_data.DataValidationError, match="missing column"):
        load_data.build_database(path, tmp_path / "x.db")


def test_failed_build_keeps_the_previous_database(tmp_path):
    db = tmp_path / "x.db"
    load_data.build_database(write_csv(tmp_path / "good.csv", [ROW]), db)
    with pytest.raises(load_data.DataValidationError):
        load_data.build_database(write_csv(tmp_path / "bad.csv", [ROW, ROW]), db)
    with sqlite3.connect(db) as c:
        assert load_data.table_counts(c)["samples"] == 1
    assert not (tmp_path / "x.db.tmp").exists()


def test_schema_enforces_constraints(tmp_path):
    db = tmp_path / "x.db"
    load_data.build_database(write_csv(tmp_path / "good.csv", [ROW]), db)
    conn = sqlite3.connect(db)
    conn.execute("PRAGMA foreign_keys = ON")
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute("UPDATE cell_counts SET cell_count = -1")
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute("UPDATE subjects SET response = 'maybe'")
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute("INSERT INTO samples (sample_code, subject_id, sample_type) VALUES ('s9', 999, 'PBMC')")
    with pytest.raises(sqlite3.IntegrityError):
        conn.execute("INSERT INTO samples (sample_code, subject_id, sample_type) VALUES ('s1', 1, 'PBMC')")
    conn.close()


def test_script_runs_without_arguments(tmp_path):
    """`python load_data.py` creates cell_counts.db next to it; -S (no site-packages) proves it is stdlib-only."""
    shutil.copy(REPO_ROOT / "load_data.py", tmp_path)
    shutil.copy(REPO_ROOT / "cell-count.csv", tmp_path)
    result = subprocess.run(
        [sys.executable, "-I", "-S", "load_data.py"], cwd=tmp_path, capture_output=True, text=True, check=True
    )
    assert (tmp_path / "cell_counts.db").is_file()
    assert "10,500" in result.stdout
