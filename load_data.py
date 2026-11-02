#!/usr/bin/env python3
"""Part 1 - Data management: build the SQLite database from ``cell-count.csv``.

Usage (from the repository root, no arguments)::

    python load_data.py

This creates ``cell_counts.db`` in the repository root, initialises the schema
defined in ``SCHEMA`` below and loads every row of ``cell-count.csv``. Re-running
the script rebuilds the database from scratch, so it is safe to run repeatedly.

Schema (third normal form)::

    projects 1──< subjects 1──< samples 1──< cell_counts >──1 cell_populations

* ``projects``          one row per project / study.
* ``subjects``          one row per patient. Attributes that are fixed for a patient
                        (indication, age, sex, treatment, response) are stored once
                        here instead of being repeated on every sample.
* ``samples``           one row per biological sample (sample type, timepoint).
* ``cell_populations``  lookup table of the measured immune populations.
* ``cell_counts``       one row per (sample, population) - a "long" table, so adding
                        a population or a new panel adds rows, not columns.

Two views cover the common read patterns:

* ``sample_metadata``   one row per sample with its subject and project attributes.
* ``cell_frequencies``  the Part 2 summary table: per-sample total cell count and
                        the relative frequency (%) of each population.

The script uses only the Python standard library, so it runs on a bare Python
install, before ``make setup`` has installed anything.
"""

from __future__ import annotations

import csv
import os
import sqlite3
import sys
from collections.abc import Iterable, Mapping
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
CSV_PATH = REPO_ROOT / "cell-count.csv"
DB_PATH = REPO_ROOT / "cell_counts.db"

# Immune populations measured in every sample, in CSV column order. The order is
# kept as ``cell_populations.population_id`` and used for display everywhere.
POPULATIONS: tuple[str, ...] = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")

REQUIRED_COLUMNS: tuple[str, ...] = (
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
    *POPULATIONS,
)

SCHEMA = """
CREATE TABLE projects (
    project_id    INTEGER PRIMARY KEY,
    project_code  TEXT    NOT NULL UNIQUE                -- e.g. 'prj1'
);

CREATE TABLE subjects (
    subject_id    INTEGER PRIMARY KEY,
    project_id    INTEGER NOT NULL REFERENCES projects (project_id),
    subject_code  TEXT    NOT NULL,                      -- e.g. 'sbj000'
    condition     TEXT    NOT NULL,                      -- indication: melanoma, carcinoma, healthy
    age           INTEGER CHECK (age >= 0),
    sex           TEXT    CHECK (sex IN ('M', 'F')),
    treatment     TEXT    NOT NULL,                      -- miraclib, phauximab, none
    response      TEXT    CHECK (response IN ('yes', 'no')),  -- NULL when not applicable
    UNIQUE (project_id, subject_code)                    -- subject codes may repeat across projects
);

CREATE TABLE samples (
    sample_id                  INTEGER PRIMARY KEY,
    sample_code                TEXT    NOT NULL UNIQUE,  -- e.g. 'sample00000'
    subject_id                 INTEGER NOT NULL REFERENCES subjects (subject_id),
    sample_type                TEXT    NOT NULL,         -- PBMC, WB
    time_from_treatment_start  INTEGER                   -- days; 0 = baseline
);

CREATE TABLE cell_populations (
    population_id    INTEGER PRIMARY KEY,
    population_name  TEXT    NOT NULL UNIQUE             -- e.g. 'cd4_t_cell'
);

CREATE TABLE cell_counts (
    sample_id      INTEGER NOT NULL REFERENCES samples (sample_id) ON DELETE CASCADE,
    population_id  INTEGER NOT NULL REFERENCES cell_populations (population_id),
    cell_count     INTEGER NOT NULL CHECK (cell_count >= 0),
    PRIMARY KEY (sample_id, population_id)
) WITHOUT ROWID;

-- Indexes on foreign keys and on the columns the analyses filter by.
CREATE INDEX idx_subjects_project        ON subjects (project_id);
CREATE INDEX idx_subjects_cohort         ON subjects (condition, treatment, response);
CREATE INDEX idx_samples_subject         ON samples (subject_id);
CREATE INDEX idx_samples_type_timepoint  ON samples (sample_type, time_from_treatment_start);
CREATE INDEX idx_cell_counts_population  ON cell_counts (population_id);

-- One row per sample with its subject and project attributes.
CREATE VIEW sample_metadata AS
SELECT s.sample_id,
       s.subject_id,
       s.sample_code    AS sample,
       p.project_code   AS project,
       sub.subject_code AS subject,
       sub.condition,
       sub.age,
       sub.sex,
       sub.treatment,
       sub.response,
       s.sample_type,
       s.time_from_treatment_start
FROM samples s
JOIN subjects sub ON sub.subject_id = s.subject_id
JOIN projects p   ON p.project_id = sub.project_id;

-- Part 2: relative frequency (%) of each population within each sample.
CREATE VIEW cell_frequencies AS
SELECT cc.sample_id,
       s.sample_code       AS sample,
       SUM(cc.cell_count) OVER (PARTITION BY cc.sample_id) AS total_count,
       cp.population_id,
       cp.population_name  AS population,
       cc.cell_count       AS count,
       100.0 * cc.cell_count / SUM(cc.cell_count) OVER (PARTITION BY cc.sample_id) AS percentage
FROM cell_counts cc
JOIN samples s           ON s.sample_id = cc.sample_id
JOIN cell_populations cp ON cp.population_id = cc.population_id;
"""


class DataValidationError(ValueError):
    """The input CSV violates an assumption the schema relies on."""


def _text(row: Mapping[str, str], column: str, line: int) -> str:
    value = (row.get(column) or "").strip()
    if not value:
        raise DataValidationError(f"line {line}: '{column}' is empty")
    return value


def _optional_text(row: Mapping[str, str], column: str) -> str | None:
    return (row.get(column) or "").strip() or None


def _integer(row: Mapping[str, str], column: str, line: int, *, required: bool = True) -> int | None:
    value = (row.get(column) or "").strip()
    if not value:
        if required:
            raise DataValidationError(f"line {line}: '{column}' is empty")
        return None
    try:
        return int(value)
    except ValueError:
        raise DataValidationError(f"line {line}: '{column}' must be an integer, got {value!r}") from None


def read_csv(csv_path: Path) -> list[dict[str, str]]:
    """Read the CSV and check that every required column is present."""
    with open(csv_path, newline="", encoding="utf-8-sig") as handle:  # -sig tolerates a BOM
        reader = csv.DictReader(handle)
        reader.fieldnames = [name.strip() for name in reader.fieldnames or []]
        missing = [column for column in REQUIRED_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise DataValidationError(f"{csv_path.name} is missing column(s): {', '.join(missing)}")
        return list(reader)


def create_schema(conn: sqlite3.Connection) -> None:
    """Create all tables, indexes and views in an empty database."""
    conn.executescript(SCHEMA)


def load_rows(conn: sqlite3.Connection, rows: Iterable[Mapping[str, str]]) -> None:
    """Normalise CSV rows into the schema, inserting everything in one transaction.

    Raises ``DataValidationError`` if a value cannot be parsed, a sample id is
    duplicated, or one subject appears with conflicting attributes (which would
    mean the subject-level columns are not really subject-level).
    """
    projects: dict[str, None] = {}  # insertion-ordered set of project codes
    subjects: dict[tuple[str, str], tuple] = {}  # (project, subject) -> attributes
    samples: dict[str, tuple] = {}  # sample code -> (subject key, sample_type, timepoint)
    counts: dict[str, list[int]] = {}  # sample code -> counts in POPULATIONS order

    for line, row in enumerate(rows, start=2):  # line 1 is the header
        project = _text(row, "project", line)
        subject_key = (project, _text(row, "subject", line))
        attributes = (
            _text(row, "condition", line),
            _integer(row, "age", line, required=False),
            _optional_text(row, "sex"),
            _text(row, "treatment", line),
            _optional_text(row, "response"),  # blank for untreated / healthy subjects
        )
        known = subjects.setdefault(subject_key, attributes)
        if known != attributes:
            raise DataValidationError(
                f"line {line}: subject {subject_key[1]!r} in {project!r} has conflicting "
                f"attributes {attributes} (first seen as {known})"
            )

        sample = _text(row, "sample", line)
        if sample in samples:
            raise DataValidationError(f"line {line}: duplicate sample id {sample!r}")
        projects.setdefault(project, None)
        samples[sample] = (
            subject_key,
            _text(row, "sample_type", line),
            _integer(row, "time_from_treatment_start", line, required=False),
        )
        counts[sample] = [_integer(row, population, line) for population in POPULATIONS]

    with conn:  # commit on success, roll back everything on error
        conn.executemany("INSERT INTO projects (project_code) VALUES (?)", ((code,) for code in projects))
        project_ids = dict(conn.execute("SELECT project_code, project_id FROM projects"))

        conn.executemany(
            "INSERT INTO subjects (project_id, subject_code, condition, age, sex, treatment, response)"
            " VALUES (?, ?, ?, ?, ?, ?, ?)",
            ((project_ids[project], code, *attrs) for (project, code), attrs in subjects.items()),
        )
        subject_ids = {
            (project, code): subject_id
            for subject_id, project, code in conn.execute(
                "SELECT sub.subject_id, p.project_code, sub.subject_code"
                " FROM subjects sub JOIN projects p ON p.project_id = sub.project_id"
            )
        }

        conn.executemany(
            "INSERT INTO samples (sample_code, subject_id, sample_type, time_from_treatment_start)"
            " VALUES (?, ?, ?, ?)",
            ((code, subject_ids[key], sample_type, time) for code, (key, sample_type, time) in samples.items()),
        )
        sample_ids = dict(conn.execute("SELECT sample_code, sample_id FROM samples"))

        conn.executemany(
            "INSERT INTO cell_populations (population_id, population_name) VALUES (?, ?)",
            enumerate(POPULATIONS, start=1),
        )
        conn.executemany(
            "INSERT INTO cell_counts (sample_id, population_id, cell_count) VALUES (?, ?, ?)",
            (
                (sample_ids[code], population_id, value)
                for code, values in counts.items()
                for population_id, value in enumerate(values, start=1)
            ),
        )


def table_counts(conn: sqlite3.Connection) -> dict[str, int]:
    """Row count of every table, in schema order."""
    tables = ("projects", "subjects", "samples", "cell_populations", "cell_counts")
    return {table: conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0] for table in tables}


def build_database(csv_path: Path = CSV_PATH, db_path: Path = DB_PATH) -> dict[str, int]:
    """(Re)build ``db_path`` from ``csv_path`` and return the row count of each table.

    The database is written to a temporary file and moved into place only after
    the load and integrity checks succeed, so a failed run never leaves a
    half-built database behind (and a running dashboard never reads one).
    """
    csv_path, db_path = Path(csv_path), Path(db_path)
    rows = read_csv(csv_path)
    tmp_path = db_path.with_name(db_path.name + ".tmp")
    tmp_path.unlink(missing_ok=True)
    try:
        conn = sqlite3.connect(tmp_path)
        try:
            conn.execute("PRAGMA foreign_keys = ON")
            create_schema(conn)
            load_rows(conn, rows)
            if conn.execute("PRAGMA foreign_key_check").fetchall():
                raise DataValidationError("foreign key check failed")
            counts = table_counts(conn)
        finally:
            conn.close()
        expected = {"samples": len(rows), "cell_counts": len(rows) * len(POPULATIONS)}
        if any(counts[table] != n for table, n in expected.items()):
            raise DataValidationError(f"row count mismatch: loaded {counts}, expected {expected}")
        os.replace(tmp_path, db_path)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise
    return counts


def main() -> int:
    try:
        counts = build_database()
    except (OSError, sqlite3.Error, DataValidationError) as exc:
        print(f"load_data.py: error: {exc}", file=sys.stderr)
        return 1
    print(f"Built {DB_PATH.name} from {CSV_PATH.name}:")
    for table, n in counts.items():
        print(f"  {table:<17} {n:>7,} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
