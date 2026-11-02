"""Part 4 - explore a subset of samples.

The default filter answers Bob's question: melanoma PBMC samples at baseline
(``time_from_treatment_start = 0``) from patients treated with miraclib. Every
subset query below starts from that same filtered set (the ``subset`` CTE), and
the aggregation is done in SQL so the database answers each question directly.

``mean_cell_count`` answers a separate question that spans *all* treatments and
sample types: the average number of B cells in male melanoma responders at
baseline.
"""

from __future__ import annotations

import sqlite3
from dataclasses import asdict, dataclass

import pandas as pd


@dataclass(frozen=True)
class SampleFilter:
    """Which samples to select; the defaults are Bob's baseline question."""

    condition: str = "melanoma"
    treatment: str = "miraclib"
    sample_type: str = "PBMC"
    timepoint: int = 0  # time_from_treatment_start, in days

    def describe(self) -> str:
        return f"{self.condition} · {self.treatment} · {self.sample_type} · day {self.timepoint}"


BASELINE = SampleFilter()

SUBSET_CTE = """\
WITH subset AS (
    SELECT s.sample_id, s.sample_code, sub.subject_id, sub.subject_code,
           sub.project_id, sub.response, sub.sex, sub.age
    FROM samples AS s
    JOIN subjects AS sub ON sub.subject_id = s.subject_id
    WHERE sub.condition = :condition
      AND sub.treatment = :treatment
      AND s.sample_type = :sample_type
      AND s.time_from_treatment_start = :timepoint
)
"""

QUERIES = {
    # Step 1: every matching sample, with its project and subject attributes.
    "samples": """\
SELECT subset.sample_code AS sample, p.project_code AS project,
       subset.subject_code AS subject, subset.response, subset.sex, subset.age
FROM subset
JOIN projects AS p ON p.project_id = subset.project_id
ORDER BY subset.sample_id
""",
    # Step 2a: matching samples per project; projects with none are listed as 0.
    "samples_per_project": """\
SELECT p.project_code AS project, COUNT(subset.sample_id) AS n_samples
FROM projects AS p
LEFT JOIN subset ON subset.project_id = p.project_id
GROUP BY p.project_id
ORDER BY p.project_id
""",
    # Step 2b: distinct subjects who are responders ('yes') / non-responders ('no').
    "subjects_by_response": """\
SELECT response, COUNT(DISTINCT subject_id) AS n_subjects
FROM subset
GROUP BY response
ORDER BY response DESC
""",
    # Step 2c: distinct subjects who are male ('M') / female ('F').
    "subjects_by_sex": """\
SELECT sex, COUNT(DISTINCT subject_id) AS n_subjects
FROM subset
GROUP BY sex
ORDER BY sex DESC
""",
}


def sql(name: str) -> str:
    """The full SQL text (filter CTE + query) behind one of the functions below."""
    return SUBSET_CTE + QUERIES[name]


def _query(conn: sqlite3.Connection, name: str, sample_filter: SampleFilter) -> pd.DataFrame:
    return pd.read_sql_query(sql(name), conn, params=asdict(sample_filter))


def subset_samples(conn: sqlite3.Connection, sample_filter: SampleFilter = BASELINE) -> pd.DataFrame:
    """Step 1: every matching sample, with its project and subject attributes."""
    return _query(conn, "samples", sample_filter)


def samples_per_project(conn: sqlite3.Connection, sample_filter: SampleFilter = BASELINE) -> pd.DataFrame:
    """Step 2a: number of matching samples in each project; projects with none are listed as 0."""
    return _query(conn, "samples_per_project", sample_filter)


def subjects_by_response(conn: sqlite3.Connection, sample_filter: SampleFilter = BASELINE) -> pd.DataFrame:
    """Step 2b: distinct subjects who are responders ('yes') and non-responders ('no')."""
    return _query(conn, "subjects_by_response", sample_filter)


def subjects_by_sex(conn: sqlite3.Connection, sample_filter: SampleFilter = BASELINE) -> pd.DataFrame:
    """Step 2c: distinct subjects who are male ('M') and female ('F')."""
    return _query(conn, "subjects_by_sex", sample_filter)


@dataclass(frozen=True)
class CountQuestion:
    """Average raw count of one population, over samples of every treatment and sample type."""

    population: str = "b_cell"
    condition: str = "melanoma"
    sex: str = "M"
    response: str = "yes"
    timepoint: int = 0  # time_from_treatment_start, in days

    def describe(self) -> str:
        sex = {"M": "male", "F": "female"}.get(self.sex, self.sex)
        group = {"yes": "responders", "no": "non-responders"}.get(self.response, self.response)
        return f"{self.condition} · {sex} · {group} · day {self.timepoint} · all treatments and sample types"


# "Considering melanoma males of all sample and treatment types, what is the
# average number of B cells for responders at time = 0?"
B_CELL_QUESTION = CountQuestion()

# Deliberately no treatment or sample-type filter, unlike the subset queries above.
MEAN_COUNT_SQL = """\
SELECT COUNT(*)                       AS n_samples,
       COUNT(DISTINCT sub.subject_id) AS n_subjects,
       AVG(cc.cell_count)             AS mean_count
FROM cell_counts AS cc
JOIN cell_populations AS cp ON cp.population_id = cc.population_id
JOIN samples AS s           ON s.sample_id = cc.sample_id
JOIN subjects AS sub        ON sub.subject_id = s.subject_id
WHERE cp.population_name = :population
  AND sub.condition = :condition
  AND sub.sex = :sex
  AND sub.response = :response
  AND s.time_from_treatment_start = :timepoint
"""


def mean_cell_count(conn: sqlite3.Connection, question: CountQuestion = B_CELL_QUESTION) -> pd.DataFrame:
    """One row: the question's filters, then ``n_samples``, ``n_subjects`` and ``mean_count``."""
    params = asdict(question)
    answer = pd.read_sql_query(MEAN_COUNT_SQL, conn, params=params)
    return pd.concat([pd.DataFrame([params]), answer], axis=1)


def filter_options(conn: sqlite3.Connection) -> dict[str, list]:
    """Distinct values of each filter field, for building selectors."""

    def distinct(sql: str) -> list:
        return [value for (value,) in conn.execute(sql)]

    return {
        "condition": distinct("SELECT DISTINCT condition FROM subjects ORDER BY condition"),
        "treatment": distinct("SELECT DISTINCT treatment FROM subjects ORDER BY treatment"),
        "sample_type": distinct("SELECT DISTINCT sample_type FROM samples ORDER BY sample_type"),
        "timepoint": distinct(
            "SELECT DISTINCT time_from_treatment_start FROM samples"
            " WHERE time_from_treatment_start IS NOT NULL ORDER BY 1"
        ),
    }
