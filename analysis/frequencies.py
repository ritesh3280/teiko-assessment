"""Part 2 - relative frequency of each immune population in each sample.

For every sample the total cell count is the sum over the five populations, and
each population's relative frequency is its share of that total, in percent.
The arithmetic lives in the ``cell_frequencies`` SQL view (see ``load_data.py``),
so the database itself answers the question and every consumer agrees on it.
"""

from __future__ import annotations

import sqlite3

import pandas as pd

FREQUENCY_COLUMNS = ["sample", "total_count", "population", "count", "percentage"]


def population_names(conn: sqlite3.Connection) -> list[str]:
    """Population names in their canonical (CSV column) order."""
    rows = conn.execute("SELECT population_name FROM cell_populations ORDER BY population_id")
    return [name for (name,) in rows]


def cell_frequency_table(conn: sqlite3.Connection) -> pd.DataFrame:
    """The Part 2 summary table: one row per (sample, population).

    Columns: ``sample``, ``total_count``, ``population``, ``count``, ``percentage``.
    Rows follow the CSV's sample order, then the canonical population order.
    """
    return pd.read_sql_query(
        """
        SELECT sample, total_count, population, count, percentage
        FROM cell_frequencies
        ORDER BY sample_id, population_id
        """,
        conn,
    )


def population_summary(frequencies: pd.DataFrame) -> pd.DataFrame:
    """Distribution of each population's relative frequency (%) across samples."""
    grouped = frequencies.groupby("population", sort=False)["percentage"]
    summary = grouped.describe(percentiles=[0.25, 0.5, 0.75])
    summary = summary.rename(columns={"25%": "q1", "50%": "median", "75%": "q3", "count": "n_samples"})
    summary["n_samples"] = summary["n_samples"].astype(int)
    return summary[["n_samples", "mean", "std", "min", "q1", "median", "q3", "max"]].reset_index()
