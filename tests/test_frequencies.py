"""Part 2 - the relative-frequency summary table."""

from __future__ import annotations

import numpy as np
import pytest

from analysis.db import open_db
from analysis.frequencies import FREQUENCY_COLUMNS, cell_frequency_table, population_names, population_summary
from load_data import POPULATIONS


@pytest.fixture(scope="module")
def table(db_path):
    with open_db(db_path) as conn:
        return cell_frequency_table(conn)


def test_columns_and_shape(table):
    assert list(table.columns) == FREQUENCY_COLUMNS == ["sample", "total_count", "population", "count", "percentage"]
    assert len(table) == 10_500 * 5
    assert table.groupby("sample").size().eq(5).all()


def test_matches_an_independent_calculation(table, csv_frame):
    long = csv_frame.melt(id_vars="sample", value_vars=list(POPULATIONS), var_name="population", value_name="count")
    long["total_count"] = long.groupby("sample")["count"].transform("sum")
    long["percentage"] = 100 * long["count"] / long["total_count"]
    merged = table.merge(long, on=["sample", "population"], suffixes=("", "_expected"), validate="one_to_one")
    assert len(merged) == len(table)
    assert (merged["count"] == merged["count_expected"]).all()
    assert (merged["total_count"] == merged["total_count_expected"]).all()
    np.testing.assert_allclose(merged["percentage"], merged["percentage_expected"], rtol=1e-12)


def test_percentages_sum_to_100_per_sample(table):
    np.testing.assert_allclose(table.groupby("sample")["percentage"].sum(), 100.0, atol=1e-9)


def test_known_sample(table):
    rows = table[table["sample"] == "sample00000"].set_index("population")
    assert rows["total_count"].unique().tolist() == [93_214]
    assert rows.loc["b_cell", "count"] == 10_908
    assert rows.loc["b_cell", "percentage"] == pytest.approx(100 * 10_908 / 93_214)


def test_population_order(db_path, table):
    with open_db(db_path) as conn:
        assert population_names(conn) == list(POPULATIONS)
    assert table["population"].iloc[:5].tolist() == list(POPULATIONS)


def test_population_summary(table):
    summary = population_summary(table)
    assert summary["population"].tolist() == list(POPULATIONS)
    assert (summary["n_samples"] == 10_500).all()
    assert summary["mean"].sum() == pytest.approx(100.0)  # mean compositions still add up to 100%
