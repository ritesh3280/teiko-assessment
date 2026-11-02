"""Part 4 - baseline subset queries, checked against an independent pandas computation."""

from __future__ import annotations

import pytest

from analysis import subsets


@pytest.fixture(scope="module")
def baseline(csv_frame):
    return csv_frame.query(
        "condition == 'melanoma' and treatment == 'miraclib' and sample_type == 'PBMC' and time_from_treatment_start == 0"
    )


def as_dict(df, key: str, value: str) -> dict:
    return dict(zip(df[key], df[value]))


def test_subset_samples(conn, baseline):
    samples = subsets.subset_samples(conn)
    assert sorted(samples["sample"]) == sorted(baseline["sample"])
    assert len(samples) == 656
    assert list(samples.columns) == ["sample", "project", "subject", "response", "sex", "age"]


def test_samples_per_project(conn, baseline):
    per_project = subsets.samples_per_project(conn)
    assert as_dict(per_project, "project", "n_samples") == {"prj1": 384, "prj2": 0, "prj3": 272}
    assert baseline["project"].value_counts().to_dict() == {"prj1": 384, "prj3": 272}


def test_subjects_by_response(conn, baseline):
    by_response = subsets.subjects_by_response(conn)
    assert as_dict(by_response, "response", "n_subjects") == {"yes": 331, "no": 325}
    assert baseline.groupby("response")["subject"].nunique().to_dict() == {"no": 325, "yes": 331}


def test_subjects_by_sex(conn, baseline):
    by_sex = subsets.subjects_by_sex(conn)
    assert as_dict(by_sex, "sex", "n_subjects") == {"M": 344, "F": 312}
    assert baseline.groupby("sex")["subject"].nunique().to_dict() == {"F": 312, "M": 344}


def test_mean_b_cells_for_male_melanoma_responders(conn, csv_frame):
    """Melanoma males of all sample and treatment types: average B cells of responders at time 0."""
    answer = subsets.mean_cell_count(conn).iloc[0]
    expected = csv_frame.query(
        "condition == 'melanoma' and sex == 'M' and response == 'yes' and time_from_treatment_start == 0"
    )
    assert set(expected["treatment"]) == {"miraclib", "phauximab"}  # every treatment ...
    assert set(expected["sample_type"]) == {"PBMC", "WB"}  # ... and every sample type is included
    assert answer["n_samples"] == len(expected) == 485
    assert answer["n_subjects"] == expected["subject"].nunique() == 485
    assert answer["mean_count"] == pytest.approx(expected["b_cell"].mean())
    assert f"{answer['mean_count']:.2f}" == "10206.15"


def test_other_filters(conn, csv_frame):
    healthy = subsets.SampleFilter(condition="healthy", treatment="none", sample_type="PBMC", timepoint=7)
    by_response = subsets.subjects_by_response(conn, healthy)
    assert by_response["response"].isna().all()  # untreated subjects have no response
    expected = csv_frame.query("condition == 'healthy' and sample_type == 'PBMC' and time_from_treatment_start == 7")
    assert by_response["n_subjects"].sum() == expected["subject"].nunique()
    assert subsets.subset_samples(conn, subsets.SampleFilter(condition="unknown")).empty


def test_filter_options(conn):
    options = subsets.filter_options(conn)
    assert options["condition"] == ["carcinoma", "healthy", "melanoma"]
    assert options["treatment"] == ["miraclib", "none", "phauximab"]
    assert options["sample_type"] == ["PBMC", "WB"]
    assert options["timepoint"] == [0, 7, 14]
