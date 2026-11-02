"""Part 3 - responders vs non-responders."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from analysis import responders
from analysis.db import open_db


@pytest.fixture(scope="module")
def cohort(db_path):
    with open_db(db_path) as conn:
        return responders.cohort_frequencies(conn)


@pytest.fixture(scope="module")
def results(cohort):
    return responders.compare_responders(cohort)


@pytest.fixture(scope="module")
def by_timepoint(cohort):
    return responders.compare_by_timepoint(cohort)


def synthetic_cohort(shift: float, n_subjects: int = 120, seed: int = 1) -> pd.DataFrame:
    """Patients sampled on days 0/7/14; 'shifted' is higher in responders by ``shift`` points."""
    rng = np.random.default_rng(seed)
    rows = []
    for subject in range(n_subjects):
        response = "yes" if subject % 2 else "no"
        patient_effect = rng.normal(0, 1)
        for timepoint in (0, 7, 14):
            for population, effect in (("shifted", shift), ("unchanged", 0.0)):
                value = 20 + patient_effect + (effect if response == "yes" else 0.0) + rng.normal(0, 2)
                rows.append(
                    {
                        "sample": f"s{subject}_{timepoint}",
                        "subject_id": subject,
                        "subject": f"p{subject}",
                        "response": response,
                        "timepoint": timepoint,
                        "population": population,
                        "percentage": value,
                    }
                )
    return pd.DataFrame(rows)


def test_cohort_matches_the_csv(cohort, csv_frame):
    expected = csv_frame.query("condition == 'melanoma' and treatment == 'miraclib' and sample_type == 'PBMC'")
    assert set(cohort["sample"]) == set(expected["sample"])
    counts = responders.cohort_counts(cohort).set_index("response")
    assert counts.loc["yes", "n_subjects"] == expected.query("response == 'yes'")["subject"].nunique() == 331
    assert counts.loc["no", "n_subjects"] == expected.query("response == 'no'")["subject"].nunique() == 325
    assert counts["n_samples"].sum() == len(expected) == 1_968


def test_cd4_t_cells_are_the_only_significant_population(results):
    """Pins the headline finding reported in the README and report."""
    assert results["population"].tolist() == ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
    assert results.loc[results["significant"], "population"].tolist() == ["cd4_t_cell"]
    cd4 = results.set_index("population").loc["cd4_t_cell"]
    assert cd4["lmm_difference"] == pytest.approx(0.636, abs=0.001)
    assert cd4["lmm_q"] == pytest.approx(0.0248, abs=0.0005)
    assert cd4["mwu_p"] == pytest.approx(0.0133, abs=0.0005)
    assert cd4["auc"] == pytest.approx(0.532, abs=0.001)


def test_mixed_model_agrees_with_a_plain_difference_in_means(results, cohort):
    """With every patient sampled at all three timepoints the design is balanced, so the
    timepoint-adjusted model estimate equals the raw difference in group means."""
    means = cohort.groupby(["population", "response"])["percentage"].mean().unstack()
    raw = (means["yes"] - means["no"]).reindex(results["population"]).to_numpy()
    np.testing.assert_allclose(results["lmm_difference"], raw, atol=1e-6)


def test_fdr_adjustment(results):
    assert (results["lmm_q"] >= results["lmm_p"] - 1e-12).all()
    order = results["lmm_p"].argsort().to_numpy()
    assert np.all(np.diff(results["lmm_q"].to_numpy()[order]) >= -1e-12)  # q is monotone in p


def test_model_fitting_does_not_leak_warnings(cohort):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        responders.compare_responders(cohort)
    assert [str(w.message) for w in caught if "statsmodels" in w.filename] == []


def test_detects_a_planted_difference():
    results = responders.compare_responders(synthetic_cohort(shift=2.0)).set_index("population")
    assert results.loc["shifted", "significant"]
    assert results.loc["shifted", "lmm_difference"] == pytest.approx(2.0, abs=0.6)
    assert results.loc["shifted", "auc"] > 0.6
    assert results.loc["shifted", "patient_variance"] > 0  # the simulated patient effect is recovered
    assert not results.loc["unchanged", "significant"]


def test_needs_both_groups():
    with pytest.raises(ValueError, match="no non-responders"):
        responders.compare_responders(synthetic_cohort(shift=1.0).query("response == 'yes'"))


def test_timepoint_comparison(by_timepoint):
    assert len(by_timepoint) == 15
    assert by_timepoint["timepoint"].unique().tolist() == [0, 7, 14]
    assert (by_timepoint["n_responders"] == 331).all()
    assert (by_timepoint["n_non_responders"] == 325).all()
    assert (by_timepoint.query("timepoint == 0")["mwu_p"] > 0.05).all()  # nothing differs at baseline


def test_timepoint_means(cohort):
    means = responders.timepoint_means(cohort)
    assert len(means) == 5 * 3 * 2
    assert ((means["ci_low"] < means["mean"]) & (means["mean"] < means["ci_high"])).all()


def test_findings_are_generated_from_the_results(results, by_timepoint):
    findings = responders.key_findings(results, by_timepoint)
    assert findings[0].startswith(
        "**cd4_t_cell** shows a statistically significant association with response under the primary"
    )
    # The secondary adjusted rank test is not significant, and the text has to say so.
    assert any("but is not significant after FDR correction (q = 0.067)" in f for f in findings)
    assert any("limited predictive value" in f for f in findings)
    assert any("no cd4_t_cell difference at baseline" in f and "none of these survives" in f for f in findings)
    assert any(f.startswith("At baseline (day 0) no population differs") for f in findings)


def test_findings_cover_a_baseline_difference():
    """Hand-made inputs exercising the branches the real data does not reach."""
    results = pd.DataFrame(
        {
            "population": ["x", "y"],
            "lmm_difference": [1.0, 0.1],
            "lmm_ci_low": [0.5, -0.2],
            "lmm_ci_high": [1.5, 0.4],
            "lmm_p": [0.001, 0.5],
            "lmm_q": [0.002, 0.5],
            "mwu_p": [0.004, 0.4],
            "mwu_q": [0.008, 0.4],
            "auc": [0.7, 0.51],
            "significant": [True, False],
        }
    )
    by_timepoint = pd.DataFrame(
        {
            "population": ["x", "x", "y", "y"],
            "timepoint": [0, 7, 0, 7],
            "mwu_p": [0.01, 0.02, 0.6, 0.001],
            "mwu_q": [0.02, 0.03, 0.6, 0.004],
        }
    )
    findings = responders.key_findings(results, by_timepoint)
    assert any("already present at baseline" in f for f in findings)
    assert any("survive FDR correction" in f and "y on day 7" in f for f in findings)
    assert not any(f.startswith("At baseline (day 0) no population differs") for f in findings)


def test_format_p():
    assert responders.format_p(0.0004) == "< 0.001"
    assert responders.format_p(0.02479) == "0.025"
