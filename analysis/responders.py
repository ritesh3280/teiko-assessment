"""Part 3 - do population frequencies differ between responders and non-responders?

Cohort (default): melanoma patients treated with miraclib, PBMC samples only.
Responders have ``response = 'yes'``, non-responders ``response = 'no'``.

Statistical approach
--------------------
Every patient contributes a sample at day 0, 7 and 14, so a patient's samples
are not independent observations. The **primary analysis** is therefore a linear
mixed-effects model fitted separately for each population::

    percentage ~ responder + C(timepoint) + (1 | patient)

Its ``responder`` coefficient is the mean difference in relative frequency
(percentage points, responders minus non-responders) adjusted for timepoint,
with a random intercept per patient to account for repeated sampling. P-values
are adjusted across the five populations with the Benjamini-Hochberg false
discovery rate (FDR), and a population is called significant when q < 0.05.

Two checks are reported alongside it:

* a two-sided **Mann-Whitney U** test on the sample-level frequencies (the data in
  the boxplots). It makes no distributional assumption, and its AUC is the
  probability that a random responder sample has a higher frequency than a
  random non-responder sample (0.5 = no separation);
* the same test **at each timepoint**, to show when a difference appears. Only a
  baseline (day 0) difference could predict response before treatment starts.
"""

from __future__ import annotations

import sqlite3
import warnings
from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.tools.sm_exceptions import ConvergenceWarning

ALPHA = 0.05
RESPONDER = "yes"
NON_RESPONDER = "no"

RESULT_COLUMNS = [
    "population",
    "n_responders",
    "n_non_responders",
    "median_responders",
    "median_non_responders",
    "lmm_difference",
    "lmm_ci_low",
    "lmm_ci_high",
    "lmm_p",
    "lmm_q",
    "mwu_p",
    "mwu_q",
    "auc",
    "patient_variance",
    "significant",
]


@dataclass(frozen=True)
class Cohort:
    """The patients whose responders and non-responders are compared."""

    condition: str = "melanoma"
    treatment: str = "miraclib"
    sample_type: str = "PBMC"

    def describe(self) -> str:
        return f"{self.condition} patients treated with {self.treatment}, {self.sample_type} samples"


DEFAULT_COHORT = Cohort()


def format_p(p: float) -> str:
    """Compact p/q-value: three decimals, or '< 0.001'."""
    return "< 0.001" if p < 0.001 else f"{p:.3f}"


def cohort_frequencies(conn: sqlite3.Connection, cohort: Cohort = DEFAULT_COHORT) -> pd.DataFrame:
    """Sample-level relative frequencies for the cohort's responders and non-responders.

    One row per (sample, population) with columns ``sample``, ``subject_id``,
    ``subject``, ``response`` ('yes'/'no'), ``timepoint`` (days), ``population``
    and ``percentage``. Samples without a response value are excluded.
    """
    return pd.read_sql_query(
        """
        SELECT f.sample, m.subject_id, m.subject, m.response,
               m.time_from_treatment_start AS timepoint,
               f.population, f.percentage
        FROM cell_frequencies AS f
        JOIN sample_metadata AS m ON m.sample_id = f.sample_id
        WHERE m.condition = :condition
          AND m.treatment = :treatment
          AND m.sample_type = :sample_type
          AND m.response IN ('yes', 'no')
        ORDER BY f.sample_id, f.population_id
        """,
        conn,
        params=asdict(cohort),
    )


def cohort_counts(freqs: pd.DataFrame) -> pd.DataFrame:
    """Patients and samples per response group."""
    samples = freqs.drop_duplicates("sample")
    counts = samples.groupby("response").agg(
        n_subjects=("subject_id", "nunique"), n_samples=("sample", "size")
    )
    return counts.reindex([RESPONDER, NON_RESPONDER], fill_value=0).rename_axis("response").reset_index()


def _check_groups(freqs: pd.DataFrame) -> None:
    present = set(freqs["response"].unique())
    missing = {RESPONDER, NON_RESPONDER} - present
    if missing:
        names = " and ".join(sorted("responders" if m == RESPONDER else "non-responders" for m in missing))
        raise ValueError(f"Cannot compare groups: the cohort has no {names}.")


def _groups(d: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    responders = d.loc[d["response"] == RESPONDER, "percentage"].to_numpy(dtype=float)
    non_responders = d.loc[d["response"] == NON_RESPONDER, "percentage"].to_numpy(dtype=float)
    return responders, non_responders


def _mann_whitney(responders: np.ndarray, non_responders: np.ndarray) -> dict[str, float]:
    test = stats.mannwhitneyu(responders, non_responders, alternative="two-sided")
    return {
        "n_responders": responders.size,
        "n_non_responders": non_responders.size,
        "median_responders": float(np.median(responders)),
        "median_non_responders": float(np.median(non_responders)),
        "mwu_p": float(test.pvalue),
        # U counts (responder > non-responder) pairs, ties as 1/2, so U / (n1 * n2) is the AUC.
        "auc": float(test.statistic) / (responders.size * non_responders.size),
    }


# Optimisers tried in order. When a patient's samples are no more alike than any
# two samples, the patient variance is estimated at 0 - the boundary of the
# parameter space - and gradient-based optimisers can hit a singular matrix there,
# so derivative-free ones go first. All of them reach the same optimum.
_OPTIMIZERS = ("powell", "nm", "bfgs")


def _fit_mixed_model(model):
    problems = []
    for method in _OPTIMIZERS:
        with warnings.catch_warnings():
            # At the boundary statsmodels warns that the MLE "may be on the boundary"
            # and that the random-effects covariance is singular. The fit is still
            # valid (it reduces to ordinary least squares); convergence is checked below.
            warnings.simplefilter("ignore", ConvergenceWarning)
            warnings.filterwarnings("ignore", message="Random effects covariance is singular")
            try:
                fit = model.fit(reml=True, method=[method])
            except np.linalg.LinAlgError as exc:
                problems.append(f"{method}: {exc}")
                continue
        if fit.converged:
            return fit
        problems.append(f"{method}: did not converge")
    raise RuntimeError(f"Mixed model could not be fitted ({'; '.join(problems)})")


def _mixed_model(d: pd.DataFrame) -> dict[str, float]:
    data = d.assign(responder=(d["response"] == RESPONDER).astype(int))
    formula = "percentage ~ responder"
    if data["timepoint"].nunique() > 1:
        formula += " + C(timepoint)"
    fit = _fit_mixed_model(smf.mixedlm(formula, data, groups=data["subject_id"]))
    ci_low, ci_high = fit.conf_int().loc["responder"]
    return {
        "lmm_difference": float(fit.params["responder"]),
        "lmm_ci_low": float(ci_low),
        "lmm_ci_high": float(ci_high),
        "lmm_p": float(fit.pvalues["responder"]),
        # Between-patient variance (%²). Near the boundary its last digits are optimiser
        # noise, so it is rounded to keep the outputs identical from machine to machine.
        "patient_variance": round(float(fit.cov_re.iloc[0, 0]), 6),
    }


def _fdr(p_values: pd.Series) -> np.ndarray:
    return stats.false_discovery_control(p_values.to_numpy(dtype=float), method="bh")


def _sort_populations(df: pd.DataFrame, order: list[str], *extra: str) -> pd.DataFrame:
    rank = {name: i for i, name in enumerate(order)}
    return df.sort_values(
        ["population", *extra],
        key=lambda col: col.map(rank) if col.name == "population" else col,
        kind="stable",
    ).reset_index(drop=True)


def compare_responders(freqs: pd.DataFrame, alpha: float = ALPHA) -> pd.DataFrame:
    """Primary analysis: one row per population (see the module docstring).

    ``lmm_*`` columns come from the mixed model (difference in percentage points
    with its 95% CI, p and FDR q); ``mwu_*`` and ``auc`` from the Mann-Whitney U
    check. ``significant`` is ``lmm_q < alpha``.
    """
    _check_groups(freqs)
    rows = []
    for population, d in freqs.groupby("population", sort=False):
        responders, non_responders = _groups(d)
        rows.append({"population": population, **_mann_whitney(responders, non_responders), **_mixed_model(d)})
    results = pd.DataFrame(rows)
    results["lmm_q"] = _fdr(results["lmm_p"])
    results["mwu_q"] = _fdr(results["mwu_p"])
    results["significant"] = results["lmm_q"] < alpha
    return results[RESULT_COLUMNS]


def compare_by_timepoint(freqs: pd.DataFrame) -> pd.DataFrame:
    """Mann-Whitney U per population at each timepoint (exploratory).

    ``mwu_q`` is FDR-adjusted across all population x timepoint tests.
    """
    _check_groups(freqs)
    rows = []
    for (population, timepoint), d in freqs.groupby(["population", "timepoint"], sort=False):
        responders, non_responders = _groups(d)
        if responders.size and non_responders.size:
            rows.append({"population": population, "timepoint": int(timepoint), **_mann_whitney(responders, non_responders)})
    out = pd.DataFrame(rows)
    out["mwu_q"] = _fdr(out["mwu_p"])
    return _sort_populations(out, list(dict.fromkeys(freqs["population"])), "timepoint")


def timepoint_means(freqs: pd.DataFrame, confidence: float = 0.95) -> pd.DataFrame:
    """Mean relative frequency with a t-based confidence interval, per population x timepoint x group."""
    grouped = freqs.groupby(["population", "timepoint", "response"], sort=False)["percentage"]
    out = grouped.agg(n="size", mean="mean", sd="std").reset_index()
    half_width = stats.t.ppf(0.5 + confidence / 2, out["n"] - 1) * out["sd"] / np.sqrt(out["n"])
    out["ci_low"] = out["mean"] - half_width
    out["ci_high"] = out["mean"] + half_width
    return _sort_populations(out, list(dict.fromkeys(freqs["population"])), "timepoint", "response")


def key_findings(results: pd.DataFrame, by_timepoint: pd.DataFrame, alpha: float = ALPHA) -> list[str]:
    """Plain-language conclusions, generated from the numbers so the text cannot drift from them."""
    findings: list[str] = []
    significant = results[results["significant"]]
    others = results[~results["significant"]]
    baseline_day = by_timepoint["timepoint"].min()
    baseline = by_timepoint[by_timepoint["timepoint"] == baseline_day]

    if significant.empty:
        closest = results.loc[results["lmm_q"].idxmin()]
        findings.append(
            "No population differs significantly between responders and non-responders after FDR "
            f"correction (closest: {closest.population}, q = {format_p(closest.lmm_q)})."
        )
    for row in significant.itertuples():
        direction = "higher" if row.lmm_difference > 0 else "lower"
        findings.append(
            f"**{row.population}** differs significantly: its relative frequency is {direction} in "
            f"responders by {abs(row.lmm_difference):.2f} percentage points (95% CI "
            f"{row.lmm_ci_low:.2f} to {row.lmm_ci_high:.2f}; mixed model p = {format_p(row.lmm_p)}, "
            f"FDR q = {format_p(row.lmm_q)})."
        )
        same_direction = (row.auc - 0.5) * row.lmm_difference > 0
        findings.append(
            f"The distribution-free Mann–Whitney U test "
            f"{'agrees on the direction' if same_direction else 'disagrees on the direction'} "
            f"(p = {format_p(row.mwu_p)}; FDR q = {format_p(row.mwu_q)}). The effect is small, though: "
            f"AUC = {row.auc:.2f} (0.5 = no separation, 1 = perfect), so {row.population} alone would "
            "be a weak predictor of response."
        )
        timeline = by_timepoint[by_timepoint["population"] == row.population]
        at_baseline = timeline[timeline["timepoint"] == baseline_day].iloc[0]
        on_treatment = timeline[timeline["timepoint"] != baseline_day]
        later = ", ".join(f"day {t.timepoint} p = {format_p(t.mwu_p)}" for t in on_treatment.itertuples())
        if at_baseline.mwu_p < alpha:
            findings.append(
                f"The {row.population} difference is already present at baseline (day {baseline_day} "
                f"p = {format_p(at_baseline.mwu_p)}), so pre-treatment {row.population} frequency could "
                "help predict response."
            )
        elif (on_treatment["mwu_p"] < alpha).any():
            findings.append(
                f"The {row.population} difference is absent at baseline (day {baseline_day} "
                f"p = {format_p(at_baseline.mwu_p)}) and appears on treatment (unadjusted {later}), so it behaves "
                "like an early on-treatment marker of response rather than a pre-treatment predictor."
            )
        else:
            findings.append(
                f"No single timepoint shows a significant {row.population} difference (day {baseline_day} "
                f"p = {format_p(at_baseline.mwu_p)}, {later}); it only emerges when timepoints are pooled."
            )
    if not others.empty and not significant.empty:
        findings.append(
            f"None of the other populations is significant: {', '.join(others['population'])} have "
            f"mixed-model q ≥ {format_p(others['lmm_q'].min())} (smallest Mann–Whitney "
            f"p = {format_p(others['mwu_p'].min())})."
        )
    isolated = by_timepoint[(by_timepoint["mwu_p"] < alpha) & ~by_timepoint["population"].isin(significant["population"])]
    for survives, group in isolated.groupby(isolated["mwu_q"] < alpha):
        items = ", ".join(f"{t.population} on day {t.timepoint}, p = {format_p(t.mwu_p)}" for t in group.itertuples())
        if survives:
            findings.append(
                f"Single-timepoint differences that survive FDR correction across the {len(by_timepoint)} "
                f"timepoint tests: {items} (q ≤ {format_p(group['mwu_q'].max())}); worth following up."
            )
        else:
            findings.append(
                f"Isolated single-timepoint differences ({items}) do not survive FDR correction across the "
                f"{len(by_timepoint)} timepoint tests (q ≥ {format_p(group['mwu_q'].min())}); treat them as "
                "hypotheses to follow up, not findings."
            )
    if (baseline["mwu_p"] >= alpha).all():
        findings.append(
            f"At baseline (day {baseline_day}) no population differs between the groups (all Mann–Whitney "
            f"p ≥ {format_p(baseline['mwu_p'].min())}), so pre-treatment frequencies of these "
            f"{len(results)} populations do not predict response."
        )
    return findings
