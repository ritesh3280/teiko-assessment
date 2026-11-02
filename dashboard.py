"""Interactive dashboard for the miraclib trial analysis (Parts 2-4).

Start it with ``make dashboard`` (or ``streamlit run dashboard.py``) and open
http://localhost:8501. It reads the SQLite database built by ``make pipeline``;
if that database does not exist yet, it is built from ``cell-count.csv`` first.
All numbers come from the same ``analysis`` functions the pipeline uses.
"""

from __future__ import annotations

import os
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.colors import hex_to_rgb

import load_data
from analysis import frequencies, responders, subsets
from analysis.db import open_db
from analysis.responders import NON_RESPONDER, RESPONDER, format_p
from analysis.theme import DARK, LIGHT, RESPONSE_LABELS, SEX_LABELS, population_label

DB_PATH = Path(os.environ.get("CELL_COUNTS_DB", load_data.DB_PATH))  # override used by the tests
GROUPS = [(RESPONDER, "responders"), (NON_RESPONDER, "non_responders")]  # fixed order = fixed colours
CHART_CONFIG = {"displaylogo": False, "modeBarButtonsToRemove": ["select2d", "lasso2d", "autoScale2d"]}

st.set_page_config(page_title="Miraclib immune profiling", page_icon=":material/biotech:", layout="wide")


# ----------------------------------------------------------------------------- data


def database_version() -> int:
    """Build the database if it is missing; return its mtime so caches refresh after a rebuild."""
    if not DB_PATH.exists():
        with st.spinner("Building the database from cell-count.csv…"):
            load_data.build_database(load_data.CSV_PATH, DB_PATH)
    return DB_PATH.stat().st_mtime_ns


@st.cache_data(show_spinner=False)
def load_overview(version: int) -> dict[str, int]:
    with open_db(DB_PATH) as conn:
        return load_data.table_counts(conn)


@st.cache_data(show_spinner="Loading cell frequencies…")
def load_frequency_table(version: int) -> pd.DataFrame:
    with open_db(DB_PATH) as conn:
        return frequencies.cell_frequency_table(conn)


@st.cache_data(show_spinner="Fitting the responder models…")
def load_responder_analysis(version: int) -> dict:
    with open_db(DB_PATH) as conn:
        freqs = responders.cohort_frequencies(conn)
    results = responders.compare_responders(freqs)
    by_timepoint = responders.compare_by_timepoint(freqs)
    return {
        "freqs": freqs,
        "results": results,
        "by_timepoint": by_timepoint,
        "means": responders.timepoint_means(freqs),
        "counts": responders.cohort_counts(freqs).set_index("response"),
        "findings": responders.key_findings(results, by_timepoint),
    }


@st.cache_data(show_spinner=False)
def load_filter_options(version: int) -> dict[str, list]:
    with open_db(DB_PATH) as conn:
        return subsets.filter_options(conn)


@st.cache_data(show_spinner=False)
def load_b_cell_answer(version: int) -> pd.DataFrame:
    with open_db(DB_PATH) as conn:
        return subsets.mean_cell_count(conn, subsets.B_CELL_QUESTION)


@st.cache_data(show_spinner=False)
def load_subset(version: int, condition: str, treatment: str, sample_type: str, timepoint: int) -> dict:
    sample_filter = subsets.SampleFilter(condition, treatment, sample_type, timepoint)
    with open_db(DB_PATH) as conn:
        return {
            "samples": subsets.subset_samples(conn, sample_filter),
            "per_project": subsets.samples_per_project(conn, sample_filter),
            "by_response": subsets.subjects_by_response(conn, sample_filter),
            "by_sex": subsets.subjects_by_sex(conn, sample_filter),
        }


# --------------------------------------------------------------------------- charts


def palette() -> dict[str, str]:
    """Chart colours for the viewer's theme (light when unknown)."""
    return DARK if st.context.theme.type == "dark" else LIGHT


def rgba(hex_color: str, alpha: float) -> str:
    r, g, b = hex_to_rgb(hex_color)
    return f"rgba({r}, {g}, {b}, {alpha})"


def legend(colors: dict[str, str]) -> None:
    """One legend for a row of small-multiple charts."""
    items = "".join(
        f'<span style="display:inline-flex;align-items:center;gap:0.4rem;margin-right:1.25rem">'
        f'<span style="width:0.8rem;height:0.8rem;border-radius:3px;border:1.5px solid {colors[key]};'
        f'background:{rgba(colors[key], 0.25)}"></span>{RESPONSE_LABELS[code]}</span>'
        for code, key in GROUPS
    )
    st.markdown(f'<div style="font-size:0.9rem;margin:0.25rem 0 0.5rem">{items}</div>', unsafe_allow_html=True)


def chart_title(title: str, subtitle: str, colors: dict[str, str], emphasise: bool = False) -> dict:
    weight = "bold" if emphasise else "normal"
    color = colors["ink"] if emphasise else colors["ink_secondary"]
    return {
        "text": f"<b>{title}</b><br><span style='font-size:12px;font-weight:{weight};color:{color}'>{subtitle}</span>",
        "x": 0.5,
        "xanchor": "center",
        "font": {"size": 15},
    }


def box_figure(d: pd.DataFrame, title: dict, colors: dict[str, str], show_y_title: bool) -> go.Figure:
    fig = go.Figure()
    for code, key in GROUPS:
        values = d.loc[d["response"] == code, "percentage"]
        fig.add_trace(
            go.Box(
                y=values,
                x=[RESPONSE_LABELS[code]] * len(values),
                name=RESPONSE_LABELS[code],
                line={"color": colors[key], "width": 1.5},
                fillcolor=rgba(colors[key], 0.14),
                marker={"color": colors[key], "size": 4, "opacity": 0.55},
                boxpoints="outliers",
                width=0.55,
                hoverlabel={"namelength": -1},
            )
        )
    fig.update_layout(
        title=title,
        height=380,
        margin={"l": 8, "r": 8, "t": 72, "b": 8},
        showlegend=False,
        yaxis_title="Relative frequency (%)" if show_y_title else None,
    )
    fig.update_xaxes(showgrid=False)
    return fig


def trajectory_figure(m: pd.DataFrame, title: dict, colors: dict[str, str], show_y_title: bool) -> go.Figure:
    fig = go.Figure()
    for offset, (code, key) in zip((-0.8, 0.8), GROUPS):
        g = m[m["response"] == code]
        fig.add_trace(
            go.Scatter(
                x=g["timepoint"] + offset,  # dodge the two groups so their intervals don't overlap
                y=g["mean"],
                customdata=g[["timepoint", "ci_low", "ci_high", "n"]].to_numpy(),
                mode="lines+markers",
                name=RESPONSE_LABELS[code],
                line={"color": colors[key], "width": 2},
                marker={"color": colors[key], "size": 9, "line": {"color": colors["surface"], "width": 2}},
                error_y={
                    "type": "data",
                    "symmetric": False,
                    "array": g["ci_high"] - g["mean"],
                    "arrayminus": g["mean"] - g["ci_low"],
                    "color": colors[key],
                    "thickness": 1.5,
                    "width": 3,
                },
                hovertemplate=(
                    "<b>%{y:.2f}%</b> mean · day %{customdata[0]}<br>"
                    "95% CI %{customdata[1]:.2f}–%{customdata[2]:.2f} · n = %{customdata[3]}"
                    f"<extra>{RESPONSE_LABELS[code]}</extra>"
                ),
            )
        )
    days = sorted(m["timepoint"].unique())
    fig.update_layout(
        title=title,
        height=320,
        margin={"l": 8, "r": 8, "t": 72, "b": 8},
        showlegend=False,
        yaxis_title="Mean relative frequency (%)" if show_y_title else None,
    )
    fig.update_xaxes(tickvals=days, ticktext=[f"Day {d}" for d in days], range=[days[0] - 2, days[-1] + 2], showgrid=False)
    return fig


def sample_context_figure(freq: pd.DataFrame, sample: str, colors: dict[str, str]) -> go.Figure:
    """Each population's distribution across all samples (gray), with one sample highlighted."""
    order = list(dict.fromkeys(freq["population"]))
    labels = [population_label(p) for p in order]
    q = freq.groupby("population", sort=False)["percentage"].quantile([0.25, 0.5, 0.75]).unstack().loc[order]
    iqr = q[0.75] - q[0.25]
    grouped = freq.groupby("population", sort=False)["percentage"]
    lower = [grouped.get_group(p)[lambda s: s >= q.loc[p, 0.25] - 1.5 * iqr[p]].min() for p in order]
    upper = [grouped.get_group(p)[lambda s: s <= q.loc[p, 0.75] + 1.5 * iqr[p]].max() for p in order]
    picked = freq[freq["sample"] == sample].set_index("population").loc[order]

    fig = go.Figure()
    fig.add_trace(
        go.Box(
            y=labels,
            q1=q[0.25],
            median=q[0.5],
            q3=q[0.75],
            lowerfence=lower,
            upperfence=upper,
            orientation="h",
            name="All samples",
            line={"color": colors["ink_muted"], "width": 1.5},
            fillcolor=rgba(colors["ink_muted"], 0.12),
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=picked["percentage"],
            y=labels,
            mode="markers",
            name=sample,
            marker={"symbol": "diamond", "size": 13, "color": colors["highlight"], "line": {"color": colors["surface"], "width": 2}},
            customdata=picked[["count"]].to_numpy(),
            hovertemplate=f"<b>%{{x:.2f}}%</b> · %{{customdata[0]:,}} cells<extra>{sample}</extra>",
        )
    )
    fig.update_layout(
        height=330,
        margin={"l": 8, "r": 8, "t": 8, "b": 8},
        legend={"orientation": "h", "y": 1.08, "x": 0},
        xaxis_title="Relative frequency (% of cells in sample)",
    )
    fig.update_yaxes(autorange="reversed", showgrid=False)
    return fig


# -------------------------------------------------------------------------- tabs


def frequencies_tab(version: int) -> None:
    freq = load_frequency_table(version)
    populations = list(dict.fromkeys(freq["population"]))
    st.subheader("What fraction of each sample is each cell population?")
    st.caption(
        "For every sample, `total_count` is the sum of the five populations and `percentage` is "
        "`count / total_count × 100`. One row per sample and population."
    )

    search_col, population_col = st.columns([1, 2])
    query = search_col.text_input("Find samples", placeholder="Sample ID, e.g. sample00042")
    chosen = population_col.multiselect("Populations", populations, default=populations, format_func=population_label)
    view = freq[freq["population"].isin(chosen)]
    if query.strip():
        view = view[view["sample"].str.contains(query.strip(), case=False, regex=False)]

    event = st.dataframe(
        view,
        hide_index=True,
        height=420,
        on_select="rerun",
        selection_mode="single-row",
        key="frequency_table",
        column_config={
            "sample": st.column_config.TextColumn("sample"),
            "total_count": st.column_config.NumberColumn("total_count", format="%d"),
            "population": st.column_config.TextColumn("population"),
            "count": st.column_config.NumberColumn("count", format="%d"),
            "percentage": st.column_config.ProgressColumn("percentage", format="%.2f%%", min_value=0, max_value=100),
        },
    )
    info_col, download_col = st.columns([3, 1], vertical_alignment="center")
    info_col.caption(f"{len(view):,} of {len(freq):,} rows · select a row to see that sample in context below.")
    download_col.download_button(
        "Download CSV",
        view.to_csv(index=False).encode(),
        file_name="cell_frequencies.csv",
        mime="text/csv",
        icon=":material/download:",
        width="stretch",
    )

    if view.empty:
        st.info("No rows match the current filters.")
        return
    selected = [row for row in (event.selection.rows if event else []) if row < len(view)]
    sample = view.iloc[selected[0] if selected else 0]["sample"]
    st.markdown(f"#### Sample `{sample}` in context")
    st.caption("Gray boxes: distribution of each population across all 10,500 samples (median, IQR, 1.5 × IQR whiskers). Diamond: the selected sample.")
    st.plotly_chart(sample_context_figure(freq, sample, palette()), config=CHART_CONFIG, key="sample_context")


def responders_tab(version: int) -> None:
    data = load_responder_analysis(version)
    results, by_timepoint, counts = data["results"], data["by_timepoint"], data["counts"]
    colors = palette()
    populations = list(results["population"])

    st.subheader("Do cell population frequencies differ between responders and non-responders?")
    st.caption(f"Cohort: {responders.DEFAULT_COHORT.describe()}; every patient was sampled on days 0, 7 and 14.")
    significant = results.loc[results["significant"], "population"].tolist()
    k1, k2, k3 = st.columns(3)
    k1.metric("Responders", f"{counts.loc[RESPONDER, 'n_subjects']:,} patients", border=True,
              help=f"{counts.loc[RESPONDER, 'n_samples']:,} PBMC samples")
    k2.metric("Non-responders", f"{counts.loc[NON_RESPONDER, 'n_subjects']:,} patients", border=True,
              help=f"{counts.loc[NON_RESPONDER, 'n_samples']:,} PBMC samples")
    k3.metric(
        "Significant in the primary model (FDR q < 0.05)",
        ", ".join(map(population_label, significant)) or "None",
        border=True,
        help="Linear mixed-effects model. The secondary Mann–Whitney test is reported in the table below.",
    )

    with st.container(border=True):
        st.markdown("**What the statistics say**")
        st.markdown("\n".join(f"- {finding}" for finding in data["findings"]))

    st.markdown("#### Relative frequencies by response")
    views = ["All timepoints", *[f"Day {t}" for t in sorted(data["freqs"]["timepoint"].unique())]]
    view = st.segmented_control("Samples shown", views, default=views[0], required=True, key="box_view")
    legend(colors)
    freqs, day = data["freqs"], None
    if view != views[0]:
        day = int(view.split()[1])
        freqs = freqs[freqs["timepoint"] == day]
    by_population = results.set_index("population")
    day_p = by_timepoint.set_index(["population", "timepoint"])["mwu_p"]
    for i, (column, population) in enumerate(zip(st.columns(len(populations)), populations)):
        if day is None:  # pooled view: the primary analysis
            row = by_population.loc[population]
            subtitle = f"q = {format_p(row.lmm_q)} · {'significant' if row.significant else 'not significant'}"
            emphasise = bool(row.significant)
        else:
            p = day_p.loc[(population, day)]
            subtitle, emphasise = f"p = {format_p(p)} (unadjusted)", bool(p < responders.ALPHA)
        title = chart_title(population_label(population), subtitle, colors, emphasise)
        with column:
            fig = box_figure(freqs[freqs["population"] == population], title, colors, show_y_title=i == 0)
            st.plotly_chart(fig, config=CHART_CONFIG, key=f"box_{population}")
    st.caption(
        "Box = median and interquartile range; whiskers = 1.5 × IQR; dots = outlying samples. "
        "All timepoints: q is the Benjamini–Hochberg-adjusted p-value of the mixed-effects model. "
        "Single days: unadjusted two-sided Mann–Whitney U p-value."
    )

    st.markdown("#### Statistics per population")
    table = pd.DataFrame(
        {
            "population": results["population"],
            "median_responders": results["median_responders"],
            "median_non_responders": results["median_non_responders"],
            "difference_pp": results["lmm_difference"],
            "ci_95": [f"{lo:+.2f} to {hi:+.2f}" for lo, hi in zip(results["lmm_ci_low"], results["lmm_ci_high"])],
            "lmm_p": results["lmm_p"],
            "lmm_q": results["lmm_q"],
            "mwu_p": results["mwu_p"],
            "mwu_q": results["mwu_q"],
            "auc": results["auc"],
            "significant": results["significant"],
        }
    )
    st.dataframe(
        table,
        hide_index=True,
        column_config={
            "population": "Population",
            "median_responders": st.column_config.NumberColumn("Median % (responders)", format="%.2f"),
            "median_non_responders": st.column_config.NumberColumn("Median % (non-responders)", format="%.2f"),
            "difference_pp": st.column_config.NumberColumn(
                "Difference (pp)", format="%+.2f", help="Mixed-model mean difference, responders − non-responders, in percentage points"
            ),
            "ci_95": st.column_config.TextColumn("95% CI (pp)"),
            "lmm_p": st.column_config.NumberColumn("Mixed-model p", format="%.4f"),
            "lmm_q": st.column_config.NumberColumn("q (FDR)", format="%.4f", help="Benjamini–Hochberg across the 5 populations"),
            "mwu_p": st.column_config.NumberColumn("Mann–Whitney p", format="%.4f"),
            "mwu_q": st.column_config.NumberColumn("Mann–Whitney q", format="%.4f"),
            "auc": st.column_config.NumberColumn(
                "AUC", format="%.3f", help="P(random responder sample > random non-responder sample); 0.5 = no separation"
            ),
            "significant": st.column_config.CheckboxColumn("Significant (primary)", help="Mixed-model q < 0.05"),
        },
    )

    st.markdown("#### When does a difference appear?")
    st.caption(
        "Mean relative frequency at each timepoint with 95% confidence intervals. "
        "A marker that could predict response before treatment would already differ on day 0."
    )
    legend(colors)
    means = data["means"]
    tp = by_timepoint.set_index(["population", "timepoint"])
    for i, (column, population) in enumerate(zip(st.columns(len(populations)), populations)):
        pvals = " · ".join(f"d{t} {format_p(tp.loc[(population, t), 'mwu_p'])}" for t in sorted(means["timepoint"].unique()))
        title = chart_title(population_label(population), f"p: {pvals}", colors)
        with column:
            fig = trajectory_figure(means[means["population"] == population], title, colors, show_y_title=i == 0)
            st.plotly_chart(fig, config=CHART_CONFIG, key=f"trajectory_{population}")
    with st.expander("Per-timepoint test results (table)"):
        st.dataframe(
            by_timepoint,
            hide_index=True,
            column_config={
                "population": "Population",
                "timepoint": st.column_config.NumberColumn("Day"),
                "n_responders": "n (responders)",
                "n_non_responders": "n (non-responders)",
                "median_responders": st.column_config.NumberColumn("Median % (responders)", format="%.2f"),
                "median_non_responders": st.column_config.NumberColumn("Median % (non-responders)", format="%.2f"),
                "mwu_p": st.column_config.NumberColumn("Mann–Whitney p", format="%.4f"),
                "auc": st.column_config.NumberColumn("AUC", format="%.3f"),
                "mwu_q": st.column_config.NumberColumn("q (FDR, 15 tests)", format="%.4f"),
            },
        )

    with st.expander("How the statistics work"):
        st.markdown(
            """
- **Primary test – linear mixed-effects model.** Each patient contributes three samples (days 0, 7, 14),
  so samples are not independent. For each population we fit
  `percentage ~ responder + C(timepoint) + (1 | patient)`: the `responder` coefficient is the difference in
  mean relative frequency (percentage points) between responders and non-responders, adjusted for timepoint,
  and the patient random intercept accounts for repeated sampling.
- **Multiple testing.** Five populations are tested, so p-values are adjusted with the Benjamini–Hochberg
  false discovery rate; a population is called significant when **q < 0.05**.
- **Secondary check – Mann–Whitney U.** A two-sided rank test on the sample-level frequencies shown in the
  boxplots; it does not assume normality and is FDR-adjusted the same way. Its **AUC** is the probability that
  a random responder sample has a higher frequency than a random non-responder sample (0.5 = no separation,
  1 = perfect), a direct measure of how useful the population would be for predicting response.
- **Reading the two together.** The tests measure different things (a difference in mean frequency vs. a
  tendency for one group's values to rank higher), so they need not agree. When only the primary test
  survives correction, as for CD4 T cells, read the evidence as moderate rather than conclusive.
- **Timing.** The same rank test at each timepoint (exploratory, FDR-adjusted across all 15 tests in the table)
  shows whether a difference exists before treatment (day 0) or only appears on treatment.
- **Note.** The estimated between-patient variance is ≈ 0 for every population: a patient's samples are no
  more alike than samples from different patients, so accounting for repeated sampling barely changes the
  estimates.
"""
        )


def b_cell_question(version: int) -> None:
    """The fixed follow-up question. It spans every treatment and sample type, so it sits
    above the subset filters (which scope everything below them) and ignores them."""
    question = subsets.B_CELL_QUESTION
    answer = load_b_cell_answer(version).iloc[0]
    with st.container(border=True):
        st.markdown(
            "**Average B cells in male melanoma responders at baseline.** Considering melanoma males of all "
            "sample and treatment types, what is the average number of B cells for responders at time = 0?"
        )
        value_col, detail_col = st.columns([1, 3], vertical_alignment="center")
        value_col.metric("Average B cells (count)", f"{answer['mean_count']:.2f}")
        detail_col.caption(
            f"Mean B-cell count of {answer['n_samples']:,} samples from {answer['n_subjects']:,} subjects: "
            f"{question.describe()}. Not affected by the filters below."
        )
        with st.expander("SQL behind this number"):
            st.code(subsets.MEAN_COUNT_SQL, language="sql")
            st.caption("Parameters: " + ", ".join(f"`{k}` = `{v!r}`" for k, v in asdict(question).items()))


def subset_tab(version: int) -> None:
    b_cell_question(version)
    options = load_filter_options(version)
    defaults = subsets.BASELINE
    keys = {"condition": "p4_condition", "treatment": "p4_treatment", "sample_type": "p4_sample_type", "timepoint": "p4_timepoint"}

    def reset() -> None:
        for field, key in keys.items():
            st.session_state[key] = getattr(defaults, field)

    for field, key in keys.items():
        if key not in st.session_state and getattr(defaults, field) in options[field]:
            st.session_state[key] = getattr(defaults, field)

    heading = st.empty()  # filled in once the filter values are known
    st.caption(
        "The filters default to Bob's question (melanoma · miraclib · PBMC · day 0). "
        "Change them to explore other subsets; every number is recomputed from the database."
    )
    c1, c2, c3, c4, c5 = st.columns([1, 1, 1, 1, 0.8], vertical_alignment="bottom")
    condition = c1.selectbox("Condition", options["condition"], key=keys["condition"])
    treatment = c2.selectbox("Treatment", options["treatment"], key=keys["treatment"])
    sample_type = c3.selectbox("Sample type", options["sample_type"], key=keys["sample_type"])
    timepoint = c4.selectbox(
        "Time from treatment start", options["timepoint"], key=keys["timepoint"], format_func=lambda d: f"Day {d}"
    )
    c5.button("Reset to Bob's question", on_click=reset, icon=":material/restart_alt:", width="stretch")

    who = "untreated patients" if treatment == "none" else f"{treatment}-treated patients"
    when = "at baseline (day 0)" if timepoint == 0 else f"on day {timepoint}"
    heading.subheader(f"{condition.capitalize()} {sample_type} samples {when} from {who}")

    data = load_subset(version, condition, treatment, sample_type, int(timepoint))
    samples = data["samples"]
    n_subjects = len(samples[["project", "subject"]].drop_duplicates())
    if samples.empty:
        st.info("No samples match these filters.")
        return

    k1, k2, k3 = st.columns(3)
    k1.metric("Matching samples", f"{len(samples):,}", border=True)
    k2.metric("Subjects", f"{n_subjects:,}", border=True)
    k3.metric("Projects with samples", f"{(data['per_project']['n_samples'] > 0).sum()} of {len(data['per_project'])}", border=True)

    def share(n: int, total: int, unit: str) -> str:
        return f"{n / total:.1%} of {unit}" if total else ""

    project_col, response_col, sex_col = st.columns(3)
    with project_col.container(border=True):
        st.markdown("**Samples per project**")
        for row in data["per_project"].itertuples():
            st.metric(row.project, f"{row.n_samples:,}")
            st.caption(share(row.n_samples, len(samples), "samples"))
    with response_col.container(border=True):
        st.markdown("**Subjects by response**")
        for row in data["by_response"].itertuples():
            st.metric(RESPONSE_LABELS.get(row.response, "Not applicable"), f"{row.n_subjects:,}")
            st.caption(share(row.n_subjects, n_subjects, "subjects"))
    with sex_col.container(border=True):
        st.markdown("**Subjects by sex**")
        for row in data["by_sex"].itertuples():
            st.metric(SEX_LABELS.get(row.sex, "Unknown"), f"{row.n_subjects:,}")
            st.caption(share(row.n_subjects, n_subjects, "subjects"))

    st.markdown("#### Matching samples")
    st.dataframe(samples, hide_index=True, height=360)
    st.download_button(
        "Download CSV",
        samples.to_csv(index=False).encode(),
        file_name="subset_samples.csv",
        mime="text/csv",
        icon=":material/download:",
    )
    with st.expander("SQL behind these numbers"):
        st.caption("Each query starts from the same filtered `subset` of samples; parameters are bound from the filters above.")
        for name in ("samples_per_project", "subjects_by_response", "subjects_by_sex"):
            st.code(subsets.sql(name), language="sql")


# --------------------------------------------------------------------------- page


def main() -> None:
    try:
        version = database_version()
    except (OSError, ValueError) as exc:
        st.error(f"Could not build the database: {exc}")
        st.stop()
    overview = load_overview(version)

    st.title("Miraclib immune profiling")
    st.caption("Loblaw Bio · immune cell populations in the miraclib clinical trial")
    k1, k2, k3, k4 = st.columns(4)
    k1.metric("Projects", f"{overview['projects']:,}")
    k2.metric("Subjects", f"{overview['subjects']:,}")
    k3.metric("Samples", f"{overview['samples']:,}")
    k4.metric("Cell populations", f"{overview['cell_populations']:,}")

    tab2, tab3, tab4 = st.tabs(
        ["Part 2 · Cell frequencies", "Part 3 · Responders vs non-responders", "Part 4 · Baseline subset"]
    )
    with tab2:
        frequencies_tab(version)
    with tab3:
        responders_tab(version)
    with tab4:
        subset_tab(version)

    built = datetime.fromtimestamp(version / 1e9).strftime("%Y-%m-%d %H:%M")
    st.divider()
    st.caption(f"Data: `{DB_PATH.name}` (built {built} by `load_data.py` from `{load_data.CSV_PATH.name}`).")


main()
