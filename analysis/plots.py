"""Static figures written by the pipeline (the dashboard draws interactive Plotly versions)."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.patches import Patch

from analysis.responders import NON_RESPONDER, RESPONDER, format_p
from analysis.theme import LIGHT, RESPONSE_LABELS, population_label


# DejaVu Sans ships with matplotlib, so the PNG renders identically on every machine.
_RC = {"font.family": "DejaVu Sans"}


def responder_boxplots(freqs: pd.DataFrame, results: pd.DataFrame, path: Path) -> Path:
    """Part 3 figure: one boxplot panel per population, responders vs non-responders.

    Each panel has its own y-axis so the within-population comparison is legible;
    the line above each panel gives the mixed-model FDR q-value.
    """
    with matplotlib.rc_context(_RC):
        return _responder_boxplots(freqs, results, Path(path))


def _responder_boxplots(freqs: pd.DataFrame, results: pd.DataFrame, path: Path) -> Path:
    c = LIGHT
    groups = [(RESPONDER, c["responders"]), (NON_RESPONDER, c["non_responders"])]
    populations = list(dict.fromkeys(freqs["population"]))
    by_population = results.set_index("population")

    fig = Figure(figsize=(2.6 * len(populations), 5.0), layout="constrained", facecolor=c["surface"])
    axes = np.atleast_1d(fig.subplots(1, len(populations)))
    for ax, population in zip(axes, populations):
        d = freqs[freqs["population"] == population]
        data = [d.loc[d["response"] == code, "percentage"].to_numpy() for code, _ in groups]
        box = ax.boxplot(
            data,
            widths=0.5,
            patch_artist=True,
            tick_labels=["Responders", "Non-\nresponders"],
        )
        for i, (_, color) in enumerate(groups):
            box["boxes"][i].set(facecolor=to_rgba(color, 0.14), edgecolor=color, linewidth=1.5)
            box["medians"][i].set(color=color, linewidth=2.5)
            for line in [*box["whiskers"][2 * i : 2 * i + 2], *box["caps"][2 * i : 2 * i + 2]]:
                line.set(color=color, linewidth=1.5)
            box["fliers"][i].set(
                marker="o", markersize=3.5, markerfacecolor=to_rgba(color, 0.45), markeredgecolor="none"
            )

        ax.set_facecolor(c["surface"])
        ax.grid(axis="y", color=c["grid"], linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(c["axis"])
        ax.tick_params(colors=c["ink_secondary"], labelsize=9, length=0, pad=6)

        row = by_population.loc[population]
        ax.set_title(population_label(population), pad=22, fontsize=11.5, fontweight="semibold", color=c["ink"])
        verdict = "significant" if row.significant else "not significant"
        ax.text(
            0.5,
            1.025,
            f"q = {format_p(row.lmm_q)} · {verdict}",
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=9,
            color=c["ink"] if row.significant else c["ink_secondary"],
            fontweight="semibold" if row.significant else "normal",
        )
    axes[0].set_ylabel("Relative frequency (% of cells in sample)", color=c["ink_secondary"], fontsize=9.5)

    samples = freqs.drop_duplicates("sample")
    n_samples = samples["response"].value_counts()
    n_patients = samples.groupby("response")["subject_id"].nunique()
    fig.suptitle(
        "Immune population frequencies: responders vs non-responders",
        x=0.01,
        ha="left",
        fontsize=13,
        fontweight="semibold",
        color=c["ink"],
    )
    fig.supxlabel(
        f"Melanoma patients on miraclib, PBMC samples from days 0, 7 and 14: {n_samples[RESPONDER]:,} samples "
        f"from {n_patients[RESPONDER]:,} responders and {n_samples[NON_RESPONDER]:,} samples from "
        f"{n_patients[NON_RESPONDER]:,} non-responders.\n"
        "Box = median and interquartile range (IQR); whiskers = 1.5 × IQR; dots = outlying samples. "
        "q = Benjamini–Hochberg-adjusted p-value of a mixed-effects model with a patient random intercept.",
        x=0.01,
        ha="left",
        fontsize=8.5,
        color=c["ink_secondary"],
        linespacing=1.6,
    )
    fig.legend(
        handles=[
            Patch(facecolor=to_rgba(color, 0.14), edgecolor=color, linewidth=1.5, label=RESPONSE_LABELS[code])
            for code, color in groups
        ],
        loc="outside upper right",
        ncols=2,
        frameon=False,
        fontsize=9.5,
        labelcolor=c["ink_secondary"],
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, facecolor=c["surface"])
    return path
