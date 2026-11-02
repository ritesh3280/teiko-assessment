"""End to end: the pipeline outputs and the dashboard, from a freshly built database."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from streamlit.testing.v1 import AppTest

import run_analysis

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = [
    "part2_cell_frequencies.csv",
    "part2_population_summary.csv",
    "part3_responder_statistics.csv",
    "part3_timepoint_statistics.csv",
    "part3_boxplots.png",
    "part4_baseline_samples.csv",
    "part4_samples_per_project.csv",
    "part4_subjects_by_response.csv",
    "part4_subjects_by_sex.csv",
    "part4_avg_b_cells_melanoma_male_responders.csv",
    "report.md",
]


def test_pipeline_writes_every_output(db_path, tmp_path):
    run_analysis.run(db_path, tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(OUTPUTS)

    frequencies = pd.read_csv(tmp_path / "part2_cell_frequencies.csv")
    assert list(frequencies.columns) == ["sample", "total_count", "population", "count", "percentage"]
    assert len(frequencies) == 52_500

    statistics = pd.read_csv(tmp_path / "part3_responder_statistics.csv")
    assert statistics.loc[statistics["significant"], "population"].tolist() == ["cd4_t_cell"]

    b_cells = pd.read_csv(tmp_path / "part4_avg_b_cells_melanoma_male_responders.csv")
    assert b_cells.loc[0, "mean_count"] == 10206.15
    assert b_cells.loc[0, "n_samples"] == 485

    report = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "**cd4_t_cell** shows a statistically significant association" in report
    assert "not significant after FDR correction (q = 0.067)" in report
    assert "| prj1 | 384 |" in report and "| prj2 | 0 |" in report
    assert "**10206.15**" in report
    assert (tmp_path / "part3_boxplots.png").stat().st_size > 50_000


def test_dashboard_renders_without_errors(db_path, monkeypatch):
    monkeypatch.setenv("CELL_COUNTS_DB", str(db_path))
    app = AppTest.from_file(str(REPO_ROOT / "dashboard.py"), default_timeout=300).run()
    assert not app.exception
    assert [tab.label for tab in app.tabs] == [
        "Part 2 · Cell frequencies",
        "Part 3 · Responders vs non-responders",
        "Part 4 · Baseline subset",
    ]
    metrics = {metric.label: metric.value for metric in app.metric}
    assert metrics["Samples"] == "10,500"
    assert metrics["Significant in the primary model (FDR q < 0.05)"] == "CD4 T cells"
    assert metrics["Matching samples"] == "656"
    assert metrics["Average B cells (count)"] == "10206.15"
