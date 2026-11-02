#!/usr/bin/env python3
"""Parts 2-4: analyse the database built by ``load_data.py`` and write ``outputs/``.

Usage (from the repository root)::

    python run_analysis.py        # ``make pipeline`` runs load_data.py first

Outputs
-------
part2_cell_frequencies.csv       Part 2 summary table (sample, total_count, population, count, percentage)
part2_population_summary.csv     distribution of each population's frequency across samples
part3_responder_statistics.csv   responders vs non-responders, one row per population
part3_timepoint_statistics.csv   the same comparison at each timepoint
part3_boxplots.png               boxplot per population, responders vs non-responders
part4_baseline_samples.csv       baseline melanoma PBMC samples from miraclib-treated patients
part4_samples_per_project.csv    ... how many per project
part4_subjects_by_response.csv   ... how many subjects were responders / non-responders
part4_subjects_by_sex.csv        ... how many subjects were male / female
report.md                        human-readable summary of all of the above
"""

from __future__ import annotations

import sys
from pathlib import Path

from analysis import frequencies, plots, report, responders, subsets
from analysis.db import open_db
from load_data import DB_PATH, REPO_ROOT

OUTPUT_DIR = REPO_ROOT / "outputs"


def run(db_path: Path = DB_PATH, output_dir: Path = OUTPUT_DIR) -> dict:
    """Run Parts 2-4 against ``db_path``, write every output into ``output_dir`` and return the results."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cohort, sample_filter = responders.DEFAULT_COHORT, subsets.BASELINE

    with open_db(db_path) as conn:
        frequency_table = frequencies.cell_frequency_table(conn)  # Part 2
        cohort_freqs = responders.cohort_frequencies(conn, cohort)  # Part 3
        subset = subsets.subset_samples(conn, sample_filter)  # Part 4
        per_project = subsets.samples_per_project(conn, sample_filter)
        by_response = subsets.subjects_by_response(conn, sample_filter)
        by_sex = subsets.subjects_by_sex(conn, sample_filter)

    population_summary = frequencies.population_summary(frequency_table)
    results = responders.compare_responders(cohort_freqs)
    by_timepoint = responders.compare_by_timepoint(cohort_freqs)
    findings = responders.key_findings(results, by_timepoint)

    def write_csv(df, name: str) -> None:
        df.to_csv(output_dir / name, index=False, float_format="%.6g")

    frequency_table.to_csv(output_dir / "part2_cell_frequencies.csv", index=False, float_format="%.4f")
    write_csv(population_summary, "part2_population_summary.csv")
    write_csv(results, "part3_responder_statistics.csv")
    write_csv(by_timepoint, "part3_timepoint_statistics.csv")
    plots.responder_boxplots(cohort_freqs, results, output_dir / "part3_boxplots.png")
    write_csv(subset, "part4_baseline_samples.csv")
    write_csv(per_project, "part4_samples_per_project.csv")
    write_csv(by_response, "part4_subjects_by_response.csv")
    write_csv(by_sex, "part4_subjects_by_sex.csv")
    markdown = report.render_report(
        n_samples=frequency_table["sample"].nunique(),
        population_summary=population_summary,
        cohort=cohort,
        cohort_counts=responders.cohort_counts(cohort_freqs),
        results=results,
        by_timepoint=by_timepoint,
        findings=findings,
        sample_filter=sample_filter,
        subset=subset,
        per_project=per_project,
        by_response=by_response,
        by_sex=by_sex,
    )
    (output_dir / "report.md").write_text(markdown, encoding="utf-8")

    return {
        "frequency_table": frequency_table,
        "results": results,
        "findings": findings,
        "subset": subset,
        "per_project": per_project,
        "by_response": by_response,
        "by_sex": by_sex,
    }


def main() -> int:
    try:
        out = run()
    except FileNotFoundError as exc:
        print(f"run_analysis.py: error: {exc}", file=sys.stderr)
        return 1

    rel = OUTPUT_DIR.relative_to(REPO_ROOT)
    significant = out["results"].loc[out["results"]["significant"], ["population", "lmm_q"]]
    names = ", ".join(f"{p} (q = {responders.format_p(q)})" for p, q in significant.itertuples(index=False))

    def counts(df, key: str, value: str) -> str:
        return ", ".join(f"{k}: {v}" for k, v in zip(df[key], df[value]))

    print(f"Part 2  {len(out['frequency_table']):,} rows -> {rel}/part2_cell_frequencies.csv")
    print(f"Part 3  significant populations (mixed model, FDR q < {responders.ALPHA}): {names or 'none'}")
    print(f"Part 4  {len(out['subset']):,} baseline samples | per project: {counts(out['per_project'], 'project', 'n_samples')}")
    print(f"        subjects by response: {counts(out['by_response'], 'response', 'n_subjects')}"
          f" | by sex: {counts(out['by_sex'], 'sex', 'n_subjects')}")
    print(f"Report  {rel}/report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
