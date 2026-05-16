"""Reporting utilities for LOXDockAudit."""

from __future__ import annotations

import math
from dataclasses import asdict
from pathlib import Path
from typing import Any, Sequence

import pandas as pd

from loxdockaudit.models import (
    ConstructConfig,
    ConstructSummary,
    ControlComparisonResult,
    PoseMetrics,
    ScreenResult,
)


POSE_COLUMNS = [
    "rank",
    "model_path",
    "active_site_to_target_distance",
    "productive",
    "lox_contacts",
    "cbd_contacts",
    "closest_active_site_resi",
    "warning",
]

CONTROL_COMPARISON_COLUMNS = [
    "construct_id",
    "control_type",
    "candidate_productive_count",
    "control_productive_count",
    "candidate_best_distance",
    "control_best_distance",
    "candidate_beats_control_count",
    "candidate_beats_control_distance",
    "conclusion",
]


def poses_to_dataframe(pose_metrics_list: list[PoseMetrics]) -> pd.DataFrame:
    """
    Convert list of PoseMetrics to a pandas DataFrame sorted by rank ascending.
    """
    rows = [asdict(metrics) for metrics in pose_metrics_list]
    frame = pd.DataFrame(rows, columns=POSE_COLUMNS)
    if frame.empty:
        return frame
    return frame.sort_values("rank", ascending=True).reset_index(drop=True)


def summaries_to_dataframe(summaries: list[ConstructSummary]) -> pd.DataFrame:
    """
    Convert list of ConstructSummary to a DataFrame with productive_fraction.
    """
    rows = [asdict(summary) for summary in summaries]
    frame = pd.DataFrame(rows)
    if frame.empty:
        frame["productive_fraction"] = pd.Series(dtype=float)
        return frame

    frame["productive_fraction"] = frame.apply(
        lambda row: (
            row["productive_count"] / row["total_poses"]
            if row["total_poses"]
            else 0.0
        ),
        axis=1,
    )
    return frame


def generate_markdown_report(
    construct_config: ConstructConfig,
    pose_df: pd.DataFrame,
    summary: ConstructSummary,
    audit_report: str,
    baseline_summary: ConstructSummary | None = None,
) -> str:
    """
    Generate a Markdown report for one construct.
    """
    decision = _executive_decision(summary)
    top_pose_table = _pose_table(pose_df)
    warnings = _collect_warnings(pose_df)

    lines = [
        f"# LOXDockAudit Report: {construct_config.construct_id}",
        "",
        "## Executive Decision",
        "",
        f"**{decision[0]}** - {decision[1]}",
        "",
        "## Input Audit",
        "",
        audit_report.strip() if audit_report.strip() else "No audit report provided.",
        "",
        "## Docking Pose Summary",
        "",
        top_pose_table,
        "",
        "## Productive Geometry",
        "",
        f"- Productive poses: {summary.productive_count} / {summary.total_poses}",
        f"- Best distance: {summary.best_distance:.3f} A",
        f"- Best productive rank: {_format_optional(summary.best_productive_rank)}",
        f"- Productive distance threshold: "
        f"{construct_config.productive_distance_threshold:.3f} A",
        f"- LOX contact frequency: {summary.lox_contact_frequency:.3f}",
        f"- CBD contact frequency: {summary.cbd_contact_frequency:.3f}",
        f"- CBD coupling: {summary.cbd_coupling:.3f}",
        "",
    ]

    if baseline_summary is not None:
        lines.extend(
            [
                "## Control / Baseline Comparison",
                "",
                _baseline_comparison_table(summary, baseline_summary),
                "",
            ]
        )

    lines.extend(
        [
            "## Warnings",
            "",
            _warnings_table(warnings),
        ]
    )

    return "\n".join(lines).rstrip() + "\n"


def save_outputs(
    out_dir: str,
    construct_id: str,
    pose_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    markdown_report: str,
) -> None:
    """
    Save pose CSV, summary CSV, and Markdown report to out_dir.
    """
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pose_df.to_csv(output_dir / f"{construct_id}_poses.csv", index=False)
    summary_df.to_csv(output_dir / f"{construct_id}_summary.csv", index=False)
    (output_dir / f"{construct_id}_report.md").write_text(
        markdown_report,
        encoding="utf-8",
    )


def control_comparisons_to_dataframe(
    comparisons: list[ControlComparisonResult],
) -> pd.DataFrame:
    """
    Convert list of ControlComparisonResult to a pandas DataFrame.
    """
    rows = [asdict(comparison) for comparison in comparisons]
    return pd.DataFrame(rows, columns=CONTROL_COMPARISON_COLUMNS)


def generate_screen_markdown_report(
    screen_result: ScreenResult,
    pose_dfs: dict[str, pd.DataFrame],
) -> str:
    """
    Generate a Markdown report for the full screen.
    """
    warnings = _collect_screen_warnings(pose_dfs)
    lines = [
        f"# LOXDockAudit Screen Report: {screen_result.screen_id}",
        "",
        "## Executive Decision",
        "",
        f"**{'PASS' if screen_result.overall_pass else 'FAIL'}**",
        "",
        screen_result.decision_text,
        "",
        "## Criteria Evaluation",
        "",
        _criteria_table(screen_result.criteria_results),
        "",
        "## Candidate Summary",
        "",
        _summary_table(screen_result.candidate_summary),
        "",
        "## Control Comparison",
        "",
        _control_comparison_table(screen_result.control_comparisons),
        "",
        "## Per-Construct Pose Details",
        "",
        _per_construct_pose_tables(pose_dfs),
        "",
        "## Warnings",
        "",
        _screen_warnings_table(warnings),
    ]
    return "\n".join(lines).rstrip() + "\n"


def save_screen_outputs(
    out_dir: str,
    screen_id: str,
    screen_result: ScreenResult,
    control_comparison_df: pd.DataFrame,
    all_summaries_df: pd.DataFrame,
    markdown_report: str,
) -> None:
    """
    Save screen-level CSV and Markdown outputs to out_dir.
    """
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    control_comparison_df.to_csv(
        output_dir / f"{screen_id}_control_comparison.csv",
        index=False,
    )
    all_summaries_df.to_csv(output_dir / f"{screen_id}_summaries.csv", index=False)
    (output_dir / f"{screen_id}_report.md").write_text(
        markdown_report,
        encoding="utf-8",
    )


def _executive_decision(summary: ConstructSummary) -> tuple[str, str]:
    if summary.strict_pass is True:
        return "PASS", "Construct meets strict pass criteria."
    if summary.strict_pass is False:
        return "FAIL", "Construct does not meet strict pass criteria."
    if summary.productive_count > 0:
        return "PASS", "At least one productive docking pose was found."
    return "FAIL", "No productive docking poses were found."


def _pose_table(pose_df: pd.DataFrame) -> str:
    columns = [
        "rank",
        "active_site_to_target_distance",
        "productive",
        "lox_contacts",
        "cbd_contacts",
    ]
    if pose_df.empty:
        return _markdown_table(columns, [])

    top_poses = pose_df.sort_values("rank", ascending=True).head(10)
    rows = [
        [
            row["rank"],
            _format_float(row["active_site_to_target_distance"]),
            row["productive"],
            row["lox_contacts"],
            row["cbd_contacts"],
        ]
        for _, row in top_poses.iterrows()
    ]
    return _markdown_table(columns, rows)


def _baseline_comparison_table(
    summary: ConstructSummary,
    baseline_summary: ConstructSummary,
) -> str:
    columns = [
        "construct_id",
        "productive_count",
        "total_poses",
        "best_distance",
        "best_productive_rank",
        "lox_contact_frequency",
        "cbd_contact_frequency",
    ]
    rows = []
    for item in (summary, baseline_summary):
        rows.append(
            [
                item.construct_id,
                item.productive_count,
                item.total_poses,
                _format_float(item.best_distance),
                _format_optional(item.best_productive_rank),
                _format_float(item.lox_contact_frequency),
                _format_float(item.cbd_contact_frequency),
            ]
        )
    return _markdown_table(columns, rows)


def _collect_warnings(pose_df: pd.DataFrame) -> list[tuple[int | str, str]]:
    if pose_df.empty or "warning" not in pose_df.columns:
        return []

    warnings: list[tuple[int | str, str]] = []
    for _, row in pose_df.iterrows():
        warning = row.get("warning")
        if warning is not None and str(warning):
            warnings.append((row.get("rank", ""), str(warning)))
    return warnings


def _warnings_table(warnings: list[tuple[int | str, str]]) -> str:
    if not warnings:
        return "No warnings."
    return _markdown_table(["rank", "warning"], warnings)


def _criteria_table(criteria_results: dict[str, Any]) -> str:
    rows = []
    skipped_notes = _criteria_skip_notes(criteria_results)
    for criterion, value in criteria_results.items():
        if criterion in {"criteria_checked", "criteria_skipped"}:
            continue
        result = _criterion_result(value)
        rows.append([criterion, result, skipped_notes.get(criterion, "")])
    return _markdown_table(["Criterion", "Result", "Notes"], rows)


def _criterion_result(value: object) -> str:
    if value is True:
        return "PASS"
    if value is False:
        return "FAIL"
    return "SKIPPED"


def _criteria_skip_notes(criteria_results: dict[str, Any]) -> dict[str, str]:
    notes: dict[str, str] = {}
    skipped = criteria_results.get("criteria_skipped", [])
    if not isinstance(skipped, list):
        return notes

    for item in skipped:
        text = str(item)
        criterion, separator, reason = text.partition(":")
        if separator:
            notes[criterion.strip()] = reason.strip()
    return notes


def _summary_table(summary: ConstructSummary) -> str:
    return _markdown_table(
        [
            "construct_id",
            "productive_count",
            "total_poses",
            "best_distance",
            "best_productive_rank",
            "lox_contact_frequency",
            "cbd_contact_frequency",
            "cbd_coupling",
            "strict_pass",
        ],
        [
            [
                summary.construct_id,
                summary.productive_count,
                summary.total_poses,
                _format_float(summary.best_distance),
                _format_optional(summary.best_productive_rank),
                _format_float(summary.lox_contact_frequency),
                _format_float(summary.cbd_contact_frequency),
                _format_float(summary.cbd_coupling),
                summary.strict_pass,
            ]
        ],
    )


def _control_comparison_table(comparisons: list[ControlComparisonResult]) -> str:
    frame = control_comparisons_to_dataframe(comparisons)
    if frame.empty:
        return _markdown_table(CONTROL_COMPARISON_COLUMNS, [])
    rows = [
        [row[column] for column in CONTROL_COMPARISON_COLUMNS]
        for _, row in frame.iterrows()
    ]
    return _markdown_table(CONTROL_COMPARISON_COLUMNS, rows)


def _per_construct_pose_tables(pose_dfs: dict[str, pd.DataFrame]) -> str:
    blocks: list[str] = []
    for construct_id, pose_df in pose_dfs.items():
        blocks.extend(
            [
                f"### {construct_id}",
                "",
                _screen_pose_table(pose_df),
                "",
            ]
        )
    return "\n".join(blocks).rstrip()


def _screen_pose_table(pose_df: pd.DataFrame) -> str:
    columns = [
        "rank",
        "active_site_to_target_distance",
        "productive",
        "lox_contacts",
        "cbd_contacts",
    ]
    if pose_df.empty:
        return _markdown_table(columns, [])

    top_poses = pose_df.sort_values("rank", ascending=True).head(5)
    rows = [
        [
            row["rank"],
            _format_float(row["active_site_to_target_distance"]),
            row["productive"],
            row["lox_contacts"],
            row["cbd_contacts"],
        ]
        for _, row in top_poses.iterrows()
    ]
    return _markdown_table(columns, rows)


def _collect_screen_warnings(
    pose_dfs: dict[str, pd.DataFrame],
) -> list[tuple[str, int | str, str]]:
    warnings: list[tuple[str, int | str, str]] = []
    for construct_id, pose_df in pose_dfs.items():
        for rank, warning in _collect_warnings(pose_df):
            warnings.append((construct_id, rank, warning))
    return warnings


def _screen_warnings_table(warnings: list[tuple[str, int | str, str]]) -> str:
    if not warnings:
        return "No warnings."
    return _markdown_table(["construct_id", "rank", "warning"], warnings)


def _markdown_table(columns: list[str], rows: Sequence[Sequence[object]]) -> str:
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    body = [
        "| " + " | ".join(_format_cell(value) for value in row) + " |"
        for row in rows
    ]
    return "\n".join([header, separator, *body])


def _format_cell(value: object) -> str:
    return str(value).replace("\n", " ").replace("|", "\\|")


def _format_float(value: Any) -> str:
    if value is None:
        return ""

    if isinstance(value, float) and math.isnan(value):
        return ""

    return f"{float(value):.3f}"


def _format_optional(value: object) -> str:
    if value is None:
        return "None"

    if isinstance(value, float) and math.isnan(value):
        return "None"

    return str(value)