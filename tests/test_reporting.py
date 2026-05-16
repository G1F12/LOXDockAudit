from __future__ import annotations

from pathlib import Path

from loxdockaudit.models import (
    ActiveSiteConfig,
    ConstructConfig,
    ConstructSummary,
    ControlComparisonResult,
    DomainRange,
    PoseMetrics,
    ScreenResult,
    StructuralQCResult,
    TargetResidue,
)
from loxdockaudit.reporting import (
    control_comparisons_to_dataframe,
    generate_markdown_report,
    generate_screen_markdown_report,
    poses_to_dataframe,
    save_outputs,
    save_screen_outputs,
    summaries_to_dataframe,
)


def test_poses_to_dataframe_empty() -> None:
    df = poses_to_dataframe([])
    assert df.empty
    assert df.columns.tolist() == [
        "rank",
        "model_path",
        "active_site_to_target_distance",
        "productive",
        "lox_contacts",
        "cbd_contacts",
        "closest_active_site_resi",
        "warning",
    ]


def test_poses_to_dataframe_three_sorted() -> None:
    df = poses_to_dataframe(
        [
            PoseMetrics(3, "p3", 7.0, True, 1, 0),
            PoseMetrics(1, "p1", 5.0, True, 1, 0),
            PoseMetrics(2, "p2", 15.0, False, 0, 0),
        ]
    )
    assert df.shape == (3, 8)
    assert df["rank"].tolist() == [1, 2, 3]


def test_summaries_to_dataframe_productive_fraction() -> None:
    df = summaries_to_dataframe([_summary(productive_count=2, total_poses=3)])
    assert df.loc[0, "productive_fraction"] == 2 / 3


def test_summaries_to_dataframe_zero_total() -> None:
    df = summaries_to_dataframe([_summary(productive_count=0, total_poses=0)])
    assert df.loc[0, "productive_fraction"] == 0.0


def test_generate_markdown_report_headers() -> None:
    report = generate_markdown_report(_config(), poses_to_dataframe([]), _summary(), "audit")
    for header in [
        "## Executive Decision",
        "## Input Audit",
        "## Structural QC",
        "## Docking Pose Summary",
        "## Productive Geometry",
        "## Warnings",
    ]:
        assert header in report


def test_generate_markdown_report_pass_or_fail() -> None:
    report = generate_markdown_report(_config(), poses_to_dataframe([]), _summary(), "audit")
    assert "PASS" in report or "FAIL" in report


def test_generate_markdown_report_pipe_table() -> None:
    report = generate_markdown_report(_config(), poses_to_dataframe([]), _summary(), "audit")
    assert "|" in report


def test_generate_markdown_report_structural_qc_not_configured() -> None:
    report = generate_markdown_report(_config(), poses_to_dataframe([]), _summary(), "audit")

    assert "Structural QC not configured." in report


def test_generate_markdown_report_includes_structural_qc_result() -> None:
    report = generate_markdown_report(
        _config(),
        poses_to_dataframe([]),
        _summary(),
        "audit",
        qc_result=StructuralQCResult(
            construct_id="construct",
            his_triad_pass=True,
            lys_tyr_pass=True,
            disulfide_pass=True,
            active_site_accessible=True,
            plddt_pass=None,
            pae_pass=None,
            fold_qc_pass=True,
            fold_corrupted=False,
            warnings=[],
            details={"disulfide": {"skipped": True}},
        ),
    )

    assert "## Structural QC: construct" in report
    assert "| His triad geometry | PASS |" in report


def test_save_outputs_creates_three_files(tmp_path: Path) -> None:
    out = tmp_path / "nested" / "out"
    save_outputs(out, "construct", poses_to_dataframe([]), summaries_to_dataframe([_summary()]), "report")
    assert sorted(path.name for path in out.iterdir()) == [
        "construct_poses.csv",
        "construct_report.md",
        "construct_summary.csv",
    ]


def test_save_outputs_creates_out_dir(tmp_path: Path) -> None:
    out = tmp_path / "missing" / "out"
    save_outputs(out, "construct", poses_to_dataframe([]), summaries_to_dataframe([_summary()]), "report")
    assert out.is_dir()


def test_control_comparisons_to_dataframe_columns() -> None:
    df = control_comparisons_to_dataframe([_comparison()])
    assert df.columns.tolist() == [
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
    assert df.loc[0, "conclusion"] == "Candidate improves on baseline."


def test_control_comparisons_to_dataframe_empty_has_columns() -> None:
    df = control_comparisons_to_dataframe([])
    assert df.empty
    assert df.columns.tolist() == [
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
    assert df.columns[-1] == "conclusion"


def test_generate_screen_markdown_report_sections_and_tables() -> None:
    pose_df = poses_to_dataframe(
        [
            PoseMetrics(1, "pose1", 5.0, True, 1, 0),
            PoseMetrics(2, "pose2", 9.0, False, 0, 0, warning="target far"),
        ]
    )
    report = generate_screen_markdown_report(
        _screen_result(),
        {"candidate": pose_df},
    )

    expected_sections = [
        "# LOXDockAudit Screen Report: screen",
        "## Executive Decision",
        "## Criteria Evaluation",
        "## Candidate Summary",
        "## Control Comparison",
        "## Per-Construct Pose Details",
        "## Warnings",
        "### candidate",
    ]
    for section in expected_sections:
        assert section in report
    assert "| Criterion | Result | Notes |" in report
    assert "| construct_id | control_type |" in report
    assert "target far" in report


def test_generate_screen_markdown_report_requested_minimal_case() -> None:
    report = generate_screen_markdown_report(_screen_result(), {})

    assert "FAIL" in report[:200]
    assert "## Criteria Evaluation" in report
    assert "SKIPPED" in report


def test_generate_screen_markdown_report_handles_empty_pose_dfs() -> None:
    report = generate_screen_markdown_report(_screen_result(), {})
    assert "## Per-Construct Pose Details" in report
    assert "No warnings." in report


def test_generate_screen_markdown_report_limits_pose_tables_to_top_five() -> None:
    pose_df = poses_to_dataframe(
        [
            PoseMetrics(rank, f"pose{rank}", float(rank), True, 1, 0)
            for rank in range(1, 7)
        ]
    )
    report = generate_screen_markdown_report(_screen_result(), {"candidate": pose_df})

    assert "| 5 | 5.000 | True | 1 | 0 |" in report
    assert "| 6 | 6.000 | True | 1 | 0 |" not in report


def test_generate_screen_markdown_report_no_warnings() -> None:
    report = generate_screen_markdown_report(
        _screen_result(),
        {"candidate": poses_to_dataframe([PoseMetrics(1, "pose1", 5.0, True, 1, 0)])},
    )
    assert "No warnings." in report


def test_save_screen_outputs_creates_three_files(tmp_path: Path) -> None:
    out = tmp_path / "screen_out"
    save_screen_outputs(
        str(out),
        "screen",
        _screen_result(),
        control_comparisons_to_dataframe([_comparison()]),
        summaries_to_dataframe([_summary()]),
        "report",
    )
    assert sorted(path.name for path in out.iterdir()) == [
        "screen_control_comparison.csv",
        "screen_report.md",
        "screen_summaries.csv",
    ]


def test_save_screen_outputs_overwrites_existing_files(tmp_path: Path) -> None:
    out = tmp_path / "screen_out"
    save_screen_outputs(
        str(out),
        "screen",
        _screen_result(),
        control_comparisons_to_dataframe([_comparison()]),
        summaries_to_dataframe([_summary()]),
        "first",
    )
    save_screen_outputs(
        str(out),
        "screen",
        _screen_result(),
        control_comparisons_to_dataframe([_comparison()]),
        summaries_to_dataframe([_summary()]),
        "second",
    )

    assert (out / "screen_report.md").read_text(encoding="utf-8") == "second"


def _config() -> ConstructConfig:
    return ConstructConfig(
        construct_id="construct",
        construct_type="lox_baseline",
        models_dir="models",
        receptor_chains=["A"],
        ligand_chains=["D"],
        active_site=ActiveSiteConfig([{"label": "His292", "resi": 124}]),
        target_residues=[TargetResidue("A", 11, "LYS")],
        lox_domain=DomainRange(1, 200),
    )


def _summary(productive_count: int = 1, total_poses: int = 2) -> ConstructSummary:
    return ConstructSummary("construct", productive_count, total_poses, 5.0, 1, 0.5, 0.0, 0.0, None)


def _comparison() -> ControlComparisonResult:
    return ControlComparisonResult(
        construct_id="baseline",
        control_type="baseline",
        candidate_productive_count=2,
        control_productive_count=1,
        candidate_best_distance=5.0,
        control_best_distance=6.0,
        candidate_beats_control_count=True,
        candidate_beats_control_distance=True,
        conclusion="Candidate improves on baseline.",
    )


def _screen_result() -> ScreenResult:
    return ScreenResult(
        screen_id="screen",
        candidate_summary=_summary(productive_count=2, total_poses=3),
        baseline_summary=_summary(productive_count=1, total_poses=3),
        control_comparisons=[_comparison()],
        criteria_results={
            "pass_count_vs_baseline": True,
            "best_distance_vs_baseline": True,
            "best_productive_rank": True,
            "cbd_contribution": None,
            "beats_scrambled": False,
            "beats_polyK": None,
            "criteria_checked": [
                "pass_count_vs_baseline",
                "best_distance_vs_baseline",
                "best_productive_rank",
                "beats_scrambled",
            ],
            "criteria_skipped": [
                "cbd_contribution: criterion disabled",
                "beats_polyK: no polyK control",
            ],
        },
        overall_pass=False,
        decision_text="candidate does not advance. Nonspecific binding concern. Candidate does not advance.",
    )
