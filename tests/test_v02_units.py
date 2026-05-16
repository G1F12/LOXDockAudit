from __future__ import annotations

from pathlib import Path

import pytest

from loxdockaudit.criteria import (
    build_control_comparisons,
    compute_overall_pass,
    evaluate_criteria,
    generate_decision_text,
)
from loxdockaudit.models import (
    ConstructSummary,
    ControlComparisonResult,
    ScreenResult,
    StrictCriteria,
)
from loxdockaudit.preprocessing import clean_pdb_lines
from loxdockaudit.reporting import (
    control_comparisons_to_dataframe,
    generate_screen_markdown_report,
    save_screen_outputs,
    summaries_to_dataframe,
)
from loxdockaudit.screen import pearson_correlation
from loxdockaudit.sorting import natural_sort_key


def test_natural_sort_order() -> None:
    filenames = [
        "model_1.pdb",
        "model_10.pdb",
        "model_2.pdb",
        "model_3.pdb",
        "model_9.pdb",
    ]
    sorted_files = sorted(filenames, key=natural_sort_key)
    assert sorted_files == [
        "model_1.pdb",
        "model_2.pdb",
        "model_3.pdb",
        "model_9.pdb",
        "model_10.pdb",
    ]


def test_hdock_preprocessing_strips_wrapper_records() -> None:
    lines = [
        "REMARK produced by HDOCK\n",
        "MODEL        1\n",
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n",
        "ENDMDL\n",
    ]
    cleaned = clean_pdb_lines(lines)
    assert not any(line.startswith(("REMARK", "MODEL", "ENDMDL")) for line in cleaned)
    assert cleaned[-1] == "END\n"


def test_hdock_preprocessing_preserves_ter_and_drops_malformed_atoms() -> None:
    lines = [
        "ATOM      1  CA  ALA A   1       bad     0.000   0.000  1.00 20.00           C\n",
        "ATOM      2  CA  GLY A   2       1.000   0.000   0.000  1.00 20.00           C\n",
        "TER\n",
    ]
    cleaned = clean_pdb_lines(lines)
    assert len([line for line in cleaned if line.startswith("ATOM")]) == 1
    assert "TER\n" in cleaned


def test_hdock_preprocessing_normalizes_blank_chain_id() -> None:
    line = "ATOM      1  CA  ALA     1       0.000   0.000   0.000  1.00 20.00           C\n"
    cleaned = clean_pdb_lines([line])
    assert cleaned[0][21] == "X"


def test_hdock_preprocessing_maps_ligand_header_chain_to_d() -> None:
    lines = [
        "HEADER lig_1.pdb\n",
        "ATOM      1  CA  ALA X   1       0.000   0.000   0.000\n",
    ]
    cleaned = clean_pdb_lines(lines)
    assert cleaned[0][21] == "D"


def test_v02_01_no_baseline_skips_baseline_count() -> None:
    result = evaluate_criteria(_summary("candidate", 1, 5.0, 1), None, [], StrictCriteria())
    assert result["pass_count_vs_baseline"] is None


def test_v02_02_candidate_beats_baseline_on_count() -> None:
    result = evaluate_criteria(
        _summary("candidate", 2, 5.0, 1),
        _summary("baseline", 1, 6.0, 1),
        [],
        StrictCriteria(must_beat_scrambled=False, must_beat_polyK=False),
    )
    assert result["pass_count_vs_baseline"] is True


def test_v02_03_best_productive_rank_above_max_fails() -> None:
    result = evaluate_criteria(_summary("candidate", 1, 5.0, 4), None, [], StrictCriteria())
    assert result["best_productive_rank"] is False


def test_v02_04_missing_best_productive_rank_fails() -> None:
    result = evaluate_criteria(_summary("candidate", 0, float("inf"), None), None, [], StrictCriteria())
    assert result["best_productive_rank"] is False


def test_v02_05_missing_scrambled_is_skipped() -> None:
    result = evaluate_criteria(_summary("candidate", 1, 5.0, 1), None, [], StrictCriteria())
    assert result["beats_scrambled"] is None


def test_v02_06_overall_pass_true_true_none() -> None:
    assert compute_overall_pass({"a": True, "b": True, "c": None}) is True


def test_v02_07_overall_pass_true_false_none() -> None:
    assert compute_overall_pass({"a": True, "b": False, "c": None}) is False


def test_v02_08_worse_than_scrambled_mentions_nonspecific() -> None:
    comparison = build_control_comparisons(
        _summary("candidate", 1, 5.0, 1),
        [("scrambled", _summary("scrambled", 2, 4.0, 1))],
    )[0]
    assert "nonspecific" in comparison.conclusion


def test_v02_09_decision_text_failed_rank_mentions_rank() -> None:
    text = generate_decision_text(
        "candidate",
        {"best_productive_rank": False, "criteria_checked": ["best_productive_rank"]},
        False,
        None,
    )
    assert "rank" in text


def test_v02_10_control_comparisons_empty_dataframe() -> None:
    df = control_comparisons_to_dataframe([])
    assert df.empty
    assert df.shape[1] == 9


def test_v02_11_screen_report_has_all_sections() -> None:
    report = generate_screen_markdown_report(_screen_result(), {})
    for section in [
        "## Executive Decision",
        "## Criteria Evaluation",
        "## Candidate Summary",
        "## Control Comparison",
        "## Per-Construct Pose Details",
        "## Warnings",
    ]:
        assert section in report


def test_v02_12_screen_report_fail_in_first_200_chars() -> None:
    report = generate_screen_markdown_report(_screen_result(), {})
    assert "FAIL" in report[:200]


def test_v02_13_save_screen_outputs_creates_three_files(tmp_path: Path) -> None:
    save_screen_outputs(
        str(tmp_path),
        "screen",
        _screen_result(),
        control_comparisons_to_dataframe([_comparison()]),
        summaries_to_dataframe([_summary("candidate", 1, 3.0, 1)]),
        "report",
    )
    assert sorted(path.name for path in tmp_path.iterdir()) == [
        "screen_control_comparison.csv",
        "screen_report.md",
        "screen_summaries.csv",
    ]


def test_v02_14_pearson_positive() -> None:
    assert pearson_correlation([1, 2, 3], [1, 2, 3]) == pytest.approx(1.0)


def test_v02_15_pearson_constant_returns_zero() -> None:
    assert pearson_correlation([1, 1, 1], [1, 2, 3]) == 0.0


def test_v02_16_pearson_short_returns_zero() -> None:
    assert pearson_correlation([1, 2], [1, 2]) == 0.0


def _summary(
    construct_id: str,
    productive_count: int,
    best_distance: float,
    rank: int | None,
) -> ConstructSummary:
    return ConstructSummary(
        construct_id=construct_id,
        productive_count=productive_count,
        total_poses=3,
        best_distance=best_distance,
        best_productive_rank=rank,
        lox_contact_frequency=1.0,
        cbd_contact_frequency=0.0,
        cbd_coupling=0.0,
        strict_pass=None,
    )


def _comparison() -> ControlComparisonResult:
    return ControlComparisonResult(
        construct_id="baseline",
        control_type="baseline",
        candidate_productive_count=2,
        control_productive_count=1,
        candidate_best_distance=6.0,
        control_best_distance=9.0,
        candidate_beats_control_count=True,
        candidate_beats_control_distance=True,
        conclusion="Candidate improves on baseline.",
    )


def _screen_result() -> ScreenResult:
    return ScreenResult(
        screen_id="screen",
        candidate_summary=_summary("candidate", 2, 6.0, 1),
        baseline_summary=_summary("baseline", 1, 9.0, 1),
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
            "criteria_skipped": ["cbd_contribution: criterion disabled"],
        },
        overall_pass=False,
        decision_text="candidate does not advance. The candidate failed the scrambled-control criterion, indicating nonspecific binding concern. Candidate does not advance.",
    )
