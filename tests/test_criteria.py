from __future__ import annotations

from loxdockaudit.criteria import (
    build_control_comparisons,
    compute_overall_pass,
    evaluate_criteria,
    generate_decision_text,
)
from loxdockaudit.models import ConstructSummary, StrictCriteria


def test_evaluate_criteria_passes_against_baseline_and_controls() -> None:
    candidate = _summary("candidate", productive_count=3, best_distance=5.0, rank=1)
    baseline = _summary("baseline", productive_count=2, best_distance=6.0, rank=1)
    controls = [
        ("scrambled", _summary("scrambled", productive_count=1, best_distance=7.0, rank=1)),
        ("polyK", _summary("polyK", productive_count=1, best_distance=8.0, rank=1)),
    ]

    result = evaluate_criteria(candidate, baseline, controls, StrictCriteria())

    assert result["pass_count_vs_baseline"] is True
    assert result["best_distance_vs_baseline"] is True
    assert result["best_productive_rank"] is True
    assert result["beats_scrambled"] is True
    assert result["beats_polyK"] is True
    assert compute_overall_pass(result) is True


def test_evaluate_criteria_skips_missing_baseline_and_optional_cbd() -> None:
    candidate = _summary("candidate", productive_count=1, best_distance=5.0, rank=1)

    result = evaluate_criteria(candidate, None, [], StrictCriteria())

    assert result["pass_count_vs_baseline"] is None
    assert result["best_distance_vs_baseline"] is None
    assert result["cbd_contribution"] is None
    assert any("no baseline summary" in item for item in result["criteria_skipped"])


def test_evaluate_criteria_no_productive_pose_fails_rank() -> None:
    candidate = _summary("candidate", productive_count=0, best_distance=float("inf"), rank=None)

    result = evaluate_criteria(candidate, None, [], StrictCriteria())

    assert result["best_productive_rank"] is False
    assert "best_productive_rank" in result["criteria_checked"]


def test_cbd_contribution_checked_when_required() -> None:
    candidate = _summary("candidate", productive_count=1, best_distance=5.0, rank=1, cbd_frequency=0.0)

    result = evaluate_criteria(
        candidate,
        None,
        [],
        StrictCriteria(cbd_must_contribute=True),
    )

    assert result["cbd_contribution"] is False
    assert compute_overall_pass(result) is False


def test_build_control_comparisons_conclusion_rules() -> None:
    candidate = _summary("candidate", productive_count=2, best_distance=5.0, rank=1)
    controls = [
        ("baseline", _summary("baseline", productive_count=1, best_distance=6.0, rank=1)),
        ("scrambled", _summary("scrambled", productive_count=2, best_distance=5.2, rank=1)),
        ("polyK", _summary("polyK", productive_count=3, best_distance=4.0, rank=1)),
        ("noCBD", _summary("noCBD", productive_count=1, best_distance=8.0, rank=1)),
        ("inactive", _summary("inactive", productive_count=0, best_distance=9.0, rank=None)),
        ("benchmark", _summary("benchmark", productive_count=2, best_distance=5.0, rank=1)),
    ]

    comparisons = build_control_comparisons(candidate, controls)

    assert comparisons[0].conclusion == "Candidate improves on baseline."
    assert comparisons[1].conclusion == (
        "Candidate matches scrambled control. Binding may be nonspecific."
    )
    assert comparisons[2].conclusion == (
        "Candidate underperforms polyK control. Binding likely nonspecific."
    )
    assert comparisons[3].conclusion == "CBD domain contributes to productive geometry."
    assert comparisons[4].conclusion == "Inactive control included for geometry reference only."
    assert comparisons[5].conclusion == "Benchmark binder included for pipeline validation."


def test_generate_decision_text_pass() -> None:
    result = {
        "pass_count_vs_baseline": True,
        "criteria_checked": ["pass_count_vs_baseline"],
        "criteria_skipped": [],
    }

    text = generate_decision_text(
        "candidate",
        result,
        True,
        _summary("baseline", productive_count=1, best_distance=6.0, rank=1),
    )

    assert text.startswith("candidate passes strict criteria.")
    assert text.endswith("Candidate advances for further evaluation.")


def test_generate_decision_text_failure_mentions_scrambled() -> None:
    result = {
        "beats_scrambled": False,
        "criteria_checked": ["beats_scrambled"],
        "criteria_skipped": [],
    }

    text = generate_decision_text("candidate", result, False, None)

    assert text.startswith("candidate does not advance.")
    assert "nonspecific binding concern" in text
    assert text.endswith("Candidate does not advance.")


def test_best_productive_rank_above_max_fails() -> None:
    candidate = _summary("candidate", productive_count=1, best_distance=5.0, rank=4)

    result = evaluate_criteria(
        candidate,
        None,
        [],
        StrictCriteria(best_productive_rank_max=3),
    )

    assert result["best_productive_rank"] is False


def test_overall_pass_ignores_none_values() -> None:
    assert compute_overall_pass(
        {
            "pass_count_vs_baseline": True,
            "best_distance_vs_baseline": None,
            "criteria_checked": ["pass_count_vs_baseline"],
            "criteria_skipped": ["best_distance_vs_baseline: no baseline summary"],
        }
    ) is True


def test_scrambled_higher_count_conclusion_mentions_nonspecific() -> None:
    candidate = _summary("candidate", productive_count=1, best_distance=4.0, rank=1)
    scrambled = _summary("scrambled", productive_count=2, best_distance=8.0, rank=1)

    comparison = build_control_comparisons(candidate, [("scrambled", scrambled)])[0]

    assert "nonspecific" in comparison.conclusion


def test_decision_text_failed_rank_mentions_rank() -> None:
    result = {
        "best_productive_rank": False,
        "criteria_checked": ["best_productive_rank"],
        "criteria_skipped": [],
    }

    text = generate_decision_text("candidate", result, False, None)

    assert "rank" in text


def test_orientation_criteria_skip_when_orientation_disabled() -> None:
    result = evaluate_criteria(
        _summary("candidate", 2, 5.0, 1),
        _summary("baseline", 1, 6.0, 1),
        [],
        StrictCriteria(
            must_beat_scrambled=False,
            must_beat_polyK=False,
            fully_productive_count_greater_than_baseline=True,
            best_orientation_angle_max_deg=75.0,
            best_fully_productive_rank_max=5,
        ),
    )

    assert result["fully_productive_count_greater_than_baseline"] is None
    assert result["best_orientation_angle_max_deg"] is None
    assert result["best_fully_productive_rank_max"] is None
    assert any("orientation disabled" in item for item in result["criteria_skipped"])


def test_orientation_criteria_pass_when_enabled_and_metrics_pass() -> None:
    candidate = _summary(
        "candidate",
        2,
        5.0,
        1,
        orientation_enabled=True,
        fully_productive_count=2,
        best_orientation_angle_deg=60.0,
        best_fully_productive_rank=2,
    )
    baseline = _summary(
        "baseline",
        1,
        6.0,
        1,
        orientation_enabled=True,
        fully_productive_count=1,
        best_orientation_angle_deg=80.0,
        best_fully_productive_rank=3,
    )

    result = evaluate_criteria(
        candidate,
        baseline,
        [],
        StrictCriteria(
            must_beat_scrambled=False,
            must_beat_polyK=False,
            fully_productive_count_greater_than_baseline=True,
            best_orientation_angle_max_deg=75.0,
            best_fully_productive_rank_max=5,
        ),
    )

    assert result["fully_productive_count_greater_than_baseline"] is True
    assert result["best_orientation_angle_max_deg"] is True
    assert result["best_fully_productive_rank_max"] is True


def test_orientation_criteria_fail_when_enabled_and_metrics_fail() -> None:
    candidate = _summary(
        "candidate",
        2,
        5.0,
        1,
        orientation_enabled=True,
        fully_productive_count=1,
        best_orientation_angle_deg=90.0,
        best_fully_productive_rank=None,
    )
    baseline = _summary(
        "baseline",
        1,
        6.0,
        1,
        orientation_enabled=True,
        fully_productive_count=1,
        best_orientation_angle_deg=80.0,
        best_fully_productive_rank=3,
    )

    result = evaluate_criteria(
        candidate,
        baseline,
        [],
        StrictCriteria(
            must_beat_scrambled=False,
            must_beat_polyK=False,
            fully_productive_count_greater_than_baseline=True,
            best_orientation_angle_max_deg=75.0,
            best_fully_productive_rank_max=5,
        ),
    )

    assert result["fully_productive_count_greater_than_baseline"] is False
    assert result["best_orientation_angle_max_deg"] is False
    assert result["best_fully_productive_rank_max"] is False


def _summary(
    construct_id: str,
    productive_count: int,
    best_distance: float,
    rank: int | None,
    cbd_frequency: float = 0.0,
    orientation_enabled: bool = False,
    fully_productive_count: int = 0,
    best_orientation_angle_deg: float | None = None,
    best_fully_productive_rank: int | None = None,
) -> ConstructSummary:
    return ConstructSummary(
        construct_id=construct_id,
        productive_count=productive_count,
        total_poses=3,
        best_distance=best_distance,
        best_productive_rank=rank,
        lox_contact_frequency=1.0,
        cbd_contact_frequency=cbd_frequency,
        cbd_coupling=0.0,
        strict_pass=None,
        orientation_enabled=orientation_enabled,
        fully_productive_count=fully_productive_count,
        best_orientation_angle_deg=best_orientation_angle_deg,
        best_fully_productive_rank=best_fully_productive_rank,
    )
