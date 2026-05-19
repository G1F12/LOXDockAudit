"""Strict criteria evaluation for LOXDockAudit screens."""

from __future__ import annotations

from typing import Any

from loxdockaudit.models import (
    ConstructSummary,
    ControlComparisonResult,
    StrictCriteria,
)


def evaluate_criteria(
    candidate: ConstructSummary,
    baseline: ConstructSummary | None,
    controls: list[tuple[str, ConstructSummary]],
    criteria: StrictCriteria,
) -> dict[str, Any]:
    """
    Evaluate candidate strict criteria against baseline and controls.
    """
    result: dict[str, Any] = {
        "pass_count_vs_baseline": None,
        "best_distance_vs_baseline": None,
        "best_productive_rank": None,
        "cbd_contribution": None,
        "beats_scrambled": None,
        "beats_polyK": None,
        "criteria_checked": [],
        "criteria_skipped": [],
    }

    if criteria.pass_count_greater_than_baseline:
        if baseline is None:
            _skip(result, "pass_count_vs_baseline", "no baseline summary")
        else:
            result["pass_count_vs_baseline"] = (
                candidate.productive_count > baseline.productive_count
            )
            _checked(result, "pass_count_vs_baseline")
    else:
        _skip(result, "pass_count_vs_baseline", "criterion disabled")

    if criteria.best_distance_lte_baseline:
        if baseline is None:
            _skip(result, "best_distance_vs_baseline", "no baseline summary")
        else:
            result["best_distance_vs_baseline"] = (
                candidate.best_distance <= baseline.best_distance
            )
            _checked(result, "best_distance_vs_baseline")
    else:
        _skip(result, "best_distance_vs_baseline", "criterion disabled")

    result["best_productive_rank"] = (
        candidate.best_productive_rank is not None
        and candidate.best_productive_rank <= criteria.best_productive_rank_max
    )
    _checked(result, "best_productive_rank")

    if criteria.cbd_must_contribute:
        result["cbd_contribution"] = candidate.cbd_contact_frequency > 0
        _checked(result, "cbd_contribution")
    else:
        _skip(result, "cbd_contribution", "criterion disabled")

    scrambled = _first_control("scrambled", controls)
    if criteria.must_beat_scrambled:
        if scrambled is None:
            _skip(result, "beats_scrambled", "no scrambled control")
        else:
            result["beats_scrambled"] = _is_better(candidate, scrambled)
            _checked(result, "beats_scrambled")
    else:
        _skip(result, "beats_scrambled", "criterion disabled")

    poly_k = _first_control("polyK", controls)
    if criteria.must_beat_polyK:
        if poly_k is None:
            _skip(result, "beats_polyK", "no polyK control")
        else:
            result["beats_polyK"] = _is_better(candidate, poly_k)
            _checked(result, "beats_polyK")
    else:
        _skip(result, "beats_polyK", "criterion disabled")

    # v0.4: optional orientation-aware strict criteria.
    _evaluate_orientation_criteria(
        result,
        candidate,
        baseline,
        controls,
        criteria,
    )

    return result


def compute_overall_pass(criteria_result: dict[str, Any]) -> bool:
    """
    Return True only if all evaluated criteria are True.
    """
    criterion_values = [
        value
        for key, value in criteria_result.items()
        if key not in {"criteria_checked", "criteria_skipped"} and value is not None
    ]
    return all(criterion_values)


def build_control_comparisons(
    candidate: ConstructSummary,
    controls: list[tuple[str, ConstructSummary]],
) -> list[ControlComparisonResult]:
    """
    Build control comparison results for each control summary.
    """
    return [
        ControlComparisonResult(
            construct_id=control.construct_id,
            control_type=control_type,
            candidate_productive_count=candidate.productive_count,
            control_productive_count=control.productive_count,
            candidate_best_distance=candidate.best_distance,
            control_best_distance=control.best_distance,
            candidate_beats_control_count=(
                candidate.productive_count > control.productive_count
            ),
            candidate_beats_control_distance=(
                candidate.best_distance < control.best_distance
            ),
            conclusion=_comparison_conclusion(candidate, control_type, control),
        )
        for control_type, control in controls
    ]


def generate_decision_text(
    candidate_id: str,
    criteria_result: dict[str, Any],
    overall_pass: bool,
    baseline_summary: ConstructSummary | None,
) -> str:
    """
    Generate a factual human-readable screen decision paragraph.
    """
    if overall_pass:
        baseline_text = _baseline_text(baseline_summary)
        return (
            f"{candidate_id} passes strict criteria. "
            f"{baseline_text} "
            "Candidate advances for further evaluation."
        )

    failing_criterion = _first_failing_criterion(criteria_result)
    failure_text = _failure_text(failing_criterion, criteria_result)
    return (
        f"{candidate_id} does not advance. "
        f"{failure_text} "
        "Candidate does not advance."
    )


def _checked(result: dict[str, Any], criterion_name: str) -> None:
    result["criteria_checked"].append(criterion_name)


def _skip(result: dict[str, Any], criterion_name: str, reason: str) -> None:
    result["criteria_skipped"].append(f"{criterion_name}: {reason}")


# v0.4: evaluate only explicitly configured orientation criteria.
def _evaluate_orientation_criteria(
    result: dict[str, Any],
    candidate: ConstructSummary,
    baseline: ConstructSummary | None,
    controls: list[tuple[str, ConstructSummary]],
    criteria: StrictCriteria,
) -> None:
    requested = {
        "fully_productive_count_greater_than_baseline": (
            criteria.fully_productive_count_greater_than_baseline is not None
        ),
        "best_orientation_angle_max_deg": (
            criteria.best_orientation_angle_max_deg is not None
        ),
        "best_fully_productive_rank_max": (
            criteria.best_fully_productive_rank_max is not None
        ),
    }
    if not any(requested.values()):
        return

    for criterion_name, enabled in requested.items():
        if enabled:
            result[criterion_name] = None

    if _any_orientation_disabled(candidate, baseline, controls):
        for criterion_name, enabled in requested.items():
            if enabled:
                _skip(result, criterion_name, "orientation disabled")
        return

    if requested["fully_productive_count_greater_than_baseline"]:
        if not criteria.fully_productive_count_greater_than_baseline:
            _skip(
                result,
                "fully_productive_count_greater_than_baseline",
                "criterion disabled",
            )
        elif baseline is None:
            _skip(
                result,
                "fully_productive_count_greater_than_baseline",
                "no baseline summary",
            )
        else:
            result["fully_productive_count_greater_than_baseline"] = (
                candidate.fully_productive_count > baseline.fully_productive_count
            )
            _checked(result, "fully_productive_count_greater_than_baseline")

    if requested["best_orientation_angle_max_deg"]:
        angle_threshold = criteria.best_orientation_angle_max_deg
        assert angle_threshold is not None
        if candidate.best_orientation_angle_deg is None:
            result["best_orientation_angle_max_deg"] = False
        else:
            result["best_orientation_angle_max_deg"] = (
                candidate.best_orientation_angle_deg <= angle_threshold
            )
        _checked(result, "best_orientation_angle_max_deg")

    if requested["best_fully_productive_rank_max"]:
        rank_threshold = criteria.best_fully_productive_rank_max
        assert rank_threshold is not None
        result["best_fully_productive_rank_max"] = (
            candidate.best_fully_productive_rank is not None
            and candidate.best_fully_productive_rank <= rank_threshold
        )
        _checked(result, "best_fully_productive_rank_max")


# v0.4: old summaries default orientation_enabled=False, causing requested
# orientation criteria to skip rather than fail.
def _any_orientation_disabled(
    candidate: ConstructSummary,
    baseline: ConstructSummary | None,
    controls: list[tuple[str, ConstructSummary]],
) -> bool:
    summaries = [candidate]
    if baseline is not None:
        summaries.append(baseline)
    summaries.extend(summary for _, summary in controls)
    return any(not summary.orientation_enabled for summary in summaries)


def _first_control(
    control_type: str,
    controls: list[tuple[str, ConstructSummary]],
) -> ConstructSummary | None:
    for current_type, summary in controls:
        if current_type == control_type:
            return summary
    return None


def _is_better(candidate: ConstructSummary, control: ConstructSummary) -> bool:
    if candidate.productive_count != control.productive_count:
        return candidate.productive_count > control.productive_count
    return candidate.best_distance < control.best_distance


def _is_equal(candidate: ConstructSummary, control: ConstructSummary) -> bool:
    return (
        candidate.productive_count == control.productive_count
        and abs(candidate.best_distance - control.best_distance) <= 0.5
    )


def _comparison_conclusion(
    candidate: ConstructSummary,
    control_type: str,
    control: ConstructSummary,
) -> str:
    better = _is_better(candidate, control)
    equal = _is_equal(candidate, control)

    if control_type == "baseline":
        if equal:
            return "Candidate matches baseline."
        if better:
            return "Candidate improves on baseline."
        return "Candidate does not improve on baseline."

    if control_type in {"scrambled", "polyK"}:
        if equal:
            return (
                f"Candidate matches {control_type} control. "
                "Binding may be nonspecific."
            )
        if better:
            return f"Candidate outperforms {control_type} control."
        return (
            f"Candidate underperforms {control_type} control. "
            "Binding likely nonspecific."
        )

    if control_type == "noCBD":
        if better:
            return "CBD domain contributes to productive geometry."
        return "CBD domain does not improve on noCBD baseline."

    if control_type == "inactive":
        return "Inactive control included for geometry reference only."

    if control_type == "benchmark":
        return "Benchmark binder included for pipeline validation."

    return f"{control_type} control included for comparison."


def _baseline_text(baseline_summary: ConstructSummary | None) -> str:
    if baseline_summary is None:
        return "No baseline summary was available for direct comparison."
    return (
        "Candidate comparison used baseline productive count "
        f"{baseline_summary.productive_count} and best distance "
        f"{baseline_summary.best_distance:.2f} A."
    )


def _first_failing_criterion(criteria_result: dict[str, Any]) -> str:
    for criterion in criteria_result.get("criteria_checked", []):
        if criteria_result.get(criterion) is False:
            return criterion
    return "no evaluated criterion passed"


def _failure_text(criterion: str, criteria_result: dict[str, Any]) -> str:
    if criterion == "best_productive_rank":
        return "The first productive pose rank exceeded the maximum allowed rank."
    if criterion == "beats_scrambled":
        return "The candidate failed the scrambled-control criterion, indicating nonspecific binding concern."
    if criterion == "pass_count_vs_baseline":
        return "The candidate did not exceed the baseline productive pose count."
    if criterion == "best_distance_vs_baseline":
        return "The candidate best distance was worse than the baseline best distance."
    if criterion == "cbd_contribution":
        return "The candidate did not show CBD contact contribution."
    if criterion == "beats_polyK":
        return "The candidate failed the polyK-control criterion."
    if criterion == "fully_productive_count_greater_than_baseline":
        return "The candidate did not exceed the baseline fully productive pose count."
    if criterion == "best_orientation_angle_max_deg":
        return "The candidate best orientation angle exceeded the maximum allowed angle."
    if criterion == "best_fully_productive_rank_max":
        return "The first fully productive pose rank exceeded the maximum allowed rank."
    return f"The candidate failed criterion {criterion}."
