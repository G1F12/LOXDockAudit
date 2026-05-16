"""Screen-level orchestration for LOXDockAudit."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import numpy as np

from loxdockaudit.config import load_construct_config
from loxdockaudit.contacts import analyze_pose_contacts, compute_contact_frequencies
from loxdockaudit.criteria import (
    build_control_comparisons,
    compute_overall_pass,
    evaluate_criteria,
    generate_decision_text,
)
from loxdockaudit.distances import active_site_to_target_distance
from loxdockaudit.models import (
    ConstructConfig,
    ConstructSummary,
    PoseMetrics,
    ScreenConfig,
    ScreenResult,
)
from loxdockaudit.pdb_parser import load_structure
from loxdockaudit.sorting import natural_sort_key


def run_single_construct(
    construct_config: ConstructConfig,
    top_n: int = 10,
    verbose: bool = False,
) -> tuple[ConstructSummary, list[PoseMetrics]]:
    """
    Run the full v0.1 analysis pipeline for one construct.
    """
    pdb_paths = sorted(
        Path(construct_config.models_dir).glob("*.pdb"),
        key=lambda path: natural_sort_key(str(path)),
    )[:top_n]
    pose_metrics: list[PoseMetrics] = []
    contact_results: list[dict[str, Any]] = []

    target_residue = construct_config.target_residues[0]
    active_site_resi = [
        int(residue["resi"])
        for residue in construct_config.active_site.residues
        if "resi" in residue
    ]

    for rank, pdb_path in enumerate(pdb_paths, start=1):
        try:
            structure = load_structure(str(pdb_path), structure_id=f"pose_{rank}")
            distance_result = active_site_to_target_distance(
                structure=structure,
                active_site_resi=active_site_resi,
                active_site_chains=construct_config.ligand_chains,
                target_resi=target_residue.resi,
                target_chain=target_residue.chain,
                threshold=construct_config.productive_distance_threshold,
                target_atom_names=[target_residue.atom],
            )
            contact_result = analyze_pose_contacts(
                structure=structure,
                receptor_chains=construct_config.receptor_chains,
                ligand_chains=construct_config.ligand_chains,
                lox_domain_resi_range=(
                    construct_config.lox_domain.start,
                    construct_config.lox_domain.end,
                ),
                cbd_domain_resi_range=_domain_tuple(construct_config.cbd_domain),
                contact_cutoff=construct_config.contact_distance_threshold,
            )
        except Exception as exc:
            message = f"{construct_config.construct_id}: skipped {pdb_path}: {exc}"
            if verbose:
                warnings.warn(message, RuntimeWarning, stacklevel=2)
            else:
                warnings.warn(message, RuntimeWarning, stacklevel=2)
            continue

        contact_results.append(contact_result)
        pose_metrics.append(
            PoseMetrics(
                rank=rank,
                model_path=str(pdb_path),
                active_site_to_target_distance=distance_result["distance_angstrom"],
                productive=distance_result["productive"],
                lox_contacts=contact_result["lox_contacts"],
                cbd_contacts=contact_result["cbd_contacts"],
                closest_active_site_resi=distance_result["closest_active_site_resi"],
                warning=distance_result["warning"],
            )
        )

    if not pose_metrics:
        raise RuntimeError(f"{construct_config.construct_id}: no poses processed")

    summary = _build_summary(construct_config, pose_metrics, contact_results)
    return summary, pose_metrics


def pearson_correlation(x: list[float], y: list[float]) -> float:
    """
    Compute Pearson correlation coefficient for two numeric vectors.
    """
    if len(x) < 3 or len(y) < 3 or len(x) != len(y):
        return 0.0

    values_x = np.array(x, dtype=float)
    values_y = np.array(y, dtype=float)
    if float(np.std(values_x)) == 0.0 or float(np.std(values_y)) == 0.0:
        return 0.0

    return float(np.corrcoef(values_x, values_y)[0, 1])


def run_screen(
    screen_config: ScreenConfig,
    top_n: int = 10,
    verbose: bool = False,
) -> ScreenResult:
    """
    Run all constructs in a screen and evaluate strict criteria.
    """
    summaries: dict[str, ConstructSummary] = {}
    control_types: dict[str, str] = {}

    for entry in screen_config.constructs:
        construct_config = load_construct_config(entry.config_path)
        try:
            summary, _ = run_single_construct(
                construct_config,
                top_n=top_n,
                verbose=verbose,
            )
        except RuntimeError as exc:
            raise RuntimeError(f"{screen_config.screen_id}: {exc}") from exc
        summaries[entry.construct_id] = summary
        control_types[entry.construct_id] = entry.control_type

    candidate_summary = summaries[screen_config.candidate_construct_id]
    baseline_summary = _first_summary_by_control_type(
        summaries,
        control_types,
        "baseline",
    )
    controls = [
        (control_types[construct_id], summary)
        for construct_id, summary in summaries.items()
        if construct_id != screen_config.candidate_construct_id
        and control_types[construct_id] != "baseline"
    ]
    controls_with_baseline = [
        (control_types[construct_id], summary)
        for construct_id, summary in summaries.items()
        if construct_id != screen_config.candidate_construct_id
    ]

    criteria_result = evaluate_criteria(
        candidate_summary,
        baseline_summary,
        controls,
        screen_config.strict_criteria,
    )
    overall_pass = compute_overall_pass(criteria_result)
    control_comparisons = build_control_comparisons(
        candidate_summary,
        controls_with_baseline,
    )
    decision_text = generate_decision_text(
        screen_config.candidate_construct_id,
        criteria_result,
        overall_pass,
        baseline_summary,
    )
    candidate_summary.strict_pass = overall_pass

    return ScreenResult(
        screen_id=screen_config.screen_id,
        candidate_summary=candidate_summary,
        baseline_summary=baseline_summary,
        control_comparisons=control_comparisons,
        criteria_results=criteria_result,
        overall_pass=overall_pass,
        decision_text=decision_text,
    )


def _build_summary(
    construct_config: ConstructConfig,
    pose_metrics: list[PoseMetrics],
    contact_results: list[dict[str, Any]],
) -> ConstructSummary:
    productive_count = sum(1 for pose in pose_metrics if pose.productive)
    best_distance = min(pose.active_site_to_target_distance for pose in pose_metrics)
    best_productive_rank = next(
        (pose.rank for pose in pose_metrics if pose.productive),
        None,
    )
    frequencies = compute_contact_frequencies(contact_results)
    cbd_coupling = pearson_correlation(
        [float(pose.cbd_contacts) for pose in pose_metrics],
        [pose.active_site_to_target_distance for pose in pose_metrics],
    )

    return ConstructSummary(
        construct_id=construct_config.construct_id,
        productive_count=productive_count,
        total_poses=len(pose_metrics),
        best_distance=best_distance,
        best_productive_rank=best_productive_rank,
        lox_contact_frequency=frequencies["lox_contact_frequency"],
        cbd_contact_frequency=frequencies["cbd_contact_frequency"],
        cbd_coupling=cbd_coupling,
        strict_pass=None,
    )


def _domain_tuple(domain: Any) -> tuple[int, int] | None:
    if domain is None:
        return None
    return domain.start, domain.end


def _first_summary_by_control_type(
    summaries: dict[str, ConstructSummary],
    control_types: dict[str, str],
    control_type: str,
) -> ConstructSummary | None:
    for construct_id, summary in summaries.items():
        if control_types[construct_id] == control_type:
            return summary
    return None
