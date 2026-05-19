"""Active-vs-inactive catalytic-control comparison utilities."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from loxdockaudit.config import load_construct_config
from loxdockaudit.contacts import analyze_pose_contacts
from loxdockaudit.distances import PoseGeometry  # v0.4
from loxdockaudit.fold_qc import format_qc_report, run_structural_qc
from loxdockaudit.input_audit import audit_screen, format_audit_report
from loxdockaudit.models import (
    ConstructConfig,
    ConstructSummary,
    OrientationConfig,
    StructuralQCResult,
)
from loxdockaudit.pdb_parser import load_structure
from loxdockaudit.reporting import (  # v0.4
    generate_geometry_scatter,
    poses_to_dataframe,
    summaries_to_dataframe,
)
from loxdockaudit.screen import run_single_construct
from loxdockaudit.sorting import natural_sort_key


@dataclass(frozen=True)
class CatalyticResidueExpectation:
    label: str
    resi: int
    expected_resn: str


@dataclass(frozen=True)
class InactiveControlSide:
    label: str
    construct_id: str
    config_path: str
    catalytic_residues: list[CatalyticResidueExpectation]


@dataclass(frozen=True)
class InactiveControlConfig:
    comparison_id: str
    active: InactiveControlSide
    inactive: InactiveControlSide
    mutations: list[str]
    top_n: int = 10
    productive_distance_threshold: float = 8.0
    orientation: OrientationConfig = field(default_factory=OrientationConfig)  # v0.4


def load_inactive_control_config(yaml_path: str) -> InactiveControlConfig:
    """Load the active-vs-inactive comparison YAML."""
    source = Path(yaml_path)
    data = _load_yaml_mapping(source)
    base_dir = source.resolve().parent

    missing = [
        field
        for field in ("comparison_id", "active", "inactive")
        if field not in data
    ]
    if missing:
        raise ValueError(f"{yaml_path}: missing required field(s): {', '.join(missing)}")

    return InactiveControlConfig(
        comparison_id=str(data["comparison_id"]),
        active=_side_from_mapping(data["active"], base_dir, f"{yaml_path}:active"),
        inactive=_side_from_mapping(data["inactive"], base_dir, f"{yaml_path}:inactive"),
        mutations=[str(mutation) for mutation in data.get("mutations", [])],
        top_n=int(data.get("top_n", 10)),
        productive_distance_threshold=float(
            data.get("productive_distance_threshold", 8.0)
        ),
        orientation=_orientation_from_mapping(data.get("orientation")),  # v0.4
    )


def run_inactive_control_analysis(
    comparison_config: InactiveControlConfig,
    out_dir: str | Path,
    top_n: int | None = None,
) -> dict[str, Any]:
    """Run active and inactive constructs through the same comparison path."""
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    n_poses = top_n if top_n is not None else comparison_config.top_n

    active_config = load_construct_config(comparison_config.active.config_path)
    inactive_config = load_construct_config(comparison_config.inactive.config_path)
    # v0.4: allow inactive-control configs to enable orientation for both sides.
    if comparison_config.orientation.enabled:
        active_config.orientation = comparison_config.orientation
        inactive_config.orientation = comparison_config.orientation

    audit_result = audit_screen(
        [
            _audit_config_from_construct(active_config),
            _audit_config_from_construct(inactive_config),
        ]
    )
    audit_report = format_audit_report(audit_result)
    (output_dir / f"{comparison_config.comparison_id}_input_audit.md").write_text(
        audit_report,
        encoding="utf-8",
    )
    if not audit_result["screen_valid"]:
        raise RuntimeError("inactive-control input audit failed")

    active_summary, active_poses = run_single_construct(active_config, top_n=n_poses)
    inactive_summary, inactive_poses = run_single_construct(inactive_config, top_n=n_poses)
    # v0.4: save distance/orientation scatter plots for each side.
    generate_geometry_scatter(
        _pose_geometries_from_metrics(active_poses),
        str(output_dir),
        active_config.construct_id,
    )
    generate_geometry_scatter(
        _pose_geometries_from_metrics(inactive_poses),
        str(output_dir),
        inactive_config.construct_id,
    )

    active_pose_df = poses_to_dataframe(active_poses)
    inactive_pose_df = poses_to_dataframe(inactive_poses)
    active_summary_df = summaries_to_dataframe([active_summary])
    inactive_summary_df = summaries_to_dataframe([inactive_summary])

    active_qc = run_configured_structural_qc(active_config)
    inactive_qc = run_configured_structural_qc(inactive_config)

    active_catalytic = evaluate_catalytic_identity(
        _first_model_path(active_config.models_dir),
        active_config.ligand_chains[0],
        comparison_config.active.catalytic_residues,
    )
    inactive_catalytic = evaluate_catalytic_identity(
        _first_model_path(inactive_config.models_dir),
        inactive_config.ligand_chains[0],
        comparison_config.inactive.catalytic_residues,
    )

    comparison_df = build_comparison_table(
        active_summary=active_summary,
        inactive_summary=inactive_summary,
        active_qc=active_qc,
        inactive_qc=inactive_qc,
        active_catalytic=active_catalytic,
        inactive_catalytic=inactive_catalytic,
    )
    distance_df = build_distance_distribution(
        active_config.construct_id,
        active_pose_df,
        inactive_config.construct_id,
        inactive_pose_df,
    )
    contact_df = build_contact_persistence_table(
        {
            active_config.construct_id: active_config,
            inactive_config.construct_id: inactive_config,
        },
        top_n=n_poses,
    )
    cluster_df = build_pose_cluster_table(
        {
            active_config.construct_id: active_config,
            inactive_config.construct_id: inactive_config,
        },
        top_n=n_poses,
    )

    interpretation = interpret_inactive_control(active_summary, inactive_summary)
    supplement = generate_inactive_control_supplement(
        comparison_config=comparison_config,
        active_config=active_config,
        inactive_config=inactive_config,
        comparison_df=comparison_df,
        interpretation=interpretation,
    )
    svg = generate_distance_histogram_svg(
        distance_df,
        threshold=comparison_config.productive_distance_threshold,
    )

    active_pose_df.to_csv(output_dir / f"{active_config.construct_id}_poses.csv", index=False)
    inactive_pose_df.to_csv(
        output_dir / f"{inactive_config.construct_id}_poses.csv",
        index=False,
    )
    active_summary_df.to_csv(
        output_dir / f"{active_config.construct_id}_summary.csv",
        index=False,
    )
    inactive_summary_df.to_csv(
        output_dir / f"{inactive_config.construct_id}_summary.csv",
        index=False,
    )
    comparison_df.to_csv(
        output_dir / f"{comparison_config.comparison_id}_comparison.csv",
        index=False,
    )
    distance_df.to_csv(
        output_dir / f"{comparison_config.comparison_id}_distance_distribution.csv",
        index=False,
    )
    contact_df.to_csv(
        output_dir / f"{comparison_config.comparison_id}_contact_persistence.csv",
        index=False,
    )
    cluster_df.to_csv(
        output_dir / f"{comparison_config.comparison_id}_pose_clusters.csv",
        index=False,
    )
    (output_dir / f"{active_config.construct_id}_structural_qc.md").write_text(
        format_qc_report(active_qc),
        encoding="utf-8",
    )
    (output_dir / f"{inactive_config.construct_id}_structural_qc.md").write_text(
        format_qc_report(inactive_qc),
        encoding="utf-8",
    )
    (output_dir / f"{comparison_config.comparison_id}_supplement.md").write_text(
        supplement,
        encoding="utf-8",
    )
    (output_dir / f"{comparison_config.comparison_id}_distance_histogram.svg").write_text(
        svg,
        encoding="utf-8",
    )

    return {
        "active_summary": active_summary,
        "inactive_summary": inactive_summary,
        "active_qc": active_qc,
        "inactive_qc": inactive_qc,
        "active_catalytic": active_catalytic,
        "inactive_catalytic": inactive_catalytic,
        "comparison_df": comparison_df,
        "distance_df": distance_df,
        "contact_df": contact_df,
        "cluster_df": cluster_df,
        "interpretation": interpretation,
    }


def evaluate_catalytic_identity(
    pdb_path: str | Path,
    chain_id: str,
    expectations: list[CatalyticResidueExpectation],
) -> dict[str, Any]:
    """Check catalytic residue identities without treating mutation as fold failure."""
    structure = load_structure(str(pdb_path))
    found: list[dict[str, Any]] = []
    all_expected = True

    model = next(structure.get_models(), None)
    chain = model[chain_id] if model is not None and chain_id in model else None
    for expectation in expectations:
        actual_resn = None
        if chain is not None:
            residue = chain.child_dict.get((" ", expectation.resi, " "))
            if residue is not None:
                actual_resn = residue.resname.strip().upper()
        expected_resn = expectation.expected_resn.strip().upper()
        matches = actual_resn == expected_resn
        all_expected = all_expected and matches
        found.append(
            {
                "label": expectation.label,
                "resi": expectation.resi,
                "expected_resn": expected_resn,
                "actual_resn": actual_resn,
                "matches_expected": matches,
            }
        )

    expected_names = {item.expected_resn.strip().upper() for item in expectations}
    if all_expected and expected_names == {"HIS"}:
        catalytic_status = "intact"
    elif all_expected:
        catalytic_status = "disrupted"
    else:
        catalytic_status = "unexpected"

    return {
        "pdb_path": str(pdb_path),
        "chain_id": chain_id,
        "expected_identity_pass": all_expected,
        "catalytic_status": catalytic_status,
        "residues": found,
    }


def build_comparison_table(
    *,
    active_summary: ConstructSummary,
    inactive_summary: ConstructSummary,
    active_qc: StructuralQCResult,
    inactive_qc: StructuralQCResult,
    active_catalytic: dict[str, Any],
    inactive_catalytic: dict[str, Any],
) -> pd.DataFrame:
    """Return a wide active-vs-inactive comparison table."""
    rows = [
        (
            "Productive poses",
            f"{active_summary.productive_count}/{active_summary.total_poses}",
            f"{inactive_summary.productive_count}/{inactive_summary.total_poses}",
        ),
        (
            "Best distance",
            _format_float(active_summary.best_distance),
            _format_float(inactive_summary.best_distance),
        ),
        (
            "Best productive rank",
            _format_optional(active_summary.best_productive_rank),
            _format_optional(inactive_summary.best_productive_rank),
        ),
        # v0.4: orientation-aware summary rows.
        (
            "Fully productive poses",
            str(active_summary.fully_productive_count),
            str(inactive_summary.fully_productive_count),
        ),
        (
            "Best orientation angle",
            _format_optional_na(active_summary.best_orientation_angle_deg),
            _format_optional_na(inactive_summary.best_orientation_angle_deg),
        ),
        (
            "Best fully productive rank",
            _format_optional_na(active_summary.best_fully_productive_rank),
            _format_optional_na(inactive_summary.best_fully_productive_rank),
        ),
        (
            "Productive fraction",
            _format_float(_productive_fraction(active_summary)),
            _format_float(_productive_fraction(inactive_summary)),
        ),
        (
            "QC status",
            "PASS" if active_qc.fold_qc_pass else "FOLD CORRUPTED",
            "PASS" if inactive_qc.fold_qc_pass else "FOLD CORRUPTED",
        ),
        (
            "Active-site access",
            "PASS" if active_qc.active_site_accessible else "FAIL",
            "PASS" if inactive_qc.active_site_accessible else "FAIL",
        ),
        (
            "Active-site geometry",
            str(active_catalytic["catalytic_status"]),
            str(inactive_catalytic["catalytic_status"]),
        ),
    ]
    frame = pd.DataFrame(rows, columns=["Metric", "Active", "Inactive"])
    # v0.4: machine-readable compatibility column for smoke checks.
    frame["fully_productive_count"] = ""
    frame.loc[
        frame["Metric"] == "Fully productive poses",
        "fully_productive_count",
    ] = (
        f"{active_summary.fully_productive_count}/"
        f"{inactive_summary.fully_productive_count}"
    )
    return frame


def build_distance_distribution(
    active_id: str,
    active_pose_df: pd.DataFrame,
    inactive_id: str,
    inactive_pose_df: pd.DataFrame,
) -> pd.DataFrame:
    """Build a normalized per-pose distance table for plotting and review."""
    frames = []
    for construct_id, pose_df in (
        (active_id, active_pose_df),
        (inactive_id, inactive_pose_df),
    ):
        current = pose_df[
            [
                "rank",
                "active_site_to_target_distance",
                "productive",
                "orientation_angle_deg",  # v0.4
                "orientation_productive",  # v0.4
                "fully_productive",  # v0.4
                "closest_active_site_resi",
            ]
        ].copy()
        current.insert(0, "construct_id", construct_id)
        frames.append(current)
    return pd.concat(frames, ignore_index=True)


def build_contact_persistence_table(
    construct_configs: dict[str, ConstructConfig],
    top_n: int,
) -> pd.DataFrame:
    """Compute how often each LOX residue contacts collagen across top poses."""
    rows: list[dict[str, Any]] = []
    for construct_id, construct_config in construct_configs.items():
        pose_contacts = _contact_fingerprints(construct_config, top_n)
        total = len(pose_contacts)
        by_residue: dict[str, list[int]] = {}
        for rank, residues in pose_contacts:
            for residue in residues:
                by_residue.setdefault(residue, []).append(rank)

        for residue, ranks in sorted(by_residue.items()):
            rows.append(
                {
                    "construct_id": construct_id,
                    "contact_residue": residue,
                    "pose_count": len(ranks),
                    "frequency": len(ranks) / total if total else 0.0,
                    "ranks": ";".join(str(rank) for rank in ranks),
                }
            )
    return pd.DataFrame(
        rows,
        columns=["construct_id", "contact_residue", "pose_count", "frequency", "ranks"],
    )


def build_pose_cluster_table(
    construct_configs: dict[str, ConstructConfig],
    top_n: int,
) -> pd.DataFrame:
    """Cluster poses by exact LOX contact-residue fingerprint."""
    rows: list[dict[str, Any]] = []
    for construct_id, construct_config in construct_configs.items():
        pose_contacts = _contact_fingerprints(construct_config, top_n)
        total = len(pose_contacts)
        clusters: dict[tuple[str, ...], list[int]] = {}
        for rank, residues in pose_contacts:
            clusters.setdefault(tuple(residues), []).append(rank)

        sorted_clusters = sorted(
            clusters.items(),
            key=lambda item: (-len(item[1]), item[1][0]),
        )
        for index, (fingerprint, ranks) in enumerate(sorted_clusters, start=1):
            rows.append(
                {
                    "construct_id": construct_id,
                    "cluster_id": f"{construct_id}_cluster_{index}",
                    "pose_count": len(ranks),
                    "frequency": len(ranks) / total if total else 0.0,
                    "ranks": ";".join(str(rank) for rank in ranks),
                    "contact_fingerprint": ";".join(fingerprint),
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            "construct_id",
            "cluster_id",
            "pose_count",
            "frequency",
            "ranks",
            "contact_fingerprint",
        ],
    )


def interpret_inactive_control(
    active_summary: ConstructSummary,
    inactive_summary: ConstructSummary,
) -> str:
    """Return the pre-specified scientific interpretation text."""
    active_fraction = _productive_fraction(active_summary)
    inactive_fraction = _productive_fraction(inactive_summary)

    if (
        active_summary.productive_count > inactive_summary.productive_count
        and active_summary.best_distance <= inactive_summary.best_distance
    ):
        return (
            "Active LOX169-417 retained more favorable productive geometry under "
            "the current docking metric. This is not proof of enzymatic activity."
        )

    similar_fraction = abs(active_fraction - inactive_fraction) <= 0.1
    similar_distance = (
        math.isfinite(active_summary.best_distance)
        and math.isfinite(inactive_summary.best_distance)
        and abs(active_summary.best_distance - inactive_summary.best_distance) <= 0.5
    )
    if similar_fraction and similar_distance:
        return (
            "Active and inactive results were similar. Docking geometry alone may "
            "not discriminate catalytic competence in this setup."
        )

    if (
        inactive_summary.productive_count > active_summary.productive_count
        or inactive_summary.best_distance < active_summary.best_distance
    ):
        return (
            "The inactive control showed equal or more favorable geometry. The "
            "current metric may overweight surface complementarity because "
            "active-site chemistry is absent from docking scoring."
        )

    return (
        "The inactive control did not produce a simple active-better pattern. "
        "Interpret the result as a structural-control outcome, not as activity proof."
    )


def generate_inactive_control_supplement(
    *,
    comparison_config: InactiveControlConfig,
    active_config: ConstructConfig,
    inactive_config: ConstructConfig,
    comparison_df: pd.DataFrame,
    interpretation: str,
) -> str:
    """Generate the one-page Markdown supplement."""
    mutation_text = ", ".join(comparison_config.mutations) or "not specified"
    lines = [
        "# Active vs Catalytically Inactive Control Comparison",
        "",
        "## Design",
        "",
        (
            "The inactive control was introduced to test whether productive "
            "geometry metrics were sensitive to catalytic-site disruption."
        ),
        "",
        f"- Active construct: {active_config.construct_id}",
        f"- Inactive construct: {inactive_config.construct_id}",
        f"- Mutations: {mutation_text}",
        (
            "- Interpretation scope: computational structural-control experiment, "
            "not enzymatic activity validation."
        ),
        "",
        "## Results",
        "",
        _markdown_table_from_dataframe(comparison_df),
        "",
        "## Orientation score interpretation",
        "",
        (
            "Orientation scoring reports whether the target Lys side chain points "
            "toward the active-site centroid. Lower angles indicate a more "
            "productive approach geometry; a fully productive pose must satisfy "
            "both the distance cutoff and the orientation cutoff when orientation "
            "scoring is enabled."
        ),
        "",
        "## Interpretation",
        "",
        interpretation,
        "",
        "## Limitations",
        "",
        (
            "The comparison uses docking-derived productive geometry and simple "
            "structural QC. HDOCK scoring does not include copper coordination, "
            "LTQ chemistry, catalytic turnover, or wet-lab activity."
        ),
    ]
    return "\n".join(lines).rstrip() + "\n"


def generate_distance_histogram_svg(
    distance_df: pd.DataFrame,
    threshold: float,
    bin_width: float = 2.0,
) -> str:
    """Generate distance and orientation histograms without plot dependencies."""
    finite = distance_df[
        distance_df["active_site_to_target_distance"].apply(math.isfinite)
    ].copy()
    if finite.empty:
        return "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"760\" height=\"620\"></svg>\n"

    max_distance = max(float(finite["active_site_to_target_distance"].max()), threshold)
    max_bin = int(math.ceil(max_distance / bin_width))
    bins = [(index * bin_width, (index + 1) * bin_width) for index in range(max_bin)]
    constructs = list(dict.fromkeys(str(value) for value in finite["construct_id"]))
    colors = ["#2563eb", "#f97316"]  # v0.4: active blue, inactive orange.
    counts: dict[tuple[str, int], int] = {}
    max_count = 1

    for construct in constructs:
        values = finite.loc[
            finite["construct_id"] == construct,
            "active_site_to_target_distance",
        ]
        for value in values:
            bin_index = min(int(float(value) // bin_width), max_bin - 1)
            key = (construct, bin_index)
            counts[key] = counts.get(key, 0) + 1
            max_count = max(max_count, counts[key])

    width = 760
    height = 620
    left = 64
    right = 24
    top = 44
    panel_gap = 76
    panel_height = 210
    plot_width = width - left - right
    plot_height = panel_height
    group_width = plot_width / max_bin
    bar_width = max(4.0, group_width / (len(constructs) + 1))

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="64" y="24" font-family="Arial" font-size="16" font-weight="700">Active-site distance distribution</text>',
        f'<line x1="{left}" y1="{top + plot_height}" x2="{width - right}" y2="{top + plot_height}" stroke="#111827"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#111827"/>',
    ]

    for index, (start, end) in enumerate(bins):
        x0 = left + index * group_width
        label = f"{int(start)}-{int(end)}"
        parts.append(
            f'<text x="{x0 + group_width / 2:.1f}" y="{height - 34}" text-anchor="middle" font-family="Arial" font-size="10">{label}</text>'
        )
        for construct_index, construct in enumerate(constructs):
            count = counts.get((construct, index), 0)
            bar_height = (count / max_count) * plot_height if count else 0
            x = x0 + 6 + construct_index * bar_width
            y = top + plot_height - bar_height
            parts.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width - 2:.1f}" height="{bar_height:.1f}" fill="{colors[construct_index % len(colors)]}"/>'
            )
            if count:
                parts.append(
                    f'<text x="{x + (bar_width - 2) / 2:.1f}" y="{y - 4:.1f}" text-anchor="middle" font-family="Arial" font-size="10">{count}</text>'
                )

    threshold_x = left + (threshold / (max_bin * bin_width)) * plot_width
    parts.extend(
        [
            f'<line x1="{threshold_x:.1f}" y1="{top}" x2="{threshold_x:.1f}" y2="{top + plot_height}" stroke="#111827" stroke-dasharray="4 4"/>',
            f'<text x="{threshold_x + 4:.1f}" y="{top + 12}" font-family="Arial" font-size="11">threshold {threshold:.1f} A</text>',
            f'<text x="{left + plot_width / 2:.1f}" y="{top + plot_height + 42}" text-anchor="middle" font-family="Arial" font-size="12">Distance bin (A)</text>',
            f'<text x="18" y="{top + plot_height / 2:.1f}" transform="rotate(-90 18 {top + plot_height / 2:.1f})" text-anchor="middle" font-family="Arial" font-size="12">Pose count</text>',
        ]
    )

    # v0.4: second panel for orientation angle distribution.
    orientation_top = top + plot_height + panel_gap
    parts.extend(
        _orientation_histogram_panel(
            distance_df,
            constructs,
            colors,
            left,
            orientation_top,
            plot_width,
            panel_height,
            width,
            right,
        )
    )

    legend_x = width - right - 260
    for index, construct in enumerate(constructs):
        y = 18 + index * 18
        parts.append(
            f'<rect x="{legend_x}" y="{y - 10}" width="10" height="10" fill="{colors[index % len(colors)]}"/>'
        )
        parts.append(
            f'<text x="{legend_x + 16}" y="{y}" font-family="Arial" font-size="11">{construct}</text>'
        )

    parts.append("</svg>")
    return "\n".join(parts) + "\n"


# v0.4: orientation angle histogram panel using 18-degree bins.
def _orientation_histogram_panel(
    distance_df: pd.DataFrame,
    constructs: list[str],
    colors: list[str],
    left: int,
    top: int,
    plot_width: int,
    plot_height: int,
    width: int,
    right: int,
) -> list[str]:
    bins = [(index * 18.0, (index + 1) * 18.0) for index in range(10)]
    group_width = plot_width / len(bins)
    bar_width = max(4.0, group_width / (len(constructs) + 1))
    counts: dict[tuple[str, int], int] = {}
    max_count = 1

    for _, row in distance_df.iterrows():
        angle = pd.to_numeric(
            pd.Series([row.get("orientation_angle_deg")]),
            errors="coerce",
        ).iloc[0]
        if pd.isna(angle):
            continue
        construct = str(row["construct_id"])
        bin_index = min(int(float(angle) // 18.0), len(bins) - 1)
        key = (construct, bin_index)
        counts[key] = counts.get(key, 0) + 1
        max_count = max(max_count, counts[key])

    parts = [
        f'<text x="64" y="{top - 16}" font-family="Arial" font-size="16" font-weight="700">Orientation angle distribution</text>',
        f'<line x1="{left}" y1="{top + plot_height}" x2="{width - right}" y2="{top + plot_height}" stroke="#111827"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#111827"/>',
    ]
    for index, (start, end) in enumerate(bins):
        x0 = left + index * group_width
        label = f"{int(start)}-{int(end)}"
        parts.append(
            f'<text x="{x0 + group_width / 2:.1f}" y="{top + plot_height + 18}" text-anchor="middle" font-family="Arial" font-size="10">{label}</text>'
        )
        for construct_index, construct in enumerate(constructs):
            count = counts.get((construct, index), 0)
            bar_height = (count / max_count) * plot_height if count else 0
            x = x0 + 6 + construct_index * bar_width
            y = top + plot_height - bar_height
            parts.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width - 2:.1f}" height="{bar_height:.1f}" fill="{colors[construct_index % len(colors)]}"/>'
            )
            if count:
                parts.append(
                    f'<text x="{x + (bar_width - 2) / 2:.1f}" y="{y - 4:.1f}" text-anchor="middle" font-family="Arial" font-size="10">{count}</text>'
                )
    parts.extend(
        [
            f'<text x="{left + plot_width / 2:.1f}" y="{top + plot_height + 42}" text-anchor="middle" font-family="Arial" font-size="12">Orientation angle bin (deg)</text>',
            f'<text x="18" y="{top + plot_height / 2:.1f}" transform="rotate(-90 18 {top + plot_height / 2:.1f})" text-anchor="middle" font-family="Arial" font-size="12">Pose count</text>',
        ]
    )
    return parts


def run_configured_structural_qc(construct_config: ConstructConfig) -> StructuralQCResult:
    """Run structural QC using the construct config's structural_qc block."""
    qc_config = construct_config.structural_qc
    if qc_config is None:
        raise ValueError(f"{construct_config.construct_id}: structural_qc is required")

    active_site_resi = [
        int(residue["resi"])
        for residue in construct_config.active_site.residues
        if "resi" in residue
    ]
    return run_structural_qc(
        pdb_path=str(_first_model_path(construct_config.models_dir)),
        construct_id=construct_config.construct_id,
        lox_chain=construct_config.ligand_chains[0],
        his_resi=qc_config.his_resi,
        lys_resi=qc_config.lys_resi,
        tyr_resi=qc_config.tyr_resi,
        disulfide_pairs=qc_config.disulfide_pairs,
        active_site_resi=active_site_resi,
        alphafold_pdb_path=qc_config.alphafold_pdb_path,
        pae_json_path=qc_config.pae_json_path,
        domain_a_resi=qc_config.domain_a_resi,
        domain_b_resi=qc_config.domain_b_resi,
        plddt_threshold=qc_config.plddt_threshold,
        his_max_ca_distance=qc_config.his_max_ca_distance,
        lys_tyr_max_cb_distance=qc_config.lys_tyr_max_cb_distance,
        disulfide_max_sg_distance=qc_config.disulfide_max_sg_distance,
    )


def _contact_fingerprints(
    construct_config: ConstructConfig,
    top_n: int,
) -> list[tuple[int, list[str]]]:
    pdb_paths = sorted(
        Path(construct_config.models_dir).glob("*.pdb"),
        key=lambda path: natural_sort_key(str(path)),
    )[:top_n]
    fingerprints: list[tuple[int, list[str]]] = []
    for rank, pdb_path in enumerate(pdb_paths, start=1):
        structure = load_structure(str(pdb_path), structure_id=f"pose_{rank}")
        contact_result = analyze_pose_contacts(
            structure=structure,
            receptor_chains=construct_config.receptor_chains,
            ligand_chains=construct_config.ligand_chains,
            lox_domain_resi_range=(
                construct_config.lox_domain.start,
                construct_config.lox_domain.end,
            ),
            cbd_domain_resi_range=(
                (construct_config.cbd_domain.start, construct_config.cbd_domain.end)
                if construct_config.cbd_domain is not None
                else None
            ),
            contact_cutoff=construct_config.contact_distance_threshold,
        )
        residues = sorted(str(residue) for residue in contact_result["lox_contact_residues"])
        fingerprints.append((rank, residues))
    return fingerprints


def _side_from_mapping(
    data: Any,
    base_dir: Path,
    source: str,
) -> InactiveControlSide:
    if not isinstance(data, dict):
        raise ValueError(f"{source}: must be a mapping")
    missing = [
        field
        for field in ("construct_id", "config_path", "catalytic_residues")
        if field not in data
    ]
    if missing:
        raise ValueError(f"{source}: missing required field(s): {', '.join(missing)}")

    raw_config_path = Path(str(data["config_path"]))
    config_path = raw_config_path if raw_config_path.is_absolute() else base_dir / raw_config_path
    return InactiveControlSide(
        label=str(data.get("label", data["construct_id"])),
        construct_id=str(data["construct_id"]),
        config_path=str(config_path),
        catalytic_residues=[
            _expectation_from_mapping(item, f"{source}:catalytic_residues[{index}]")
            for index, item in enumerate(data["catalytic_residues"])
        ],
    )


def _orientation_from_mapping(data: Any) -> OrientationConfig:
    """Parse optional v0.4 inactive-control orientation override."""
    if data is None:
        return OrientationConfig()
    if not isinstance(data, dict):
        raise ValueError("orientation must be a mapping")
    threshold_deg = float(data.get("threshold_deg", 90.0))
    if not 0.0 <= threshold_deg <= 180.0:
        raise ValueError("orientation.threshold_deg must be between 0.0 and 180.0")
    return OrientationConfig(
        enabled=bool(data.get("enabled", False)),
        threshold_deg=threshold_deg,
        sidechain_atoms=[str(value) for value in data.get("sidechain_atoms", ["CB", "NZ"])],
        activesite_centroid_atoms=[
            str(value)
            for value in data.get("activesite_centroid_atoms", ["CA"])
        ],
    )


def _expectation_from_mapping(data: Any, source: str) -> CatalyticResidueExpectation:
    if not isinstance(data, dict):
        raise ValueError(f"{source}: must be a mapping")
    missing = [field for field in ("label", "resi", "expected_resn") if field not in data]
    if missing:
        raise ValueError(f"{source}: missing required field(s): {', '.join(missing)}")
    return CatalyticResidueExpectation(
        label=str(data["label"]),
        resi=int(data["resi"]),
        expected_resn=str(data["expected_resn"]).strip().upper(),
    )


def _load_yaml_mapping(path: Path) -> dict[str, Any]:
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ValueError(f"Could not read YAML file {path}: {exc}") from exc
    except yaml.YAMLError as exc:
        raise ValueError(f"Could not parse YAML file {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError(f"YAML file must contain a mapping: {path}")
    return data


def _audit_config_from_construct(construct_config: ConstructConfig) -> dict[str, Any]:
    return {
        "construct_id": construct_config.construct_id,
        "models_dir": construct_config.models_dir,
        "ligand_chains": construct_config.ligand_chains,
        "receptor_chains": construct_config.receptor_chains,
        "expected_n_terminus": None,
    }


def _first_model_path(models_dir: str) -> Path:
    pdb_paths = sorted(
        Path(models_dir).glob("*.pdb"),
        key=lambda path: natural_sort_key(str(path)),
    )
    if not pdb_paths:
        raise ValueError(f"No PDB files found in models_dir: {models_dir}")
    return pdb_paths[0]


def _productive_fraction(summary: ConstructSummary) -> float:
    return summary.productive_count / summary.total_poses if summary.total_poses else 0.0


def _pose_geometries_from_metrics(poses: list[Any]) -> list[PoseGeometry]:
    return [
        PoseGeometry(
            pose_rank=pose.rank,
            distance_A=pose.active_site_to_target_distance,
            distance_productive=pose.productive,
            orientation_angle_deg=pose.orientation_angle_deg,
            orientation_productive=pose.orientation_productive,
            fully_productive=(
                bool(pose.fully_productive)
                if pose.fully_productive is not None
                else pose.productive
            ),
        )
        for pose in poses
    ]


def _format_float(value: float) -> str:
    if not math.isfinite(value):
        return "inf"
    return f"{value:.3f}"


def _format_optional(value: object) -> str:
    return "None" if value is None else str(value)


def _format_optional_na(value: object) -> str:
    if value is None:
        return "N/A"
    if isinstance(value, float) and math.isnan(value):
        return "N/A"
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def _markdown_table_from_dataframe(frame: pd.DataFrame) -> str:
    columns = [str(column) for column in frame.columns]
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    rows = [
        "| " + " | ".join(str(value).replace("|", "\\|") for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *rows])
