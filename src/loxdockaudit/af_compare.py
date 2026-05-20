"""AlphaFold-Multimer comparison workflow for LOXDockAudit.

The workflow compares productive HDOCK poses against an independent
AlphaFold-Multimer or ColabFold model. All labels are heuristic and
computational; they do not establish enzymatic activity.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from numbers import Real
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from loxdockaudit.af_geometry import AFGeometryResult, GeometryConfig, score_af_model
from loxdockaudit.af_parser import AFMultimerModel, parse_af_json, parse_af_pdb
from loxdockaudit.convergence import ConvergenceReport, compute_interface_overlap
from loxdockaudit.distances import active_site_to_target_distance
from loxdockaudit.interface_extractor import (
    InterfaceResidues,
    contact_pairs_between_chains,
    get_interface_residues,
)
from loxdockaudit.models import OrientationConfig
from loxdockaudit.pdb_parser import load_structure
from loxdockaudit.sorting import natural_sort_key


@dataclass(frozen=True)
class AFCompareTargetResidue:
    """Target substrate residue used for productive geometry scoring."""

    chain: str
    resi: int
    resn: str = ""
    atom: str = "NZ"


@dataclass(frozen=True)
class AFCompareConfig:
    """Configuration for HDOCK versus AF-Multimer comparison."""

    comparison_id: str
    hdock_models_dir: Path
    alphafold_model_path: Path
    target_residue: AFCompareTargetResidue
    active_site_residues: list[int]
    contact_distance_threshold: float
    orientation_enabled: bool
    orientation_threshold_deg: float
    active_site_chain: str = "D"
    productive_distance_threshold: float = 8.0
    top_n: int | None = None
    compare_contact_fingerprints: bool = True
    compare_interface_overlap: bool = True
    compare_productive_geometry: bool = True


@dataclass(frozen=True)
class HDOCKPoseAFMetrics:
    """Per-pose HDOCK metrics used by AF comparison."""

    rank: int
    model_path: str
    distance_to_active_site: float
    productive: bool
    orientation_angle_deg: float | None
    orientation_productive: bool | None
    fully_productive: bool
    interface: InterfaceResidues
    contact_pairs: frozenset[tuple[int, int]]


@dataclass(frozen=True)
class AFCompareResult:
    """Structured result for one HDOCK versus AF-Multimer comparison."""

    comparison_id: str
    hdock_poses: list[HDOCKPoseAFMetrics]
    af_model: AFMultimerModel
    af_geometry: AFGeometryResult
    hdock_productive_interface: InterfaceResidues
    af_interface: InterfaceResidues
    convergence: ConvergenceReport
    convergence_category: str
    contact_persistence: pd.DataFrame


def load_af_compare_config(yaml_path: str | Path) -> AFCompareConfig:
    """Load and validate an ``af-compare`` YAML config."""
    source = Path(yaml_path)
    data = _load_yaml_mapping(source)
    base_dir = source.resolve().parent

    required = [
        "comparison_id",
        "hdock_models_dir",
        "alphafold_model_path",
        "target_residue",
        "active_site_residues",
        "contact_distance_threshold",
        "orientation_enabled",
        "orientation_threshold_deg",
    ]
    missing = [field for field in required if field not in data]
    if missing:
        raise ValueError(f"{source}: missing required field(s): {', '.join(missing)}")

    hdock_models_dir = _resolve_path(base_dir, data["hdock_models_dir"])
    alphafold_model_path = _resolve_path(base_dir, data["alphafold_model_path"])
    if not hdock_models_dir.is_dir():
        raise ValueError(f"{source}: hdock_models_dir is not readable: {hdock_models_dir}")
    if not alphafold_model_path.is_file():
        raise ValueError(
            f"{source}: alphafold_model_path is not readable: {alphafold_model_path}"
        )

    return AFCompareConfig(
        comparison_id=str(data["comparison_id"]),
        hdock_models_dir=hdock_models_dir,
        alphafold_model_path=alphafold_model_path,
        target_residue=_target_from_mapping(data["target_residue"], source),
        active_site_residues=_residue_list(data["active_site_residues"], source),
        contact_distance_threshold=float(data["contact_distance_threshold"]),
        orientation_enabled=bool(data["orientation_enabled"]),
        orientation_threshold_deg=float(data["orientation_threshold_deg"]),
        active_site_chain=str(data.get("active_site_chain", "D")),
        productive_distance_threshold=float(
            data.get("productive_distance_threshold", 8.0)
        ),
        top_n=(
            int(data["top_n"])
            if data.get("top_n") is not None
            else None
        ),
        compare_contact_fingerprints=bool(
            data.get("compare_contact_fingerprints", True)
        ),
        compare_interface_overlap=bool(data.get("compare_interface_overlap", True)),
        compare_productive_geometry=bool(data.get("compare_productive_geometry", True)),
    )


def run_af_compare(
    config: AFCompareConfig,
    out_dir: str | Path,
) -> AFCompareResult:
    """Run HDOCK versus AF-Multimer comparison and write report outputs."""
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    hdock_poses = analyze_hdock_poses(config)
    af_model = _load_af_model(config.alphafold_model_path)
    af_geometry = score_af_model(
        af_model,
        GeometryConfig(
            active_site_residues=config.active_site_residues,
            active_site_chain=config.active_site_chain,
            target_residue=config.target_residue.resi,
            target_chain=config.target_residue.chain,
            target_atom=config.target_residue.atom,
            distance_threshold=config.productive_distance_threshold,
            orientation_enabled=config.orientation_enabled,
            orientation_threshold_deg=config.orientation_threshold_deg,
        ),
    )
    af_interface = get_interface_residues(
        config.alphafold_model_path,
        chain_a=config.target_residue.chain,
        chain_b=config.active_site_chain,
        cutoff_angstrom=config.contact_distance_threshold,
    )
    hdock_interface = _productive_interface(hdock_poses)
    convergence = (
        compute_interface_overlap(hdock_interface, af_interface)
        if config.compare_interface_overlap
        else compute_interface_overlap(
            InterfaceResidues(frozenset(), frozenset()),
            InterfaceResidues(frozenset(), frozenset()),
        )
    )
    category = convergence_category(convergence)
    contact_persistence = build_contact_persistence_table(config, hdock_poses)

    result = AFCompareResult(
        comparison_id=config.comparison_id,
        hdock_poses=hdock_poses,
        af_model=af_model,
        af_geometry=af_geometry,
        hdock_productive_interface=hdock_interface,
        af_interface=af_interface,
        convergence=convergence,
        convergence_category=category,
        contact_persistence=contact_persistence,
    )
    _write_outputs(config, result, output_dir)
    return result


def analyze_hdock_poses(config: AFCompareConfig) -> list[HDOCKPoseAFMetrics]:
    """Score HDOCK poses using the same distance/orientation and interface rules."""
    pdb_paths = sorted(
        config.hdock_models_dir.glob("*.pdb"),
        key=lambda path: natural_sort_key(str(path)),
    )
    if config.top_n is not None:
        pdb_paths = pdb_paths[: config.top_n]
    if not pdb_paths:
        raise ValueError(f"No HDOCK PDB files found in {config.hdock_models_dir}")

    orientation = OrientationConfig(
        enabled=config.orientation_enabled,
        threshold_deg=config.orientation_threshold_deg,
    )
    poses = []
    for rank, pdb_path in enumerate(pdb_paths, start=1):
        structure = load_structure(str(pdb_path), structure_id=f"hdock_{rank}")
        distance_result = active_site_to_target_distance(
            structure=structure,
            active_site_resi=config.active_site_residues,
            active_site_chains=[config.active_site_chain],
            target_resi=config.target_residue.resi,
            target_chain=config.target_residue.chain,
            threshold=config.productive_distance_threshold,
            target_atom_names=[config.target_residue.atom],
            orientation_config=orientation,
            active_site_residues=[
                {"resi": resi, "chain": config.active_site_chain}
                for resi in config.active_site_residues
            ],
        )
        interface = get_interface_residues(
            pdb_path,
            chain_a=config.target_residue.chain,
            chain_b=config.active_site_chain,
            cutoff_angstrom=config.contact_distance_threshold,
        )
        contact_pairs = (
            contact_pairs_between_chains(
                pdb_path,
                chain_a=config.target_residue.chain,
                chain_b=config.active_site_chain,
                cutoff_angstrom=config.contact_distance_threshold,
            )
            if config.compare_contact_fingerprints
            else frozenset()
        )
        productive = bool(distance_result["productive"])
        fully_productive = (
            bool(distance_result.get("fully_productive"))
            if config.orientation_enabled
            else productive
        )
        poses.append(
            HDOCKPoseAFMetrics(
                rank=rank,
                model_path=str(pdb_path),
                distance_to_active_site=float(distance_result["distance_angstrom"]),
                productive=productive,
                orientation_angle_deg=distance_result.get("orientation_angle_deg"),
                orientation_productive=distance_result.get("orientation_productive"),
                fully_productive=fully_productive,
                interface=interface,
                contact_pairs=contact_pairs,
            )
        )
    return poses


def convergence_category(convergence: ConvergenceReport) -> str:
    """Return the heuristic interface convergence category."""
    if convergence.jaccard_index >= 0.5 and convergence.overlap_count > 0:
        return "strong convergence"
    if convergence.overlap_count > 0:
        return "partial convergence"
    return "divergent interfaces"


def build_contact_persistence_table(
    config: AFCompareConfig,
    hdock_poses: list[HDOCKPoseAFMetrics],
) -> pd.DataFrame:
    """Compute residue-pair contact persistence across productive HDOCK poses."""
    productive_poses = [pose for pose in hdock_poses if pose.productive]
    af_pairs = contact_pairs_between_chains(
        config.alphafold_model_path,
        chain_a=config.target_residue.chain,
        chain_b=config.active_site_chain,
        cutoff_angstrom=config.contact_distance_threshold,
    )
    counts: Counter[tuple[int, int]] = Counter()
    for pose in productive_poses:
        counts.update(pose.contact_pairs)

    all_pairs = sorted(set(counts) | set(af_pairs))
    denominator = len(productive_poses)
    rows = [
        {
            "receptor_resnum": receptor_resnum,
            "ligand_resnum": ligand_resnum,
            "productive_pose_count": counts.get((receptor_resnum, ligand_resnum), 0),
            "productive_pose_frequency": (
                counts.get((receptor_resnum, ligand_resnum), 0) / denominator
                if denominator
                else 0.0
            ),
            "present_in_af": (receptor_resnum, ligand_resnum) in af_pairs,
            "shared_contact": (
                counts.get((receptor_resnum, ligand_resnum), 0) > 0
                and (receptor_resnum, ligand_resnum) in af_pairs
            ),
        }
        for receptor_resnum, ligand_resnum in all_pairs
    ]
    return pd.DataFrame(
        rows,
        columns=[
            "receptor_resnum",
            "ligand_resnum",
            "productive_pose_count",
            "productive_pose_frequency",
            "present_in_af",
            "shared_contact",
        ],
    )


def _write_outputs(
    config: AFCompareConfig,
    result: AFCompareResult,
    output_dir: Path,
) -> None:
    summary_df = af_compare_summary_to_dataframe(result)
    overlap_df = interface_overlap_to_dataframe(result)
    summary_df.to_csv(output_dir / f"{config.comparison_id}_af_compare.csv", index=False)
    overlap_df.to_csv(
        output_dir / f"{config.comparison_id}_interface_overlap.csv",
        index=False,
    )
    result.contact_persistence.to_csv(
        output_dir / f"{config.comparison_id}_contact_persistence.csv",
        index=False,
    )
    markdown_report = render_af_compare_markdown(result)
    (output_dir / f"{config.comparison_id}_report.md").write_text(
        markdown_report,
        encoding="utf-8",
    )


def af_compare_summary_to_dataframe(result: AFCompareResult) -> pd.DataFrame:
    """Return one-row summary table for the AF comparison."""
    productive_poses = [pose for pose in result.hdock_poses if pose.productive]
    fully_productive_poses = [
        pose for pose in result.hdock_poses if pose.fully_productive
    ]
    best_distance = min(pose.distance_to_active_site for pose in result.hdock_poses)
    rows = [
        {
            "comparison_id": result.comparison_id,
            "hdock_total_poses": len(result.hdock_poses),
            "hdock_productive_count": len(productive_poses),
            "hdock_fully_productive_count": len(fully_productive_poses),
            "hdock_best_distance": best_distance,
            "hdock_best_productive_rank": _first_rank(productive_poses),
            "hdock_best_fully_productive_rank": _first_rank(fully_productive_poses),
            "af_rank": result.af_geometry.rank,
            "af_iptm": result.af_geometry.iptm,
            "af_distance_to_active_site": result.af_geometry.distance_to_active_site,
            "af_orientation_angle_deg": result.af_geometry.orientation_angle,
            "af_productive": result.af_geometry.is_productive,
            "af_fully_productive": result.af_geometry.is_fully_productive,
            "interface_overlap_count": result.convergence.overlap_count,
            "interface_jaccard_index": result.convergence.jaccard_index,
            "convergence_category": result.convergence_category,
        }
    ]
    return pd.DataFrame(rows)


def interface_overlap_to_dataframe(result: AFCompareResult) -> pd.DataFrame:
    """Return chain-aware residue overlap table."""
    rows = []
    for chain_label, hdock_residues, af_residues in (
        (
            "target_chain",
            result.hdock_productive_interface.chain_a_residues,
            result.af_interface.chain_a_residues,
        ),
        (
            "active_site_chain",
            result.hdock_productive_interface.chain_b_residues,
            result.af_interface.chain_b_residues,
        ),
    ):
        for residue in sorted(hdock_residues | af_residues):
            in_hdock = residue in hdock_residues
            in_af = residue in af_residues
            rows.append(
                {
                    "chain_role": chain_label,
                    "residue": residue,
                    "in_productive_hdock": in_hdock,
                    "in_af": in_af,
                    "classification": _residue_classification(in_hdock, in_af),
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            "chain_role",
            "residue",
            "in_productive_hdock",
            "in_af",
            "classification",
        ],
    )


def render_af_compare_markdown(result: AFCompareResult) -> str:
    """Render a conservative Markdown report for AF comparison."""
    summary = af_compare_summary_to_dataframe(result).iloc[0].to_dict()
    lines = [
        f"# HDOCK vs AF-Multimer Comparison: {result.comparison_id}",
        "",
        "## Executive Summary",
        "",
        (
            f"The comparison shows {result.convergence_category} between "
            "productive HDOCK interfaces and the AF-Multimer interface under "
            "the current geometric metric."
        ),
        "",
        (
            "This is a computational agreement check. It is consistent with a "
            "cross-method structural comparison, but it does not establish "
            "enzymatic activity, catalytic turnover, collagen oxidation, or "
            "biological efficacy."
        ),
        "",
        "## Productive Geometry",
        "",
        f"- HDOCK productive poses: {summary['hdock_productive_count']} / {summary['hdock_total_poses']}",
        f"- HDOCK best productive rank: {_format_optional(summary['hdock_best_productive_rank'])}",
        f"- AF active-site distance: {_format_float(summary['af_distance_to_active_site'])} A",
        f"- AF productive: {summary['af_productive']}",
        f"- AF fully productive: {summary['af_fully_productive']}",
        "",
        "## Interface Convergence",
        "",
        f"- Shared interface residues: {result.convergence.overlap_count}",
        f"- Jaccard index: {result.convergence.jaccard_index:.3f}",
        f"- Category: {result.convergence_category}",
        f"- Narrative: {result.convergence.convergence_narrative}",
        "",
        "## Contact Persistence",
        "",
        _markdown_table_from_dataframe(result.contact_persistence),
        "",
        "## Interpretation Limits",
        "",
        (
            "AlphaFold-Multimer and docking have different objective functions "
            "and failure modes. Partial convergence is useful as independent "
            "computational support for an interface hypothesis, but divergent "
            "interfaces do not by themselves disprove docking and convergent "
            "interfaces do not prove mechanism."
        ),
    ]
    return "\n".join(lines).rstrip() + "\n"


def _productive_interface(poses: list[HDOCKPoseAFMetrics]) -> InterfaceResidues:
    chain_a: set[int] = set()
    chain_b: set[int] = set()
    for pose in poses:
        if not pose.productive:
            continue
        chain_a.update(pose.interface.chain_a_residues)
        chain_b.update(pose.interface.chain_b_residues)
    return InterfaceResidues(frozenset(chain_a), frozenset(chain_b))


def _load_af_model(path: Path) -> AFMultimerModel:
    model = parse_af_pdb(path)
    json_path = path.with_suffix(".json")
    if not json_path.is_file():
        return model
    scores = parse_af_json(json_path)
    return AFMultimerModel(
        pdb_path=model.pdb_path,
        rank=model.rank,
        iptm=scores.iptm,
        ptm=scores.ptm,
        chain_a_residues=model.chain_a_residues,
        chain_b_residues=model.chain_b_residues,
    )


def _target_from_mapping(data: Any, source: Path) -> AFCompareTargetResidue:
    if not isinstance(data, dict):
        raise ValueError(f"{source}: target_residue must be a mapping")
    missing = [field for field in ("chain", "resi") if field not in data]
    if missing:
        raise ValueError(
            f"{source}: missing required field(s) in target_residue: {', '.join(missing)}"
        )
    return AFCompareTargetResidue(
        chain=str(data["chain"]),
        resi=int(data["resi"]),
        resn=str(data.get("resn", "")),
        atom=str(data.get("atom", "NZ")),
    )


def _residue_list(data: Any, source: Path) -> list[int]:
    if not isinstance(data, list) or not data:
        raise ValueError(f"{source}: active_site_residues must be a non-empty list")
    residues = []
    for item in data:
        if isinstance(item, dict):
            if "resi" not in item:
                raise ValueError(f"{source}: active_site_residues entries need resi")
            residues.append(int(item["resi"]))
        else:
            residues.append(int(item))
    return residues


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


def _resolve_path(base_dir: Path, value: Any) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else (base_dir / path).resolve()


def _first_rank(poses: list[HDOCKPoseAFMetrics]) -> int | None:
    return min((pose.rank for pose in poses), default=None)


def _residue_classification(in_hdock: bool, in_af: bool) -> str:
    if in_hdock and in_af:
        return "converged"
    if in_hdock:
        return "hdock_only"
    return "af_only"


def _format_float(value: object) -> str:
    if value is None:
        return "N/A"
    if isinstance(value, Real) and not np.isfinite(float(value)):
        return "inf"
    if isinstance(value, Real):
        return f"{float(value):.3f}"
    return str(value)


def _format_optional(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "N/A"
    return str(value)


def _markdown_table_from_dataframe(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "No contact-persistence rows were generated."
    columns = [str(column) for column in frame.columns]
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    rows = [
        "| " + " | ".join(str(value).replace("|", "\\|") for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *rows])
