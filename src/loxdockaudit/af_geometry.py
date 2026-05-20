"""Apply LOXDockAudit geometry metrics to AF-Multimer models."""

from __future__ import annotations

from dataclasses import dataclass

from loxdockaudit.af_parser import AFMultimerModel
from loxdockaudit.distances import active_site_to_target_distance
from loxdockaudit.models import OrientationConfig
from loxdockaudit.pdb_parser import load_structure


@dataclass(frozen=True)
class GeometryConfig:
    """Configuration for scoring AF-Multimer productive geometry."""

    active_site_residues: list[int]
    active_site_chain: str
    target_residue: int
    target_chain: str
    target_atom: str = "NZ"
    distance_threshold: float = 8.0
    orientation_enabled: bool = False
    orientation_threshold_deg: float = 90.0


@dataclass(frozen=True)
class AFGeometryResult:
    """Productive-geometry result for one AF-Multimer model."""

    rank: int
    iptm: float | None
    distance_to_active_site: float | None
    orientation_angle: float | None
    is_productive: bool
    is_fully_productive: bool


def score_af_model(
    model: AFMultimerModel,
    config: GeometryConfig,
) -> AFGeometryResult:
    """Score one AF-Multimer model using LOXDockAudit geometry metrics."""
    structure = load_structure(str(model.pdb_path), structure_id=f"af_{model.rank}")
    orientation = OrientationConfig(
        enabled=config.orientation_enabled,
        threshold_deg=config.orientation_threshold_deg,
    )
    distance_result = active_site_to_target_distance(
        structure=structure,
        active_site_resi=config.active_site_residues,
        active_site_chains=[config.active_site_chain],
        target_resi=config.target_residue,
        target_chain=config.target_chain,
        threshold=config.distance_threshold,
        target_atom_names=[config.target_atom],
        orientation_config=orientation,
        active_site_residues=[
            {"resi": resi, "chain": config.active_site_chain}
            for resi in config.active_site_residues
        ],
    )
    productive = bool(distance_result["productive"])
    fully_productive = (
        bool(distance_result.get("fully_productive"))
        if config.orientation_enabled
        else productive
    )
    return AFGeometryResult(
        rank=model.rank,
        iptm=model.iptm,
        distance_to_active_site=float(distance_result["distance_angstrom"]),
        orientation_angle=distance_result.get("orientation_angle_deg"),
        is_productive=productive,
        is_fully_productive=fully_productive,
    )


def score_af_run(
    models: list[AFMultimerModel],
    config: GeometryConfig,
) -> list[AFGeometryResult]:
    """Score all AF-Multimer models from a run."""
    return [score_af_model(model, config) for model in sorted(models, key=lambda item: item.rank)]
