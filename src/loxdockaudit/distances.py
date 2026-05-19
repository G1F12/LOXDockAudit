"""Distance calculation utilities for LOXDockAudit."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from loxdockaudit.models import OrientationConfig
from loxdockaudit.orientation import compute_orientation_score
from loxdockaudit.pdb_parser import get_atoms_by_selection, load_structure


# v0.4: combined distance/orientation result for a ranked docking pose.
@dataclass
class PoseGeometry:
    """Productive-geometry summary for one docking pose."""

    pose_rank: int
    distance_A: float
    distance_productive: bool
    orientation_angle_deg: float | None
    orientation_productive: bool | None
    fully_productive: bool


def min_heavy_atom_distance(
    atoms_a: list[Any],
    atoms_b: list[Any],
) -> float:
    """
    Return minimum pairwise distance between any atom in atoms_a and atoms_b.

    Returns ``np.inf`` if either atom list is empty.
    """
    if not atoms_a or not atoms_b:
        return float(np.inf)

    coords_a = np.array([atom.coord for atom in atoms_a], dtype=float)
    coords_b = np.array([atom.coord for atom in atoms_b], dtype=float)
    deltas = coords_a[:, np.newaxis, :] - coords_b[np.newaxis, :, :]
    distances = np.linalg.norm(deltas, axis=2)
    return float(np.min(distances))


def active_site_to_target_distance(
    structure: Any,
    active_site_resi: list[int],
    active_site_chains: list[str],
    target_resi: int,
    target_chain: str,
    threshold: float = 8.0,
    target_atom_names: list[str] | None = None,
    orientation_config: OrientationConfig | None = None,  # v0.4
    active_site_residues: list[dict[str, Any]] | None = None,  # v0.4
) -> dict[str, Any]:
    """
    Compute the minimum heavy-atom distance from LOX active-site residues to target.

    Missing active-site or target selections return structured output with a
    warning instead of raising. v0.4 orientation scoring is opt-in through
    ``orientation_config``; when it is omitted, the returned dictionary is
    identical to the pre-v0.4 distance-only output.
    """
    # v0.4: preserve exact v0.3 output unless an orientation config is supplied.
    include_orientation = orientation_config is not None
    orientation = orientation_config or OrientationConfig(enabled=False)

    result: dict[str, Any] = {
        "distance_angstrom": float(np.inf),
        "productive": False,
        "closest_active_site_resi": None,
        "closest_target_atom": None,
        "n_active_site_atoms": 0,
        "n_target_atoms": 0,
        "warning": None,
    }
    if include_orientation:
        result.update(
            {
                "orientation_angle_deg": None,
                "orientation_productive": None,
                "fully_productive": False,
                "orientation_error": None,
            }
        )

    try:
        active_site_atoms = get_atoms_by_selection(
            structure,
            chains=active_site_chains,
            resi_list=active_site_resi,
        )
        target_atoms = get_atoms_by_selection(
            structure,
            chains=[target_chain],
            resi_list=[target_resi],
            atom_names=target_atom_names,
        )
    except Exception as exc:
        result["warning"] = f"atom selection failed: {exc}"
        if include_orientation:
            result["orientation_error"] = result["warning"]
        return result

    result["n_active_site_atoms"] = len(active_site_atoms)
    result["n_target_atoms"] = len(target_atoms)

    warnings: list[str] = []
    if not active_site_atoms:
        warnings.append("active site residues not found")
    if not target_atoms:
        warnings.append("target residue not found")
    if warnings:
        result["warning"] = "; ".join(warnings)
        if include_orientation:
            result["orientation_error"] = result["warning"]
        return result

    distance, active_index, target_index = _closest_atom_pair(
        active_site_atoms,
        target_atoms,
    )
    closest_active_atom = active_site_atoms[active_index]
    closest_target_atom = target_atoms[target_index]
    closest_active_residue = closest_active_atom.get_parent()

    result["distance_angstrom"] = distance
    result["productive"] = distance <= threshold
    result["closest_active_site_resi"] = (
        closest_active_residue.id[1] if closest_active_residue is not None else None
    )
    result["closest_target_atom"] = closest_target_atom.name

    # v0.4: add orientation on top of the existing distance calculation.
    if include_orientation:
        if orientation.enabled:
            pose_atoms = _structure_to_atom_dicts(structure)
            centroid_residues = active_site_residues or _active_site_residue_dicts(
                active_site_resi,
                active_site_chains,
            )
            orientation_result = compute_orientation_score(
                pose_atoms=pose_atoms,
                active_site_residues=centroid_residues,
                target_resi=target_resi,
                target_chain=target_chain,
                threshold_deg=orientation.threshold_deg,
                enabled=True,
            )
            result["orientation_angle_deg"] = orientation_result.angle_deg
            result["orientation_productive"] = orientation_result.productive
            result["orientation_error"] = orientation_result.error
            result["fully_productive"] = bool(
                result["productive"] and orientation_result.productive
            )
        else:
            result["fully_productive"] = result["productive"]

    return result


# v0.4: dataclass-facing pose analysis API for combined geometry scoring.
def analyze_pose_geometry(
    structure: Any,
    pose_rank: int,
    active_site_resi: list[int],
    active_site_chains: list[str],
    target_resi: int,
    target_chain: str,
    threshold: float = 8.0,
    target_atom_names: list[str] | None = None,
    orientation_config: OrientationConfig | None = None,
    active_site_residues: list[dict[str, Any]] | None = None,
) -> PoseGeometry:
    """
    Analyze distance and optional orientation geometry for one ranked pose.

    If ``orientation_config`` is omitted or disabled, orientation fields are
    ``None`` and ``fully_productive`` is identical to ``distance_productive``.
    The underlying distance calculation is delegated to
    ``active_site_to_target_distance`` unchanged.
    """
    orientation = orientation_config or OrientationConfig(enabled=False)
    result = active_site_to_target_distance(
        structure=structure,
        active_site_resi=active_site_resi,
        active_site_chains=active_site_chains,
        target_resi=target_resi,
        target_chain=target_chain,
        threshold=threshold,
        target_atom_names=target_atom_names,
        orientation_config=orientation,
        active_site_residues=active_site_residues,
    )

    distance_productive = bool(result["productive"])
    fully_productive = (
        bool(result.get("fully_productive"))
        if orientation.enabled
        else distance_productive
    )

    return PoseGeometry(
        pose_rank=pose_rank,
        distance_A=float(result["distance_angstrom"]),
        distance_productive=distance_productive,
        orientation_angle_deg=result.get("orientation_angle_deg"),
        orientation_productive=result.get("orientation_productive"),
        fully_productive=fully_productive,
    )


def score_all_poses(
    pdb_paths: list[str],
    active_site_resi: list[int],
    active_site_chains: list[str],
    target_resi: int,
    target_chain: str,
    threshold: float = 8.0,
    target_atom_names: list[str] | None = None,
    orientation_config: OrientationConfig | None = None,  # v0.4
    active_site_residues: list[dict[str, Any]] | None = None,  # v0.4
) -> list[dict[str, Any]]:
    """
    Run active_site_to_target_distance for each PDB in pdb_paths, preserving order.

    v0.4 orientation fields are included only when ``orientation_config`` is
    passed, preserving pre-v0.4 output for existing callers.
    """
    scores: list[dict[str, Any]] = []
    for rank, pdb_path in enumerate(pdb_paths, start=1):
        try:
            structure = load_structure(pdb_path, structure_id=f"pose_{rank}")
            distance_result = active_site_to_target_distance(
                structure=structure,
                active_site_resi=active_site_resi,
                active_site_chains=active_site_chains,
                target_resi=target_resi,
                target_chain=target_chain,
                threshold=threshold,
                target_atom_names=target_atom_names,
                orientation_config=orientation_config,  # v0.4
                active_site_residues=active_site_residues,  # v0.4
            )
        except Exception as exc:
            distance_result = _empty_distance_result(
                warning=f"failed to load or score pose: {exc}",
                include_orientation=orientation_config is not None,  # v0.4
            )

        scores.append(
            {
                "model_path": pdb_path,
                "rank": rank,
                **distance_result,
            }
        )

    return scores


def _closest_atom_pair(atoms_a: list[Any], atoms_b: list[Any]) -> tuple[float, int, int]:
    coords_a = np.array([atom.coord for atom in atoms_a], dtype=float)
    coords_b = np.array([atom.coord for atom in atoms_b], dtype=float)
    deltas = coords_a[:, np.newaxis, :] - coords_b[np.newaxis, :, :]
    distances = np.linalg.norm(deltas, axis=2)
    flat_index = int(np.argmin(distances))
    active_index, target_index = np.unravel_index(flat_index, distances.shape)
    return (
        float(distances[active_index, target_index]),
        int(active_index),
        int(target_index),
    )


def _empty_distance_result(
    warning: str,
    include_orientation: bool = False,  # v0.4
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "distance_angstrom": float(np.inf),
        "productive": False,
        "closest_active_site_resi": None,
        "closest_target_atom": None,
        "n_active_site_atoms": 0,
        "n_target_atoms": 0,
        "warning": warning,
    }
    # v0.4: failed opt-in orientation scoring still has the v0.4 fields.
    if include_orientation:
        result.update(
            {
                "orientation_angle_deg": None,
                "orientation_productive": None,
                "fully_productive": False,
                "orientation_error": warning,
            }
        )
    return result


# v0.4: convert Biopython structures into the atom-dict shape used by orientation.py.
def _structure_to_atom_dicts(structure: Any) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atoms.append(
                        {
                            "atom_name": atom.name,
                            "resname": residue.resname,
                            "chain_id": chain.id,
                            "resseq": residue.id[1],
                            "x": float(atom.coord[0]),
                            "y": float(atom.coord[1]),
                            "z": float(atom.coord[2]),
                        }
                    )
    return atoms


# v0.4: synthesize active-site residue dictionaries when callers only pass v0.3 lists.
def _active_site_residue_dicts(
    active_site_resi: list[int],
    active_site_chains: list[str],
) -> list[dict[str, Any]]:
    residues: list[dict[str, Any]] = []
    for resi in active_site_resi:
        if active_site_chains:
            for chain in active_site_chains:
                residues.append({"resi": resi, "chain": chain})
        else:
            residues.append({"resi": resi})
    return residues
