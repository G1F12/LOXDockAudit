"""Distance calculation utilities for LOXDockAudit."""

from __future__ import annotations

from typing import Any

import numpy as np

from loxdockaudit.pdb_parser import get_atoms_by_selection, load_structure


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
) -> dict[str, Any]:
    """
    Compute the minimum heavy-atom distance from LOX active-site residues to target.

    Missing active-site or target selections return structured output with a
    warning instead of raising.
    """
    result: dict[str, Any] = {
        "distance_angstrom": float(np.inf),
        "productive": False,
        "closest_active_site_resi": None,
        "closest_target_atom": None,
        "n_active_site_atoms": 0,
        "n_target_atoms": 0,
        "warning": None,
    }

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
    return result


def score_all_poses(
    pdb_paths: list[str],
    active_site_resi: list[int],
    active_site_chains: list[str],
    target_resi: int,
    target_chain: str,
    threshold: float = 8.0,
    target_atom_names: list[str] | None = None,
) -> list[dict[str, Any]]:
    """
    Run active_site_to_target_distance for each PDB in pdb_paths, preserving order.
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
            )
        except Exception as exc:
            distance_result = _empty_distance_result(
                warning=f"failed to load or score pose: {exc}",
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
    return float(distances[active_index, target_index]), int(active_index), int(target_index)


def _empty_distance_result(warning: str) -> dict[str, Any]:
    return {
        "distance_angstrom": float(np.inf),
        "productive": False,
        "closest_active_site_resi": None,
        "closest_target_atom": None,
        "n_active_site_atoms": 0,
        "n_target_atoms": 0,
        "warning": warning,
    }
