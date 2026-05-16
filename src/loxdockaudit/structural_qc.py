"""Structural quality-control checks for LOX active-site models."""

from __future__ import annotations

from itertools import combinations
from typing import Any

import numpy as np

from loxdockaudit import pdb_parser


def check_his_triad_geometry(
    structure: Any,
    his_resi: list[int],
    chain: str,
    reference_ca_rmsd_threshold: float = 2.0,
) -> dict:
    """
    Check His triad geometry using CA positions.

    The reference_ca_rmsd_threshold argument is reserved for compatibility with
    future reference-structure comparisons.
    """
    del reference_ca_rmsd_threshold

    result: dict[str, Any] = {
        "his_resi": list(his_resi),
        "ca_coords_found": [],
        "ca_coords_missing": list(his_resi),
        "pairwise_ca_distances": {},
        "max_ca_distance": 0.0,
        "geometry_plausible": False,
        "warning": None,
    }

    try:
        coords: dict[int, np.ndarray] = {}
        for resi in his_resi:
            atom = _get_atom(structure, chain, resi, "CA")
            if atom is not None:
                coords[resi] = _coord(atom)

        found = [resi for resi in his_resi if resi in coords]
        missing = [resi for resi in his_resi if resi not in coords]
        result["ca_coords_found"] = found
        result["ca_coords_missing"] = missing

        distances: dict[str, float] = {}
        for resi_a, resi_b in combinations(his_resi, 2):
            if resi_a in coords and resi_b in coords:
                distances[f"{resi_a}-{resi_b}"] = _distance(coords[resi_a], coords[resi_b])

        result["pairwise_ca_distances"] = distances
        result["max_ca_distance"] = max(distances.values(), default=0.0)

        if len(found) != 3 or len(his_resi) != 3:
            result["warning"] = "His triad CA atoms missing."
        elif result["max_ca_distance"] > 10.0:
            result["warning"] = "His triad CA distances exceed expected range."
        else:
            result["geometry_plausible"] = True

    except Exception as exc:  # pragma: no cover - defensive API boundary
        result["geometry_plausible"] = False
        result["warning"] = f"His triad geometry check failed: {exc}"

    return result


def check_lys_tyr_geometry(
    structure: Any,
    lys_resi: int,
    tyr_resi: int,
    chain: str,
    max_cb_distance: float = 12.0,
) -> dict:
    """Check Lys320-Tyr355 spatial relationship using CB atoms."""
    result: dict[str, Any] = {
        "lys_resi": lys_resi,
        "tyr_resi": tyr_resi,
        "lys_cb_found": False,
        "tyr_cb_found": False,
        "cb_distance": None,
        "geometry_plausible": False,
        "warning": None,
    }

    try:
        lys_atom, lys_cb_found = _get_atom_with_fallback(structure, chain, lys_resi, "CB", "CA")
        tyr_atom, tyr_cb_found = _get_atom_with_fallback(structure, chain, tyr_resi, "CB", "CA")
        result["lys_cb_found"] = lys_cb_found
        result["tyr_cb_found"] = tyr_cb_found

        if lys_atom is None or tyr_atom is None:
            missing = []
            if lys_atom is None:
                missing.append(f"Lys {lys_resi}")
            if tyr_atom is None:
                missing.append(f"Tyr {tyr_resi}")
            result["warning"] = "Missing CB and CA atom for " + ", ".join(missing) + "."
            return result

        distance = _distance(_coord(lys_atom), _coord(tyr_atom))
        result["cb_distance"] = distance
        result["geometry_plausible"] = distance <= max_cb_distance

        warnings = []
        if not lys_cb_found or not tyr_cb_found:
            warnings.append("CB missing; used CA fallback.")
        if not result["geometry_plausible"]:
            warnings.append("Lys-Tyr distance exceeds expected range.")
        result["warning"] = " ".join(warnings) if warnings else None

    except Exception as exc:  # pragma: no cover - defensive API boundary
        result["geometry_plausible"] = False
        result["warning"] = f"Lys-Tyr geometry check failed: {exc}"

    return result


def check_disulfide_geometry(
    structure: Any,
    disulfide_pairs: list[tuple[int, int]],
    chain: str,
    max_sg_distance: float = 2.5,
) -> dict:
    """Check disulfide bond geometry using SG atom distances."""
    result: dict[str, Any] = {
        "pairs_checked": len(disulfide_pairs),
        "pairs_intact": 0,
        "pair_results": [],
        "all_intact": False,
        "warning": None,
    }

    try:
        pair_results = []
        for cys_a, cys_b in disulfide_pairs:
            atom_a, sg_a_found = _get_atom_with_fallback(structure, chain, cys_a, "SG", "CB")
            atom_b, sg_b_found = _get_atom_with_fallback(structure, chain, cys_b, "SG", "CB")
            pair_result: dict[str, Any] = {
                "cys_a": cys_a,
                "cys_b": cys_b,
                "sg_distance": None,
                "intact": False,
                "warning": None,
            }

            if atom_a is None or atom_b is None:
                missing = []
                if atom_a is None:
                    missing.append(f"Cys {cys_a}")
                if atom_b is None:
                    missing.append(f"Cys {cys_b}")
                pair_result["warning"] = "Missing SG and CB atom for " + ", ".join(missing) + "."
                pair_results.append(pair_result)
                continue

            distance = _distance(_coord(atom_a), _coord(atom_b))
            pair_result["sg_distance"] = distance
            pair_result["intact"] = distance <= max_sg_distance

            warnings = []
            if not sg_a_found or not sg_b_found:
                warnings.append("SG missing; used CB fallback.")
            if not pair_result["intact"]:
                warnings.append("Disulfide distance exceeds expected range.")
            pair_result["warning"] = " ".join(warnings) if warnings else None
            pair_results.append(pair_result)

        pairs_intact = sum(1 for pair in pair_results if pair["intact"])
        result["pair_results"] = pair_results
        result["pairs_intact"] = pairs_intact
        result["all_intact"] = pairs_intact == len(disulfide_pairs)
        if not result["all_intact"]:
            result["warning"] = "One or more disulfide pairs are broken or missing."

    except Exception as exc:  # pragma: no cover - defensive API boundary
        result["all_intact"] = False
        result["warning"] = f"Disulfide geometry check failed: {exc}"

    return result


def check_active_site_sasa_proxy(
    structure: Any,
    active_site_resi: list[int],
    chain: str,
    contact_cutoff: float = 5.0,
) -> dict:
    """
    Compute a simple burial proxy from non-active-site heavy-atom contacts.

    This is not true SASA.
    """
    result: dict[str, Any] = {
        "active_site_resi": list(active_site_resi),
        "n_active_site_atoms": 0,
        "n_surrounding_atoms_within_cutoff": 0,
        "burial_score": 0.0,
        "accessibility_flag": "accessible",
        "warning": None,
    }

    try:
        active_set = set(active_site_resi)
        active_atoms = pdb_parser.get_atoms_by_selection(
            structure,
            chains=[chain],
            resi_list=list(active_site_resi),
        )
        surrounding_atoms = [
            atom
            for atom in pdb_parser.get_atoms_by_selection(structure, chains=[chain])
            if _residue_number(atom) not in active_set
        ]

        result["n_active_site_atoms"] = len(active_atoms)
        if not active_atoms:
            result["warning"] = "No active-site heavy atoms found."
            return result

        if not surrounding_atoms:
            result["warning"] = "No surrounding heavy atoms found in chain."
            return result

        active_coords = np.array([_coord(atom) for atom in active_atoms], dtype=float)
        surrounding_coords = np.array([_coord(atom) for atom in surrounding_atoms], dtype=float)
        deltas = surrounding_coords[:, np.newaxis, :] - active_coords[np.newaxis, :, :]
        distances = np.linalg.norm(deltas, axis=2)
        surrounding_count = int(np.any(distances <= contact_cutoff, axis=1).sum())
        burial_score = surrounding_count / len(active_atoms)

        result["n_surrounding_atoms_within_cutoff"] = surrounding_count
        result["burial_score"] = float(burial_score)
        result["accessibility_flag"] = _accessibility_flag(burial_score)

    except Exception as exc:  # pragma: no cover - defensive API boundary
        result["warning"] = f"Active-site SASA proxy check failed: {exc}"

    return result


def _get_atom(structure: Any, chain_id: str, resi: int, atom_name: str) -> Any | None:
    for atom in pdb_parser.get_atoms_by_selection(
        structure,
        chains=[chain_id],
        resi_list=[resi],
        atom_names=[atom_name],
    ):
        return atom
    return None


def _get_atom_with_fallback(
    structure: Any,
    chain_id: str,
    resi: int,
    primary_atom_name: str,
    fallback_atom_name: str,
) -> tuple[Any | None, bool]:
    primary_atom = _get_atom(structure, chain_id, resi, primary_atom_name)
    if primary_atom is not None:
        return primary_atom, True
    return _get_atom(structure, chain_id, resi, fallback_atom_name), False


def _coord(atom: Any) -> np.ndarray:
    return np.asarray(atom.get_coord(), dtype=float)


def _residue_number(atom: Any) -> int:
    return int(atom.get_parent().id[1])


def _distance(coord_a: np.ndarray, coord_b: np.ndarray) -> float:
    return float(np.linalg.norm(coord_a - coord_b))


def _accessibility_flag(burial_score: float) -> str:
    if burial_score < 5.0:
        return "accessible"
    if burial_score < 15.0:
        return "partially_buried"
    return "buried"
