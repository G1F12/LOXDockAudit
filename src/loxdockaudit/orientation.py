"""Orientation scoring utilities for LOX-collagen docking poses."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass
class OrientationResult:
    """Structured result for a target-Lys orientation calculation."""

    angle_deg: float | None
    productive: bool | None
    error: str | None


def get_atom_coords(atoms: list[dict[str, Any]], atom_name: str) -> np.ndarray:
    """Return xyz coordinates for the first atom matching ``atom_name``.

    Parameters
    ----------
    atoms
        Atom dictionaries containing an atom-name field and either ``x``,
        ``y``, ``z`` fields or a three-value coordinate field.
    atom_name
        PDB atom name to find, such as ``"CA"``, ``"CB"``, or ``"NZ"``.

    Returns
    -------
    np.ndarray
        A length-three floating-point coordinate array.

    Raises
    ------
    ValueError
        If no matching atom is present or if the matching atom has invalid
        coordinates.
    """
    for atom in atoms:
        if _atom_name(atom) == atom_name:
            return _coords(atom)
    raise ValueError(f"atom not found: {atom_name}")


def compute_activesite_centroid(
    pose_atoms: list[dict[str, Any]],
    active_site_residues: list[dict[str, Any]],
    centroid_atom: str = "CA",
) -> np.ndarray:
    """Return the centroid of configured active-site residue atoms.

    The active-site residues are taken from configuration dictionaries. When
    labeled histidine entries such as His292, His294, and His296 are present,
    those entries define the centroid. Otherwise all configured residue entries
    are used, which keeps the function usable with minimal fixtures. The
    calculation uses the first matching ``centroid_atom`` from each configured
    residue. If a residue dictionary also contains a chain field, that chain is
    respected.

    Parameters
    ----------
    pose_atoms
        Atom dictionaries for the pose.
    active_site_residues
        Configured active-site residue dictionaries. Each entry must contain a
        residue number under ``resi`` or ``resseq``.
    centroid_atom
        Atom name to use for each residue centroid point, defaulting to ``CA``.

    Returns
    -------
    np.ndarray
        A length-three floating-point centroid coordinate.

    Raises
    ------
    ValueError
        If no configured active-site residue atom can be found.
    """
    centroid_residues = _centroid_residue_entries(active_site_residues)
    residue_specs = [_residue_spec(residue) for residue in centroid_residues]
    coords: list[np.ndarray] = []

    for atom in pose_atoms:
        if _atom_name(atom) != centroid_atom:
            continue
        atom_resi = _resi(atom)
        atom_chain = _chain(atom)
        for spec_resi, spec_chain in residue_specs:
            if atom_resi == spec_resi and (spec_chain is None or atom_chain == spec_chain):
                coords.append(_coords(atom))
                break

    if not coords:
        raise ValueError("fewer than 1 active-site atom found")

    return np.mean(np.stack(coords, axis=0), axis=0)


def compute_lys_sidechain_vector(
    pose_atoms: list[dict[str, Any]],
    target_resi: int,
    target_chain: str,
) -> np.ndarray:
    """Return the unit side-chain vector from target Lys/Hyl C-beta to NZ.

    Parameters
    ----------
    pose_atoms
        Atom dictionaries for the pose.
    target_resi
        Residue number of the target Lys or hydroxy-Lys.
    target_chain
        Chain identifier for the target residue.

    Returns
    -------
    np.ndarray
        Length-three unit vector pointing from ``CB`` to ``NZ``.

    Raises
    ------
    ValueError
        If either ``CB`` or ``NZ`` is missing, coordinates are invalid, or the
        vector has zero length.
    """
    target_atoms = [
        atom
        for atom in pose_atoms
        if _resi(atom) == target_resi and _chain(atom) == target_chain
    ]
    cb_coords = get_atom_coords(target_atoms, "CB")
    nz_coords = get_atom_coords(target_atoms, "NZ")
    return _unit_vector(nz_coords - cb_coords, "Lys side-chain vector")


def compute_approach_vector(centroid: np.ndarray, nz_coords: np.ndarray) -> np.ndarray:
    """Return the unit approach vector from active-site centroid to target NZ.

    Parameters
    ----------
    centroid
        Active-site centroid coordinates.
    nz_coords
        Target Lys/Hyl ``NZ`` coordinates.

    Returns
    -------
    np.ndarray
        Length-three unit vector pointing from ``centroid`` to ``nz_coords``.

    Raises
    ------
    ValueError
        If the coordinates are invalid or the vector has zero length.
    """
    return _unit_vector(_as_xyz(nz_coords) - _as_xyz(centroid), "approach vector")


def compute_orientation_angle(sidechain_vec: np.ndarray, approach_vec: np.ndarray) -> float:
    """Return the angle in degrees between two vectors.

    Both input vectors are normalized before computing the angle. The dot
    product is clamped to ``[-1, 1]`` before ``arccos`` to avoid NaN results
    from small floating-point roundoff errors.

    Parameters
    ----------
    sidechain_vec
        Vector describing the target side-chain direction.
    approach_vec
        Vector describing the active-site-to-target approach direction.

    Returns
    -------
    float
        Angle in degrees in the inclusive range ``0.0`` to ``180.0``.

    Raises
    ------
    ValueError
        If either vector has invalid coordinates or zero length.
    """
    unit_sidechain = _unit_vector(sidechain_vec, "side-chain vector")
    unit_approach = _unit_vector(approach_vec, "approach vector")
    dot_product = float(np.dot(unit_sidechain, unit_approach))
    clamped_dot = float(np.clip(dot_product, -1.0, 1.0))
    return float(np.degrees(np.arccos(clamped_dot)))


def is_orientation_productive(angle_deg: float, threshold_deg: float = 90.0) -> bool:
    """Return whether ``angle_deg`` is productive under an inclusive threshold.

    Parameters
    ----------
    angle_deg
        Orientation angle in degrees.
    threshold_deg
        Maximum productive angle in degrees. The comparison is inclusive.

    Returns
    -------
    bool
        ``True`` when ``angle_deg <= threshold_deg``.
    """
    return angle_deg <= threshold_deg


def compute_orientation_score(
    pose_atoms: list[dict[str, Any]],
    active_site_residues: list[dict[str, Any]],
    target_resi: int,
    target_chain: str,
    threshold_deg: float = 90.0,
    centroid_atom: str = "CA",
    enabled: bool = True,
) -> OrientationResult:
    """Compute the configured target-Lys orientation score for one pose.

    This convenience wrapper returns an ``OrientationResult`` instead of
    raising. When ``enabled`` is ``False``, all result fields are ``None`` so
    callers can cleanly propagate disabled orientation scoring from config.

    Parameters
    ----------
    pose_atoms
        Atom dictionaries for the pose.
    active_site_residues
        Configured active-site residue dictionaries used for the centroid.
    target_resi
        Residue number of the target Lys or hydroxy-Lys.
    target_chain
        Chain identifier for the target residue.
    threshold_deg
        Inclusive productive-orientation threshold in degrees.
    centroid_atom
        Active-site atom name used to build the centroid.
    enabled
        Whether orientation scoring should be performed.

    Returns
    -------
    OrientationResult
        Angle/productive fields when scoring succeeds, an error message when it
        fails, or all ``None`` values when scoring is disabled.
    """
    if not enabled:
        return OrientationResult(angle_deg=None, productive=None, error=None)

    try:
        centroid = compute_activesite_centroid(
            pose_atoms,
            active_site_residues,
            centroid_atom=centroid_atom,
        )
        target_atoms = [
            atom
            for atom in pose_atoms
            if _resi(atom) == target_resi and _chain(atom) == target_chain
        ]
        nz_coords = get_atom_coords(target_atoms, "NZ")
        sidechain_vec = compute_lys_sidechain_vector(
            pose_atoms,
            target_resi=target_resi,
            target_chain=target_chain,
        )
        approach_vec = compute_approach_vector(centroid, nz_coords)
        angle_deg = compute_orientation_angle(sidechain_vec, approach_vec)
    except ValueError as exc:
        return OrientationResult(angle_deg=None, productive=None, error=str(exc))

    return OrientationResult(
        angle_deg=angle_deg,
        productive=is_orientation_productive(angle_deg, threshold_deg),
        error=None,
    )


def _atom_name(atom: dict[str, Any]) -> str | None:
    name = atom.get("atom_name", atom.get("name"))
    return str(name).strip() if name is not None else None


def _chain(atom: dict[str, Any]) -> str | None:
    chain = atom.get("chain_id", atom.get("chain"))
    return str(chain).strip() if chain is not None else None


def _resi(atom: dict[str, Any]) -> int | None:
    resi = atom.get("resseq", atom.get("resi"))
    return int(resi) if resi is not None else None


def _residue_spec(residue: dict[str, Any]) -> tuple[int, str | None]:
    resi = residue.get("resi", residue.get("resseq"))
    if resi is None:
        raise ValueError("active-site residue is missing resi")

    chain = residue.get("chain_id", residue.get("chain"))
    chain_id = str(chain).strip() if chain is not None else None
    return int(resi), chain_id or None


def _centroid_residue_entries(
    active_site_residues: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    histidine_entries = [
        residue
        for residue in active_site_residues
        if str(residue.get("label", "")).strip().upper().startswith("HIS")
    ]
    return histidine_entries if histidine_entries else active_site_residues


def _coords(atom: dict[str, Any]) -> np.ndarray:
    if all(key in atom for key in ("x", "y", "z")):
        return _as_xyz(np.array([atom["x"], atom["y"], atom["z"]], dtype=float))

    for key in ("coord", "coords", "xyz"):
        if key in atom:
            return _as_xyz(atom[key])

    raise ValueError(f"atom has no xyz coordinates: {_atom_name(atom)}")


def _as_xyz(coords: Any) -> np.ndarray:
    array = np.asarray(coords, dtype=float)
    if array.shape != (3,):
        raise ValueError("coordinates must be a length-3 vector")
    if not np.all(np.isfinite(array)):
        raise ValueError("coordinates must be finite")
    return array


def _unit_vector(vector: np.ndarray, label: str) -> np.ndarray:
    vector = _as_xyz(vector)
    norm = float(np.linalg.norm(vector))
    if norm == 0.0:
        raise ValueError(f"{label} has zero length")
    return vector / norm
