"""Interface residue extraction for HDOCK and AF-Multimer PDB structures."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from loxdockaudit.pdb_parser import get_atoms_by_selection, load_structure


@dataclass(frozen=True)
class InterfaceResidues:
    """Residue numbers participating in a two-chain interface."""

    chain_a_residues: frozenset[int]
    chain_b_residues: frozenset[int]


@dataclass(frozen=True)
class ContactPair:
    """Closest active-site/substrate residue contact pair."""

    receptor_resnum: int
    ligand_resnum: int
    distance: float


def get_interface_residues(
    pdb_path: Path,
    chain_a: str,
    chain_b: str,
    cutoff_angstrom: float = 5.0,
) -> InterfaceResidues:
    """Return interface residues between two chains within ``cutoff_angstrom``."""
    structure = load_structure(str(pdb_path))
    atoms_a = get_atoms_by_selection(structure, chains=[chain_a])
    atoms_b = get_atoms_by_selection(structure, chains=[chain_b])
    residues_a, residues_b = _interface_sets(atoms_a, atoms_b, cutoff_angstrom)
    return InterfaceResidues(
        chain_a_residues=frozenset(residues_a),
        chain_b_residues=frozenset(residues_b),
    )


def get_active_site_contacts(
    pdb_path: Path,
    active_site_residues: list[int],
    substrate_chain: str,
) -> list[ContactPair]:
    """Return active-site contacts to a substrate chain."""
    structure = load_structure(str(pdb_path))
    substrate_atoms = get_atoms_by_selection(structure, chains=[substrate_chain])
    active_atoms = [
        atom
        for atom in get_atoms_by_selection(
            structure,
            chains=_non_substrate_chains(structure, substrate_chain),
            resi_list=active_site_residues,
        )
    ]
    contacts = []
    for active_atom in active_atoms:
        closest_atom, distance = _closest_atom(active_atom, substrate_atoms)
        if closest_atom is None:
            continue
        contacts.append(
            ContactPair(
                receptor_resnum=_resi(active_atom),
                ligand_resnum=_resi(closest_atom),
                distance=distance,
            )
        )
    return sorted(
        contacts,
        key=lambda contact: (
            contact.receptor_resnum,
            contact.ligand_resnum,
            contact.distance,
        ),
    )


def contact_pairs_between_chains(
    pdb_path: Path,
    chain_a: str,
    chain_b: str,
    cutoff_angstrom: float = 5.0,
) -> frozenset[tuple[int, int]]:
    """Return residue-pair contacts between two chains within ``cutoff_angstrom``."""
    structure = load_structure(str(pdb_path))
    atoms_a = get_atoms_by_selection(structure, chains=[chain_a])
    atoms_b = get_atoms_by_selection(structure, chains=[chain_b])
    if not atoms_a or not atoms_b:
        return frozenset()

    distances = _distance_matrix(atoms_a, atoms_b)
    pairs = {
        (_resi(atoms_a[index_a]), _resi(atoms_b[index_b]))
        for index_a, index_b in zip(*np.where(distances <= cutoff_angstrom))
    }
    return frozenset(sorted(pairs))


def _interface_sets(
    atoms_a: list[Any],
    atoms_b: list[Any],
    cutoff_angstrom: float,
) -> tuple[set[int], set[int]]:
    if not atoms_a or not atoms_b:
        return set(), set()
    distances = _distance_matrix(atoms_a, atoms_b)
    indices_a, indices_b = np.where(distances <= cutoff_angstrom)
    return (
        {_resi(atoms_a[index]) for index in indices_a},
        {_resi(atoms_b[index]) for index in indices_b},
    )


def _distance_matrix(atoms_a: list[Any], atoms_b: list[Any]) -> np.ndarray:
    coords_a = np.array([atom.coord for atom in atoms_a], dtype=float)
    coords_b = np.array([atom.coord for atom in atoms_b], dtype=float)
    deltas = coords_a[:, np.newaxis, :] - coords_b[np.newaxis, :, :]
    return np.linalg.norm(deltas, axis=2)


def _closest_atom(atom: Any, candidates: list[Any]) -> tuple[Any | None, float]:
    if not candidates:
        return None, float("inf")
    coords = np.array([candidate.coord for candidate in candidates], dtype=float)
    distances = np.linalg.norm(coords - np.asarray(atom.coord, dtype=float), axis=1)
    index = int(np.argmin(distances))
    return candidates[index], float(distances[index])


def _non_substrate_chains(structure: Any, substrate_chain: str) -> list[str]:
    model = next(structure.get_models(), None)
    if model is None:
        return []
    return [chain.id for chain in model if chain.id != substrate_chain]


def _resi(atom: Any) -> int:
    return int(atom.get_parent().id[1])
