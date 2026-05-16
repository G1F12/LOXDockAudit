"""Contact analysis utilities for LOXDockAudit."""

from __future__ import annotations

from typing import Any

import numpy as np

from loxdockaudit.pdb_parser import get_atoms_by_selection


def count_contacts(
    atoms_a: list[Any],
    atoms_b: list[Any],
    cutoff: float = 4.0,
) -> int:
    """
    Count unique atom pairs between atoms_a and atoms_b within cutoff distance.

    Returns 0 if either atom list is empty.
    """
    if not atoms_a or not atoms_b:
        return 0

    contacts = _contact_mask(atoms_a, atoms_b, cutoff)
    return int(np.count_nonzero(contacts))


def get_contact_residues(
    atoms_a: list[Any],
    atoms_b: list[Any],
    cutoff: float = 4.0,
) -> list[str]:
    """
    Return unique residue identifiers from atoms_a contacting any atom in atoms_b.
    """
    if not atoms_a or not atoms_b:
        return []

    contacts = _contact_mask(atoms_a, atoms_b, cutoff)
    contacting_atom_indices = np.flatnonzero(np.any(contacts, axis=1))
    residues = {_residue_identifier(atoms_a[index]) for index in contacting_atom_indices}
    return sorted(residues)


def analyze_pose_contacts(
    structure: Any,
    receptor_chains: list[str],
    ligand_chains: list[str],
    lox_domain_resi_range: tuple[int, int],
    cbd_domain_resi_range: tuple[int, int] | None,
    contact_cutoff: float = 4.0,
) -> dict[str, Any]:
    """
    Compute receptor contacts for ligand, LOX-domain, and CBD/FMOD-domain atoms.
    """
    receptor_atoms = get_atoms_by_selection(structure, chains=receptor_chains)
    ligand_atoms = get_atoms_by_selection(structure, chains=ligand_chains)
    lox_atoms = get_atoms_by_selection(
        structure,
        chains=ligand_chains,
        resi_range=lox_domain_resi_range,
    )

    if cbd_domain_resi_range is None:
        cbd_atoms: list[Any] = []
    else:
        cbd_atoms = get_atoms_by_selection(
            structure,
            chains=ligand_chains,
            resi_range=cbd_domain_resi_range,
        )

    return {
        "total_ligand_contacts": count_contacts(
            ligand_atoms,
            receptor_atoms,
            cutoff=contact_cutoff,
        ),
        "lox_contacts": count_contacts(
            lox_atoms,
            receptor_atoms,
            cutoff=contact_cutoff,
        ),
        "cbd_contacts": count_contacts(
            cbd_atoms,
            receptor_atoms,
            cutoff=contact_cutoff,
        ),
        "lox_contact_residues": get_contact_residues(
            lox_atoms,
            receptor_atoms,
            cutoff=contact_cutoff,
        ),
        "cbd_contact_residues": get_contact_residues(
            cbd_atoms,
            receptor_atoms,
            cutoff=contact_cutoff,
        ),
        "contact_cutoff_angstrom": contact_cutoff,
    }


def compute_contact_frequencies(
    pose_contact_list: list[dict[str, Any]],
    keys: list[str] = ["lox_contacts", "cbd_contacts"],
) -> dict[str, float]:
    """
    Compute fraction of poses where each requested contact count is greater than 0.
    """
    pose_count = len(pose_contact_list)
    frequencies: dict[str, float] = {}

    for key in keys:
        output_key = f"{key.removesuffix('_contacts')}_contact_frequency"
        if pose_count == 0:
            frequencies[output_key] = 0.0
            continue

        contact_pose_count = sum(1 for pose in pose_contact_list if pose.get(key, 0) > 0)
        frequencies[output_key] = contact_pose_count / pose_count

    return frequencies


def _contact_mask(atoms_a: list[Any], atoms_b: list[Any], cutoff: float) -> np.ndarray:
    coords_a = np.array([atom.coord for atom in atoms_a], dtype=float)
    coords_b = np.array([atom.coord for atom in atoms_b], dtype=float)
    deltas = coords_a[:, np.newaxis, :] - coords_b[np.newaxis, :, :]
    distances = np.linalg.norm(deltas, axis=2)
    return distances <= cutoff


def _residue_identifier(atom: Any) -> str:
    residue = atom.get_parent()
    chain = residue.get_parent()
    return f"{chain.id}_{residue.id[1]}_{residue.resname}"
