from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from loxdockaudit.contacts import (
    analyze_pose_contacts,
    compute_contact_frequencies,
    count_contacts,
    get_contact_residues,
)
from loxdockaudit.pdb_parser import load_structure


@dataclass
class Chain:
    id: str


@dataclass
class Residue:
    id: tuple[str, int, str]
    resname: str
    chain: Chain

    def get_parent(self) -> Chain:
        return self.chain


@dataclass
class Atom:
    coord: np.ndarray
    residue: Residue

    def get_parent(self) -> Residue:
        return self.residue


def atom(x: float, resi: int = 124, resname: str = "HIS") -> Atom:
    return Atom(
        coord=np.array([x, 0.0, 0.0], dtype=float),
        residue=Residue((" ", resi, " "), resname, Chain("D")),
    )


def test_count_contacts_within_cutoff() -> None:
    assert count_contacts([atom(0)], [atom(3)], cutoff=4.0) == 1


def test_count_contacts_outside_cutoff() -> None:
    assert count_contacts([atom(0)], [atom(5)], cutoff=4.0) == 0


def test_count_contacts_empty_atoms_a() -> None:
    assert count_contacts([], [atom(0)]) == 0


def test_get_contact_residues_format() -> None:
    assert get_contact_residues([atom(0)], [atom(3)], cutoff=4.0) == ["D_124_HIS"]


def test_get_contact_residues_no_duplicates() -> None:
    atoms_a = [atom(0), atom(1)]
    assert get_contact_residues(atoms_a, [atom(3)], cutoff=4.0) == ["D_124_HIS"]


def test_analyze_pose_contacts_fixture(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = analyze_pose_contacts(
        structure,
        receptor_chains=["A", "B"],
        ligand_chains=["D"],
        lox_domain_resi_range=(1, 200),
        cbd_domain_resi_range=None,
    )
    assert set(result) == {
        "total_ligand_contacts",
        "lox_contacts",
        "cbd_contacts",
        "lox_contact_residues",
        "cbd_contact_residues",
        "contact_cutoff_angstrom",
    }
    assert isinstance(result["total_ligand_contacts"], int)
    assert isinstance(result["lox_contacts"], int)
    assert isinstance(result["cbd_contacts"], int)
    assert isinstance(result["lox_contact_residues"], list)
    assert isinstance(result["cbd_contact_residues"], list)
    assert isinstance(result["contact_cutoff_angstrom"], float)


def test_compute_contact_frequencies_fraction() -> None:
    result = compute_contact_frequencies(
        [{"lox_contacts": 2}, {"lox_contacts": 0}, {"lox_contacts": 1}],
        keys=["lox_contacts"],
    )
    assert result["lox_contact_frequency"] == pytest.approx(0.667, abs=0.001)


def test_compute_contact_frequencies_all_zeros() -> None:
    result = compute_contact_frequencies(
        [{"lox_contacts": 0}, {"lox_contacts": 0}],
        keys=["lox_contacts"],
    )
    assert result["lox_contact_frequency"] == 0.0


def test_compute_contact_frequencies_empty() -> None:
    assert compute_contact_frequencies([], keys=["lox_contacts"])["lox_contact_frequency"] == 0.0
