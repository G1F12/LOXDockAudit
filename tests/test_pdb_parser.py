from __future__ import annotations

import re
from pathlib import Path

import pytest
from Bio.PDB.Structure import Structure

from loxdockaudit.pdb_parser import (
    compute_sha256,
    get_atoms_by_selection,
    get_chain_ids,
    get_sequence_from_chain,
    load_structure,
)


def test_load_structure_returns_structure(fixture_pdb_path: Path) -> None:
    assert isinstance(load_structure(str(fixture_pdb_path)), Structure)


def test_get_chain_ids(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    assert get_chain_ids(structure) == ["A", "B", "D"]


def test_get_atoms_by_selection_residue_list(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    atoms = get_atoms_by_selection(structure, chains=["D"], resi_list=[124])
    assert atoms
    assert {atom.get_parent().resname for atom in atoms} == {"HIS"}
    assert {atom.get_parent().id[1] for atom in atoms} == {124}


def test_get_atoms_by_selection_residue_range(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    atoms = get_atoms_by_selection(structure, chains=["A"], resi_range=(1, 5))
    assert atoms
    assert {atom.get_parent().id[1] for atom in atoms} == {1, 2, 3, 4, 5}


def test_get_atoms_by_selection_excludes_hydrogens(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    atoms = get_atoms_by_selection(structure, chains=["D"], resi_list=[124])
    assert all(atom.element.strip().upper() != "H" for atom in atoms)


def test_get_sequence_from_chain(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    assert get_sequence_from_chain(structure, "A") == "AAAAAAAAAAK"


def test_get_sequence_from_missing_chain(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    assert get_sequence_from_chain(structure, "Z") == ""


def test_compute_sha256_is_hex_string(fixture_pdb_path: Path) -> None:
    digest = compute_sha256(str(fixture_pdb_path))
    assert len(digest) == 64
    assert re.fullmatch(r"[0-9a-f]{64}", digest)


def test_compute_sha256_differs_for_different_files(
    fixture_pdb_path: Path,
    tmp_path: Path,
) -> None:
    other = tmp_path / "other.pdb"
    other.write_text("END\n", encoding="utf-8")
    assert compute_sha256(str(fixture_pdb_path)) != compute_sha256(str(other))


def test_load_structure_nonexistent_path_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_structure(str(tmp_path / "missing.pdb"))
