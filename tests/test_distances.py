from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from loxdockaudit.distances import (
    active_site_to_target_distance,
    min_heavy_atom_distance,
    score_all_poses,
)
from loxdockaudit.pdb_parser import load_structure


@dataclass
class Residue:
    id: tuple[str, int, str]


@dataclass
class Atom:
    coord: np.ndarray
    element: str = "C"
    name: str = "CA"

    def get_parent(self) -> Residue:
        return Residue((" ", 1, " "))


def atom(x: float, y: float = 0.0, z: float = 0.0) -> Atom:
    return Atom(np.array([x, y, z], dtype=float))


def test_min_heavy_atom_distance_5_angstrom() -> None:
    assert min_heavy_atom_distance([atom(0)], [atom(5)]) == pytest.approx(5.0)


def test_min_heavy_atom_distance_empty_a() -> None:
    assert np.isinf(min_heavy_atom_distance([], [atom(5)]))


def test_min_heavy_atom_distance_empty_b() -> None:
    assert np.isinf(min_heavy_atom_distance([atom(0)], []))


def test_min_heavy_atom_distance_multiple_atoms_closest_3() -> None:
    assert min_heavy_atom_distance([atom(0), atom(20)], [atom(10), atom(23)]) == pytest.approx(3.0)


def test_active_site_to_target_distance_fixture(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = active_site_to_target_distance(
        structure,
        active_site_resi=[124, 126, 128, 152, 187],
        active_site_chains=["D"],
        target_resi=11,
        target_chain="A",
    )
    assert np.isfinite(result["distance_angstrom"])
    assert result["productive"] is True


def test_active_site_to_target_distance_missing_target(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = active_site_to_target_distance(
        structure,
        active_site_resi=[124],
        active_site_chains=["D"],
        target_resi=999,
        target_chain="A",
    )
    assert result["productive"] is False
    assert result["warning"] is not None


def test_score_all_poses_repeated_fixture(fixture_pdb_path: Path) -> None:
    scores = score_all_poses(
        [str(fixture_pdb_path), str(fixture_pdb_path), str(fixture_pdb_path)],
        active_site_resi=[124, 126, 128, 152, 187],
        active_site_chains=["D"],
        target_resi=11,
        target_chain="A",
    )
    assert len(scores) == 3
    assert [score["rank"] for score in scores] == [1, 2, 3]
    assert len({score["distance_angstrom"] for score in scores}) == 1
