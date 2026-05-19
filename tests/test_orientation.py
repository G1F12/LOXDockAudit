from __future__ import annotations

import numpy as np
import pytest

from loxdockaudit.orientation import (
    compute_approach_vector,
    compute_lys_sidechain_vector,
    compute_orientation_angle,
    get_atom_coords,
    is_orientation_productive,
)


def test_parallel_vectors_zero_angle() -> None:
    sidechain_vec = np.array([1.0, 0.0, 0.0])
    approach_vec = np.array([1.0, 0.0, 0.0])

    angle = compute_orientation_angle(sidechain_vec, approach_vec)

    assert angle == pytest.approx(0.0, abs=1e-6)


def test_antiparallel_vectors_180_angle() -> None:
    sidechain_vec = np.array([1.0, 0.0, 0.0])
    approach_vec = np.array([-1.0, 0.0, 0.0])

    angle = compute_orientation_angle(sidechain_vec, approach_vec)

    assert angle == pytest.approx(180.0, abs=1e-6)


def test_perpendicular_vectors_90_angle() -> None:
    sidechain_vec = np.array([1.0, 0.0, 0.0])
    approach_vec = np.array([0.0, 1.0, 0.0])

    angle = compute_orientation_angle(sidechain_vec, approach_vec)

    assert angle == pytest.approx(90.0, abs=1e-6)


def test_is_orientation_productive_boundary() -> None:
    assert is_orientation_productive(89.9) is True
    assert is_orientation_productive(90.0) is True
    assert is_orientation_productive(90.1) is False


def test_is_orientation_productive_custom_threshold() -> None:
    assert is_orientation_productive(74.0, threshold_deg=75.0) is True
    assert is_orientation_productive(76.0, threshold_deg=75.0) is False


def test_compute_orientation_angle_clamps_dot_product() -> None:
    sidechain_vec = np.array([1.0, 0.0, 0.0])
    approach_vec = np.array([1.0 + 1e-12, 0.0, 0.0])

    angle = compute_orientation_angle(sidechain_vec, approach_vec)

    assert isinstance(angle, float)
    assert np.isfinite(angle)
    assert 0.0 <= angle <= 180.0


def test_get_atom_coords_found() -> None:
    atoms = [
        {
            "atom_name": "CA",
            "x": 1.0,
            "y": 2.0,
            "z": 3.0,
        }
    ]

    coords = get_atom_coords(atoms, "CA")

    np.testing.assert_allclose(coords, np.array([1.0, 2.0, 3.0]))


def test_get_atom_coords_not_found() -> None:
    atoms = [
        {
            "atom_name": "CB",
            "x": 1.0,
            "y": 2.0,
            "z": 3.0,
        }
    ]

    with pytest.raises(ValueError):
        get_atom_coords(atoms, "CA")


def test_compute_approach_vector_zero_length() -> None:
    centroid = np.array([1.0, 2.0, 3.0])
    nz_coords = np.array([1.0, 2.0, 3.0])

    with pytest.raises(ValueError):
        compute_approach_vector(centroid, nz_coords)


def test_compute_lys_sidechain_vector_missing_atom() -> None:
    pose_atoms = [
        {
            "atom_name": "CB",
            "chain_id": "A",
            "resseq": 11,
            "x": 1.0,
            "y": 2.0,
            "z": 3.0,
        }
    ]

    with pytest.raises(ValueError):
        compute_lys_sidechain_vector(
            pose_atoms,
            target_resi=11,
            target_chain="A",
        )
