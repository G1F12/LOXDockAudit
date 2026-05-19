from __future__ import annotations

import numpy as np
import pytest
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from loxdockaudit.distances import (
    OrientationConfig,
    PoseGeometry,
    analyze_pose_geometry,
)


def test_pose_geometry_orientation_disabled() -> None:
    result = analyze_pose_geometry(
        structure=_make_pose(cb=(1.0, 0.0, 0.0), nz=(2.0, 0.0, 0.0)),
        pose_rank=1,
        active_site_resi=[100],
        active_site_chains=["D"],
        target_resi=11,
        target_chain="A",
        orientation_config=OrientationConfig(enabled=False),
        active_site_residues=[{"label": "His292", "resi": 100, "chain": "D"}],
    )

    assert result.orientation_angle_deg is None
    assert result.orientation_productive is None
    assert result.fully_productive == result.distance_productive


def test_pose_geometry_orientation_enabled_productive() -> None:
    result = analyze_pose_geometry(
        structure=_make_pose(cb=(1.0, 0.0, 0.0), nz=(2.0, 0.0, 0.0)),
        pose_rank=1,
        active_site_resi=[100],
        active_site_chains=["D"],
        target_resi=11,
        target_chain="A",
        threshold=8.0,
        target_atom_names=["NZ"],
        orientation_config=OrientationConfig(enabled=True, threshold_deg=90.0),
        active_site_residues=[{"label": "His292", "resi": 100, "chain": "D"}],
    )

    assert result.distance_productive is True
    assert result.orientation_angle_deg == pytest.approx(0.0, abs=1e-6)
    assert result.orientation_productive is True
    assert result.fully_productive is True


def test_pose_geometry_orientation_enabled_antiproductive() -> None:
    result = analyze_pose_geometry(
        structure=_make_pose(cb=(3.0, 0.0, 0.0), nz=(2.0, 0.0, 0.0)),
        pose_rank=1,
        active_site_resi=[100],
        active_site_chains=["D"],
        target_resi=11,
        target_chain="A",
        threshold=8.0,
        target_atom_names=["NZ"],
        orientation_config=OrientationConfig(enabled=True, threshold_deg=90.0),
        active_site_residues=[{"label": "His292", "resi": 100, "chain": "D"}],
    )

    assert result.distance_productive is True
    assert result.orientation_angle_deg == pytest.approx(180.0, abs=1e-6)
    assert result.orientation_productive is False
    assert result.fully_productive is False


def test_fully_productive_requires_both() -> None:
    assert _pose_geometry(True, False).fully_productive is False
    assert _pose_geometry(False, True).fully_productive is False
    assert _pose_geometry(True, True).fully_productive is True


def test_backward_compat_no_orientation_config() -> None:
    result = analyze_pose_geometry(
        structure=_make_pose(cb=(1.0, 0.0, 0.0), nz=(2.0, 0.0, 0.0)),
        pose_rank=7,
        active_site_resi=[100],
        active_site_chains=["D"],
        target_resi=11,
        target_chain="A",
        threshold=8.0,
        target_atom_names=["NZ"],
        active_site_residues=[{"label": "His292", "resi": 100, "chain": "D"}],
    )

    assert result.pose_rank == 7
    assert result.distance_A == pytest.approx(2.0, abs=1e-6)
    assert result.distance_productive is True
    assert result.orientation_angle_deg is None
    assert result.orientation_productive is None
    assert result.fully_productive == result.distance_productive


def _pose_geometry(
    distance_productive: bool,
    orientation_productive: bool,
) -> PoseGeometry:
    return PoseGeometry(
        pose_rank=1,
        distance_A=1.0 if distance_productive else 10.0,
        distance_productive=distance_productive,
        orientation_angle_deg=45.0 if orientation_productive else 135.0,
        orientation_productive=orientation_productive,
        fully_productive=distance_productive and orientation_productive,
    )


def _make_pose(
    cb: tuple[float, float, float],
    nz: tuple[float, float, float],
) -> Structure:
    structure = Structure("pose")
    model = Model(0)
    structure.add(model)

    active_chain = Chain("D")
    active_residue = Residue((" ", 100, " "), "HIS", "")
    active_residue.add(_atom("CA", (0.0, 0.0, 0.0), element="C"))
    active_chain.add(active_residue)
    model.add(active_chain)

    target_chain = Chain("A")
    target_residue = Residue((" ", 11, " "), "LYS", "")
    target_residue.add(_atom("CB", cb, element="C"))
    target_residue.add(_atom("NZ", nz, element="N"))
    target_chain.add(target_residue)
    model.add(target_chain)

    return structure


def _atom(
    name: str,
    coord: tuple[float, float, float],
    element: str,
) -> Atom:
    return Atom(
        name=name,
        coord=np.array(coord, dtype=float),
        bfactor=0.0,
        occupancy=1.0,
        altloc=" ",
        fullname=f"{name:<4}",
        serial_number=1,
        element=element,
    )
