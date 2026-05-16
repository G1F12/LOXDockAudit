from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture
def fixture_pdb_path() -> Path:
    return Path(__file__).parent / "fixtures" / "test_model.pdb"


def pdb_text_for_distance(distance: float) -> str:
    return (
        f"ATOM      1  NZ  LYS A  11       0.000   0.000   0.000  1.00 20.00           N  \n"
        f"ATOM      2  CA  GLY B   1      30.000   0.000   0.000  1.00 20.00           C  \n"
        f"ATOM      3  CA  HIS D 124      {distance:6.3f}   0.000   0.000  1.00 20.00           C  \n"
        f"ATOM      4  CA  HIS D 126      {distance:6.3f}   1.000   0.000  1.00 20.00           C  \n"
        f"ATOM      5  CA  HIS D 128      {distance:6.3f}  -1.000   0.000  1.00 20.00           C  \n"
        f"ATOM      6  CA  LYS D 152      {distance:6.3f}   0.000   1.000  1.00 20.00           C  \n"
        f"ATOM      7  CA  TYR D 187      {distance:6.3f}   0.000  -1.000  1.00 20.00           C  \n"
        "TER\n"
        "END\n"
    )


def write_pose(path: Path, distance: float) -> Path:
    path.write_text(pdb_text_for_distance(distance), encoding="utf-8")
    return path
