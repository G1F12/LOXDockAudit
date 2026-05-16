from __future__ import annotations

import math
from pathlib import Path

from loxdockaudit.pdb_parser import load_structure
from loxdockaudit.structural_qc import (
    check_active_site_sasa_proxy,
    check_disulfide_geometry,
    check_his_triad_geometry,
    check_lys_tyr_geometry,
)


def test_check_his_triad_geometry_plausible(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = check_his_triad_geometry(structure, [124, 126, 128], chain="D")

    assert result["ca_coords_found"] == [124, 126, 128]
    assert result["ca_coords_missing"] == []
    assert result["geometry_plausible"] is True
    assert result["warning"] is None
    assert result["max_ca_distance"] < 10.0


def test_check_his_triad_geometry_chain_d_all_pairwise_distances_are_finite(
    fixture_pdb_path: Path,
) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = check_his_triad_geometry(structure, [124, 126, 128], chain="D")

    assert len(result["ca_coords_found"]) == 3
    assert set(result["pairwise_ca_distances"]) == {"124-126", "124-128", "126-128"}
    assert all(math.isfinite(distance) for distance in result["pairwise_ca_distances"].values())


def test_check_his_triad_geometry_missing_ca(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = check_his_triad_geometry(structure, [124, 126, 999], chain="D")

    assert result["ca_coords_found"] == [124, 126]
    assert result["ca_coords_missing"] == [999]
    assert result["geometry_plausible"] is False
    assert result["warning"] is not None


def test_check_his_triad_geometry_all_missing(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = check_his_triad_geometry(structure, [999, 998, 997], chain="D")

    assert result["geometry_plausible"] is False
    assert result["ca_coords_missing"] == [999, 998, 997]


def test_check_lys_tyr_geometry_uses_ca_fallback(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = check_lys_tyr_geometry(structure, lys_resi=152, tyr_resi=187, chain="D")

    assert result["lys_cb_found"] is False
    assert result["tyr_cb_found"] is False
    assert result["cb_distance"] == 2.0
    assert result["geometry_plausible"] is True
    assert result["warning"] is not None


def test_check_lys_tyr_geometry_returns_required_keys(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = check_lys_tyr_geometry(structure, lys_resi=152, tyr_resi=187, chain="D")

    assert set(result) == {
        "lys_resi",
        "tyr_resi",
        "lys_cb_found",
        "tyr_cb_found",
        "cb_distance",
        "geometry_plausible",
        "warning",
    }


def test_check_disulfide_geometry(tmp_path: Path) -> None:
    pdb_path = tmp_path / "disulfide.pdb"
    pdb_path.write_text(
        "\n".join(
            [
                "ATOM      1  SG  CYS A   1       0.000   0.000   0.000  1.00 20.00           S  ",
                "ATOM      2  SG  CYS A   2       2.050   0.000   0.000  1.00 20.00           S  ",
                "ATOM      3  SG  CYS A   3      10.000   0.000   0.000  1.00 20.00           S  ",
                "ATOM      4  SG  CYS A   4      14.000   0.000   0.000  1.00 20.00           S  ",
                "TER",
                "END",
            ]
        ),
        encoding="utf-8",
    )
    structure = load_structure(str(pdb_path))

    result = check_disulfide_geometry(structure, [(1, 2), (3, 4)], chain="A")

    assert result["pairs_checked"] == 2
    assert result["pairs_intact"] == 1
    assert result["all_intact"] is False
    assert result["pair_results"][0]["intact"] is True
    assert result["pair_results"][1]["intact"] is False
    assert result["warning"] is not None


def test_check_disulfide_geometry_empty_pairs(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = check_disulfide_geometry(structure, [], chain="D")

    assert result["pairs_checked"] == 0
    assert result["all_intact"] is True
    assert result["pair_results"] == []


def test_check_active_site_sasa_proxy(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = check_active_site_sasa_proxy(structure, [124, 126, 128], chain="D")

    assert result["n_active_site_atoms"] == 4
    assert result["n_surrounding_atoms_within_cutoff"] == 4
    assert result["burial_score"] == 1.0
    assert result["accessibility_flag"] == "accessible"
    assert result["warning"] is None


def test_check_active_site_sasa_proxy_missing_active_site(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = check_active_site_sasa_proxy(structure, [999], chain="D")

    assert result["n_active_site_atoms"] == 0
    assert result["warning"] is not None
