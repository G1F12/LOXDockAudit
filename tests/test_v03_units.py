from __future__ import annotations

import json
from pathlib import Path

from loxdockaudit.alphafold_qc import (
    check_active_site_plddt,
    parse_pae_from_json,
    parse_plddt_from_pdb,
    summarize_alphafold_qc,
)
from loxdockaudit.fold_qc import StructuralQCResult, format_qc_report, run_structural_qc
from loxdockaudit.pdb_parser import load_structure
from loxdockaudit.structural_qc import (
    check_active_site_sasa_proxy,
    check_disulfide_geometry,
    check_his_triad_geometry,
    check_lys_tyr_geometry,
)


REQUIRED_PLDDT_KEYS = {
    "chain",
    "residues_requested",
    "residues_found",
    "plddt_per_residue",
    "mean_plddt",
    "min_plddt",
    "max_plddt",
    "n_low_confidence",
    "n_medium_confidence",
    "n_high_confidence",
    "n_very_high_confidence",
    "overall_confidence",
    "warning",
}


def test_v03_01_his_triad_three_found_plausible(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = check_his_triad_geometry(structure, [124, 126, 128], chain="D")

    assert result["geometry_plausible"] is True


def test_v03_02_his_triad_one_missing_not_plausible(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = check_his_triad_geometry(structure, [124, 126, 999], chain="D")

    assert result["geometry_plausible"] is False


def test_v03_03_his_triad_all_missing_all_keys_present(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = check_his_triad_geometry(structure, [999, 998, 997], chain="D")

    assert set(result) == {
        "his_resi",
        "ca_coords_found",
        "ca_coords_missing",
        "pairwise_ca_distances",
        "max_ca_distance",
        "geometry_plausible",
        "warning",
    }
    assert result["ca_coords_missing"] == [999, 998, 997]


def test_v03_04_lys_tyr_both_found_distance_float(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = check_lys_tyr_geometry(structure, 152, 187, chain="D")

    assert isinstance(result["cb_distance"], float)


def test_v03_05_lys_tyr_both_missing_not_plausible(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = check_lys_tyr_geometry(structure, 999, 998, chain="D")

    assert result["geometry_plausible"] is False


def test_v03_06_disulfide_empty_pairs_intact(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = check_disulfide_geometry(structure, [], chain="D")

    assert result["all_intact"] is True
    assert result["pairs_checked"] == 0


def test_v03_07_disulfide_sg_missing_pair_not_intact(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = check_disulfide_geometry(structure, [(124, 126)], chain="D")

    assert result["pair_results"][0]["intact"] is False


def test_v03_08_sasa_proxy_valid_resi_burial_score_float(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = check_active_site_sasa_proxy(structure, [124, 126, 128], chain="D")

    assert isinstance(result["burial_score"], float)


def test_v03_09_sasa_proxy_empty_resi_warning(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = check_active_site_sasa_proxy(structure, [], chain="D")

    assert result["warning"] is not None


def test_v03_10_sasa_proxy_accessibility_flag_allowed(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = check_active_site_sasa_proxy(structure, [124, 126, 128], chain="D")

    assert result["accessibility_flag"] in ["accessible", "partially_buried", "buried"]


def test_v03_11_parse_plddt_returns_required_keys(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = parse_plddt_from_pdb(structure, chain="D")

    assert set(result) == REQUIRED_PLDDT_KEYS


def test_v03_12_plddt_bucket_counts_sum_to_total(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = parse_plddt_from_pdb(structure, chain="D")
    bucket_total = (
        result["n_low_confidence"]
        + result["n_medium_confidence"]
        + result["n_high_confidence"]
        + result["n_very_high_confidence"]
    )

    assert bucket_total == len(result["residues_found"])


def test_v03_13_parse_pae_colabfold_format(tmp_path: Path) -> None:
    path = tmp_path / "pae.json"
    path.write_text(json.dumps({"pae": [[5.0, 10.0], [3.0, 7.0]]}), encoding="utf-8")

    result = parse_pae_from_json(str(path), [1], [2])

    assert result["format_detected"] == "colabfold"


def test_v03_14_parse_pae_missing_file_warns(tmp_path: Path) -> None:
    result = parse_pae_from_json(str(tmp_path / "missing.json"), [1], [2])

    assert result["warning"] is not None


def test_v03_15_parse_pae_unknown_format(tmp_path: Path) -> None:
    path = tmp_path / "unknown.json"
    path.write_text(json.dumps({"not_pae": []}), encoding="utf-8")

    result = parse_pae_from_json(str(path), [1], [2])

    assert result["format_detected"] == "unknown"


def test_v03_16_active_site_plddt_qc_pass_bool(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))
    result = check_active_site_plddt(structure, chain="D", active_site_resi=[124, 126, 128])

    assert isinstance(result["qc_pass"], bool)


def test_v03_17_summary_without_pae_uses_plddt_only() -> None:
    result = summarize_alphafold_qc(
        {"overall_confidence": "high"},
        {"qc_pass": True, "mean_plddt": 80.0},
        pae_result=None,
    )

    assert result["fold_qc_pass"] is True


def test_v03_18_summary_active_site_fail_fails_fold_qc() -> None:
    result = summarize_alphafold_qc(
        {"overall_confidence": "high"},
        {"qc_pass": False, "mean_plddt": 40.0},
        pae_result=None,
    )

    assert result["fold_qc_pass"] is False


def test_v03_19_run_structural_qc_returns_result(fixture_pdb_path: Path) -> None:
    result = _run_fixture_qc(fixture_pdb_path)

    assert isinstance(result, StructuralQCResult)


def test_v03_20_run_structural_qc_empty_disulfides_pass(fixture_pdb_path: Path) -> None:
    result = _run_fixture_qc(fixture_pdb_path)

    assert result.disulfide_pass is True


def test_v03_21_run_structural_qc_no_alphafold_plddt_none(fixture_pdb_path: Path) -> None:
    result = _run_fixture_qc(fixture_pdb_path)

    assert result.plddt_pass is None


def test_v03_22_run_structural_qc_fold_corrupted_inverse(fixture_pdb_path: Path) -> None:
    result = _run_fixture_qc(fixture_pdb_path)

    assert result.fold_corrupted is (not result.fold_qc_pass)


def test_v03_23_format_qc_report_contains_header(fixture_pdb_path: Path) -> None:
    report = format_qc_report(_run_fixture_qc(fixture_pdb_path))

    assert "## Structural QC" in report


def test_v03_24_format_qc_report_corrupted_label(fixture_pdb_path: Path) -> None:
    result = run_structural_qc(
        pdb_path=str(fixture_pdb_path),
        construct_id="fixture",
        lox_chain="D",
        his_resi=[999, 998, 997],
        lys_resi=152,
        tyr_resi=187,
        disulfide_pairs=[],
        active_site_resi=[124, 126, 128, 152, 187],
    )

    assert "FOLD CORRUPTED" in format_qc_report(result)


def test_v03_25_format_qc_report_empty_disulfides_skipped(fixture_pdb_path: Path) -> None:
    report = format_qc_report(_run_fixture_qc(fixture_pdb_path))

    assert "SKIPPED" in report


def _run_fixture_qc(fixture_pdb_path: Path) -> StructuralQCResult:
    return run_structural_qc(
        pdb_path=str(fixture_pdb_path),
        construct_id="fixture",
        lox_chain="D",
        his_resi=[124, 126, 128],
        lys_resi=152,
        tyr_resi=187,
        disulfide_pairs=[],
        active_site_resi=[124, 126, 128, 152, 187],
    )
