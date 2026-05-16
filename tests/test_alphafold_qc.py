from __future__ import annotations

import json
from pathlib import Path

from loxdockaudit.alphafold_qc import (
    check_active_site_plddt,
    parse_pae_from_json,
    parse_plddt_from_pdb,
    summarize_alphafold_qc,
)
from loxdockaudit.pdb_parser import load_structure


def test_parse_plddt_from_pdb_extracts_ca_bfactor_only(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = parse_plddt_from_pdb(structure, chain="D", resi_list=[124, 126, 128])

    assert result["residues_requested"] == [124, 126, 128]
    assert result["residues_found"] == [124, 126, 128]
    assert result["plddt_per_residue"] == {124: 20.0, 126: 20.0, 128: 20.0}
    assert result["mean_plddt"] == 20.0
    assert result["overall_confidence"] == "low"
    assert result["n_low_confidence"] == 3


def test_parse_plddt_from_pdb_returns_required_keys(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = parse_plddt_from_pdb(structure, chain="D")

    assert set(result) == {
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
    assert isinstance(result["mean_plddt"], float)
    assert result["overall_confidence"] == "low"


def test_parse_plddt_from_pdb_no_residues_found(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = parse_plddt_from_pdb(structure, chain="D", resi_list=[999])

    assert result["residues_found"] == []
    assert result["mean_plddt"] is None
    assert result["overall_confidence"] == "unknown"
    assert result["warning"] is not None


def test_check_active_site_plddt_passes_when_all_found_and_mean_above_threshold(
    tmp_path: Path,
) -> None:
    structure = _load_plddt_fixture(tmp_path, [80.0, 90.0, 100.0])

    result = check_active_site_plddt(structure, chain="A", active_site_resi=[1, 2, 3])

    assert result["mean_plddt"] == 90.0
    assert result["all_above_threshold"] is True
    assert result["below_threshold_resi"] == []
    assert result["qc_pass"] is True
    assert result["warning"] is None


def test_check_active_site_plddt_fails_when_residue_missing(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = check_active_site_plddt(structure, chain="D", active_site_resi=[124, 999])

    assert result["qc_pass"] is False
    assert result["mean_plddt"] == 20.0
    assert result["warning"] is not None


def test_parse_pae_from_json_colabfold_format(tmp_path: Path) -> None:
    path = tmp_path / "pae_colabfold.json"
    path.write_text(json.dumps({"pae": [[1, 2, 3], [4, 5, 6], [7, 8, 9]]}), encoding="utf-8")

    result = parse_pae_from_json(str(path), domain_a_resi=[1, 2], domain_b_resi=[2, 3])

    assert result["format_detected"] == "colabfold"
    assert result["matrix_size"] == 3
    assert result["inter_domain_pae_mean"] == 4.0
    assert result["inter_domain_pae_max"] == 6.0
    assert result["inter_domain_confidence"] == "very_high"
    assert result["warning"] is None


def test_parse_pae_from_json_review_colabfold_matrix_size(tmp_path: Path) -> None:
    path = tmp_path / "pae_colabfold_review.json"
    path.write_text(json.dumps({"pae": [[5.0, 10.0], [3.0, 7.0]]}), encoding="utf-8")

    result = parse_pae_from_json(str(path), domain_a_resi=[1], domain_b_resi=[2])

    assert result["format_detected"] == "colabfold"
    assert result["matrix_size"] == 2


def test_parse_pae_from_json_alphafold_db_format(tmp_path: Path) -> None:
    path = tmp_path / "pae_afdb.json"
    path.write_text(
        json.dumps({"predicted_aligned_error": [[10, 12], [14, 16]]}),
        encoding="utf-8",
    )

    result = parse_pae_from_json(str(path), domain_a_resi=[1, 2], domain_b_resi=[1, 2])

    assert result["format_detected"] == "alphafold_db"
    assert result["inter_domain_pae_mean"] == 13.0
    assert result["inter_domain_confidence"] == "medium"


def test_parse_pae_from_json_af2_scores_format(tmp_path: Path) -> None:
    path = tmp_path / "pae_scores.json"
    path.write_text(
        json.dumps([{"predicted_aligned_error": [[25, 25], [25, 25]]}]),
        encoding="utf-8",
    )

    result = parse_pae_from_json(str(path), domain_a_resi=[1, 2], domain_b_resi=[1, 2])

    assert result["format_detected"] == "af2_scores"
    assert result["inter_domain_pae_mean"] == 25.0
    assert result["inter_domain_confidence"] == "low"


def test_parse_pae_from_json_missing_file_never_raises(tmp_path: Path) -> None:
    result = parse_pae_from_json(str(tmp_path / "missing.json"), [1], [2])

    assert result["format_detected"] == "unknown"
    assert result["matrix_size"] is None
    assert result["inter_domain_pae_mean"] is None
    assert result["warning"] is not None


def test_parse_pae_from_json_unknown_format_warns(tmp_path: Path) -> None:
    path = tmp_path / "unknown.json"
    path.write_text(json.dumps({"not_pae": [[1.0]]}), encoding="utf-8")

    result = parse_pae_from_json(str(path), [1], [1])

    assert result["format_detected"] == "unknown"
    assert result["inter_domain_pae_mean"] is None
    assert result["warning"] is not None


def test_check_active_site_plddt_returns_required_keys(fixture_pdb_path: Path) -> None:
    structure = load_structure(str(fixture_pdb_path))

    result = check_active_site_plddt(
        structure,
        chain="D",
        active_site_resi=[124, 126, 128, 152, 187],
    )

    assert set(result) == {
        "active_site_resi",
        "plddt_per_residue",
        "mean_plddt",
        "all_above_threshold",
        "below_threshold_resi",
        "threshold_used",
        "qc_pass",
        "warning",
    }
    assert isinstance(result["below_threshold_resi"], list)


def test_summarize_alphafold_qc_passes_without_pae() -> None:
    result = summarize_alphafold_qc(
        {"overall_confidence": "high"},
        {"qc_pass": True, "mean_plddt": 82.0},
    )

    assert result == {
        "overall_plddt_confidence": "high",
        "active_site_plddt_pass": True,
        "active_site_mean_plddt": 82.0,
        "inter_domain_pae_confidence": None,
        "fold_qc_pass": True,
        "fold_qc_warning": None,
    }


def test_summarize_alphafold_qc_fails_on_low_pae() -> None:
    result = summarize_alphafold_qc(
        {"overall_confidence": "high"},
        {"qc_pass": True, "mean_plddt": 82.0},
        {"inter_domain_confidence": "low"},
    )

    assert result["fold_qc_pass"] is False
    assert result["fold_qc_warning"] is not None


def _load_plddt_fixture(tmp_path: Path, scores: list[float]):
    lines = []
    for index, score in enumerate(scores, start=1):
        lines.append(
            f"ATOM  {index:5d}  CA  ALA A{index:4d}    "
            f"{float(index):8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{score:6.2f}           C  "
        )
    lines.extend(["TER", "END"])
    pdb_path = tmp_path / "plddt.pdb"
    pdb_path.write_text("\n".join(lines), encoding="utf-8")
    return load_structure(str(pdb_path))
