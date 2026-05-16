from __future__ import annotations

import shutil
from pathlib import Path

from loxdockaudit.input_audit import audit_screen, audit_single_file, format_audit_report


def test_audit_single_file_valid_fixture(fixture_pdb_path: Path) -> None:
    result = audit_single_file(str(fixture_pdb_path), expected_chains=["A", "B", "D"])
    assert result["valid"] is True


def test_audit_single_file_nonexistent(tmp_path: Path) -> None:
    result = audit_single_file(str(tmp_path / "missing.pdb"), expected_chains=["A"])
    assert result["exists"] is False
    assert result["valid"] is False


def test_audit_single_file_empty(tmp_path: Path) -> None:
    pdb_path = tmp_path / "empty.pdb"
    pdb_path.write_bytes(b"")
    result = audit_single_file(str(pdb_path), expected_chains=["A"])
    assert result["parseable"] is False
    assert result["valid"] is False


def test_audit_single_file_correct_n_terminus(fixture_pdb_path: Path) -> None:
    result = audit_single_file(
        str(fixture_pdb_path),
        expected_chains=["A", "B", "D"],
        expected_n_terminus="AAAAAAAAAA",
    )
    assert result["n_terminus_ok"] is True


def test_audit_single_file_wrong_n_terminus(fixture_pdb_path: Path) -> None:
    result = audit_single_file(
        str(fixture_pdb_path),
        expected_chains=["A", "B", "D"],
        expected_n_terminus="KKKKKKKKKK",
    )
    assert result["n_terminus_ok"] is False
    assert result["valid"] is False


def test_audit_screen_different_files_no_duplicates(
    fixture_pdb_path: Path,
    tmp_path: Path,
) -> None:
    c1 = tmp_path / "c1"
    c2 = tmp_path / "c2"
    c1.mkdir()
    c2.mkdir()
    shutil.copyfile(fixture_pdb_path, c1 / "model_1.pdb")
    (c2 / "model_1.pdb").write_text("END\n", encoding="utf-8")
    result = audit_screen(
        [
            _audit_config("c1", c1),
            {"construct_id": "c2", "models_dir": str(c2), "ligand_chains": [], "receptor_chains": []},
        ]
    )
    assert result["constructs"]["c1"]["duplicate_sha256"] == []
    assert result["screen_valid"] is True


def test_audit_screen_cross_construct_duplicate(
    fixture_pdb_path: Path,
    tmp_path: Path,
) -> None:
    c1 = tmp_path / "candidate"
    c2 = tmp_path / "control"
    c1.mkdir()
    c2.mkdir()
    shutil.copyfile(fixture_pdb_path, c1 / "candidate.pdb")
    shutil.copyfile(fixture_pdb_path, c2 / "control.pdb")
    result = audit_screen([_audit_config("candidate", c1), _audit_config("control", c2)])
    assert result["screen_valid"] is False
    assert result["constructs"]["candidate"]["cross_construct_duplicates"]


def test_format_audit_report_contains_markdown_headers(fixture_pdb_path: Path, tmp_path: Path) -> None:
    models = tmp_path / "models"
    models.mkdir()
    shutil.copyfile(fixture_pdb_path, models / "model_1.pdb")
    report = format_audit_report(audit_screen([_audit_config("c1", models)]))
    assert "##" in report


def test_format_audit_report_invalid_contains_fail(tmp_path: Path) -> None:
    result = audit_screen([_audit_config("missing", tmp_path / "missing")])
    report = format_audit_report(result)
    assert "ERROR" in report.upper() or "FAIL" in report


def _audit_config(construct_id: str, models_dir: Path) -> dict:
    return {
        "construct_id": construct_id,
        "models_dir": str(models_dir),
        "ligand_chains": ["D"],
        "receptor_chains": ["A", "B"],
        "expected_n_terminus": None,
    }
