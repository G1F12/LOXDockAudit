from __future__ import annotations

import json
from dataclasses import fields
from pathlib import Path

from loxdockaudit.fold_qc import (
    StructuralQCResult,
    format_qc_report,
    run_structural_qc,
)


def test_run_structural_qc_without_alphafold_passes_fixture(fixture_pdb_path: Path) -> None:
    result = run_structural_qc(
        pdb_path=str(fixture_pdb_path),
        construct_id="construct",
        lox_chain="D",
        his_resi=[124, 126, 128],
        lys_resi=152,
        tyr_resi=187,
        disulfide_pairs=[],
        active_site_resi=[124, 126, 128, 152, 187],
    )

    assert isinstance(result, StructuralQCResult)
    assert result.his_triad_pass is True
    assert result.lys_tyr_pass is True
    assert result.disulfide_pass is True
    assert result.active_site_accessible is True
    assert result.plddt_pass is None
    assert result.pae_pass is None
    assert result.fold_qc_pass is True
    assert result.fold_corrupted is False
    assert result.details["disulfide"]["skipped"] is True


def test_run_structural_qc_returns_result_with_expected_field_types(
    fixture_pdb_path: Path,
) -> None:
    result = run_structural_qc(
        pdb_path=str(fixture_pdb_path),
        construct_id="construct",
        lox_chain="D",
        his_resi=[124, 126, 128],
        lys_resi=152,
        tyr_resi=187,
        disulfide_pairs=[],
        active_site_resi=[124, 126, 128, 152, 187],
    )

    assert [field.name for field in fields(StructuralQCResult)] == [
        "construct_id",
        "his_triad_pass",
        "lys_tyr_pass",
        "disulfide_pass",
        "active_site_accessible",
        "plddt_pass",
        "pae_pass",
        "fold_qc_pass",
        "fold_corrupted",
        "warnings",
        "details",
    ]
    assert isinstance(result.construct_id, str)
    assert isinstance(result.his_triad_pass, bool)
    assert isinstance(result.lys_tyr_pass, bool)
    assert isinstance(result.disulfide_pass, bool)
    assert isinstance(result.active_site_accessible, bool)
    assert result.plddt_pass is None or isinstance(result.plddt_pass, bool)
    assert result.pae_pass is None or isinstance(result.pae_pass, bool)
    assert isinstance(result.fold_qc_pass, bool)
    assert result.fold_corrupted is (not result.fold_qc_pass)
    assert isinstance(result.warnings, list)
    assert all(isinstance(warning, str) for warning in result.warnings)
    assert isinstance(result.details, dict)


def test_run_structural_qc_empty_disulfide_pairs_are_skipped(
    fixture_pdb_path: Path,
) -> None:
    result = run_structural_qc(
        pdb_path=str(fixture_pdb_path),
        construct_id="construct",
        lox_chain="D",
        his_resi=[124, 126, 128],
        lys_resi=152,
        tyr_resi=187,
        disulfide_pairs=[],
        active_site_resi=[124, 126, 128, 152, 187],
    )

    assert result.disulfide_pass is True
    assert result.details["disulfide"]["skipped"] is True


def test_run_structural_qc_without_alphafold_skips_af_checks(
    fixture_pdb_path: Path,
) -> None:
    result = run_structural_qc(
        pdb_path=str(fixture_pdb_path),
        construct_id="construct",
        lox_chain="D",
        his_resi=[124, 126, 128],
        lys_resi=152,
        tyr_resi=187,
        disulfide_pairs=[],
        active_site_resi=[124, 126, 128, 152, 187],
        alphafold_pdb_path=None,
    )

    assert result.plddt_pass is None
    assert result.pae_pass is None
    assert result.fold_qc_pass is (
        result.his_triad_pass
        and result.lys_tyr_pass
        and result.disulfide_pass
        and result.active_site_accessible
    )


def test_run_structural_qc_with_alphafold_and_pae_passes_when_threshold_low(
    fixture_pdb_path: Path,
    tmp_path: Path,
) -> None:
    pae_path = tmp_path / "pae.json"
    pae_path.write_text(json.dumps({"pae": [[5.0, 10.0], [3.0, 7.0]]}), encoding="utf-8")

    result = run_structural_qc(
        pdb_path=str(fixture_pdb_path),
        construct_id="construct",
        lox_chain="D",
        his_resi=[124, 126, 128],
        lys_resi=152,
        tyr_resi=187,
        disulfide_pairs=[],
        active_site_resi=[124, 126, 128],
        alphafold_pdb_path=str(fixture_pdb_path),
        pae_json_path=str(pae_path),
        domain_a_resi=[1],
        domain_b_resi=[2],
        plddt_threshold=10.0,
    )

    assert result.plddt_pass is True
    assert result.pae_pass is True
    assert result.fold_qc_pass is True
    assert result.details["pae"]["format_detected"] == "colabfold"


def test_run_structural_qc_with_alphafold_fails_on_low_plddt(fixture_pdb_path: Path) -> None:
    result = run_structural_qc(
        pdb_path=str(fixture_pdb_path),
        construct_id="construct",
        lox_chain="D",
        his_resi=[124, 126, 128],
        lys_resi=152,
        tyr_resi=187,
        disulfide_pairs=[],
        active_site_resi=[124, 126, 128],
        alphafold_pdb_path=str(fixture_pdb_path),
        plddt_threshold=70.0,
    )

    assert result.plddt_pass is False
    assert result.fold_qc_pass is False
    assert result.fold_corrupted is True
    assert result.warnings


def test_run_structural_qc_missing_pdb_returns_corrupted(tmp_path: Path) -> None:
    result = run_structural_qc(
        pdb_path=str(tmp_path / "missing.pdb"),
        construct_id="missing",
        lox_chain="D",
        his_resi=[124, 126, 128],
        lys_resi=152,
        tyr_resi=187,
        disulfide_pairs=[],
        active_site_resi=[124, 126, 128],
    )

    assert result.fold_qc_pass is False
    assert result.fold_corrupted is True
    assert result.warnings


def test_format_qc_report_contains_required_sections(fixture_pdb_path: Path) -> None:
    result = run_structural_qc(
        pdb_path=str(fixture_pdb_path),
        construct_id="construct",
        lox_chain="D",
        his_resi=[124, 126, 128],
        lys_resi=152,
        tyr_resi=187,
        disulfide_pairs=[],
        active_site_resi=[124, 126, 128, 152, 187],
    )

    report = format_qc_report(result)

    assert "## Structural QC: construct" in report
    assert "**Overall: PASS**" in report
    assert "| His triad geometry | PASS |" in report
    assert "| Disulfide geometry | SKIPPED |" in report
    assert "| pLDDT (active site) | N/A |" in report
    assert "### Warnings" in report


def test_format_qc_report_returns_markdown_table(fixture_pdb_path: Path) -> None:
    result = run_structural_qc(
        pdb_path=str(fixture_pdb_path),
        construct_id="construct",
        lox_chain="D",
        his_resi=[124, 126, 128],
        lys_resi=152,
        tyr_resi=187,
        disulfide_pairs=[],
        active_site_resi=[124, 126, 128, 152, 187],
    )

    report = format_qc_report(result)

    assert isinstance(report, str)
    assert "## Structural QC" in report
    assert "|" in report


def test_format_qc_report_fold_corrupted_label(fixture_pdb_path: Path) -> None:
    result = run_structural_qc(
        pdb_path=str(fixture_pdb_path),
        construct_id="construct",
        lox_chain="D",
        his_resi=[124, 126, 999],
        lys_resi=152,
        tyr_resi=187,
        disulfide_pairs=[],
        active_site_resi=[124, 126, 128, 152, 187],
    )

    report = format_qc_report(result)

    assert result.fold_qc_pass is False
    assert "FOLD CORRUPTED" in report
