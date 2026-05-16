"""Aggregate structural QC checks into a per-construct fold verdict."""

from __future__ import annotations

from typing import Any, Callable

from loxdockaudit.alphafold_qc import (
    check_active_site_plddt,
    parse_pae_from_json,
    parse_plddt_from_pdb,
    summarize_alphafold_qc,
)
from loxdockaudit.models import StructuralQCResult
from loxdockaudit.pdb_parser import load_structure
from loxdockaudit.structural_qc import (
    check_active_site_sasa_proxy,
    check_disulfide_geometry,
    check_his_triad_geometry,
    check_lys_tyr_geometry,
)


def run_structural_qc(
    pdb_path: str,
    construct_id: str,
    lox_chain: str,
    his_resi: list[int],
    lys_resi: int,
    tyr_resi: int,
    disulfide_pairs: list[tuple[int, int]],
    active_site_resi: list[int],
    alphafold_pdb_path: str | None = None,
    pae_json_path: str | None = None,
    domain_a_resi: list[int] | None = None,
    domain_b_resi: list[int] | None = None,
    plddt_threshold: float = 70.0,
    his_max_ca_distance: float = 10.0,
    lys_tyr_max_cb_distance: float = 12.0,
    disulfide_max_sg_distance: float = 2.5,
) -> StructuralQCResult:
    """Run all structural QC checks and aggregate the result."""
    details: dict[str, Any] = {}
    warnings: list[str] = []

    structure = None
    try:
        structure = load_structure(pdb_path)
    except Exception as exc:
        warnings.append(f"Structure load failed: {exc}")
        details["structure_load"] = {"pdb_path": pdb_path, "warning": str(exc)}

    if structure is None:
        return _result(
            construct_id=construct_id,
            his_triad_pass=False,
            lys_tyr_pass=False,
            disulfide_pass=False,
            active_site_accessible=False,
            plddt_pass=None,
            pae_pass=None,
            warnings=warnings,
            details=details,
        )

    his_result = _run_check(
        "His triad",
        warnings,
        check_his_triad_geometry,
        structure,
        his_resi,
        lox_chain,
    )
    his_triad_pass = bool(his_result.get("geometry_plausible", False))
    if float(his_result.get("max_ca_distance", 0.0)) > his_max_ca_distance:
        his_triad_pass = False
        warning = "His triad CA distances exceed configured maximum."
        his_result["warning"] = _append_warning(his_result.get("warning"), warning)
        warnings.append(f"His triad: {warning}")
    details["his_triad"] = his_result

    lys_tyr_result = _run_check(
        "Lys-Tyr",
        warnings,
        check_lys_tyr_geometry,
        structure,
        lys_resi,
        tyr_resi,
        lox_chain,
        lys_tyr_max_cb_distance,
    )
    lys_tyr_pass = bool(lys_tyr_result.get("geometry_plausible", False))
    details["lys_tyr"] = lys_tyr_result

    if disulfide_pairs:
        disulfide_result = _run_check(
            "Disulfide",
            warnings,
            check_disulfide_geometry,
            structure,
            disulfide_pairs,
            lox_chain,
            disulfide_max_sg_distance,
        )
        disulfide_pass = bool(disulfide_result.get("all_intact", False))
    else:
        disulfide_result = {
            "pairs_checked": 0,
            "pairs_intact": 0,
            "pair_results": [],
            "all_intact": True,
            "warning": None,
            "skipped": True,
        }
        disulfide_pass = True
    details["disulfide"] = disulfide_result

    sasa_result = _run_check(
        "Active-site access",
        warnings,
        check_active_site_sasa_proxy,
        structure,
        active_site_resi,
        lox_chain,
    )
    active_site_accessible = sasa_result.get("accessibility_flag") != "buried"
    details["active_site_sasa_proxy"] = sasa_result

    plddt_pass: bool | None = None
    pae_pass: bool | None = None
    if alphafold_pdb_path is not None:
        plddt_pass, pae_pass = _run_alphafold_qc(
            alphafold_pdb_path=alphafold_pdb_path,
            pae_json_path=pae_json_path,
            lox_chain=lox_chain,
            active_site_resi=active_site_resi,
            domain_a_resi=domain_a_resi,
            domain_b_resi=domain_b_resi,
            plddt_threshold=plddt_threshold,
            warnings=warnings,
            details=details,
        )

    return _result(
        construct_id=construct_id,
        his_triad_pass=his_triad_pass,
        lys_tyr_pass=lys_tyr_pass,
        disulfide_pass=disulfide_pass,
        active_site_accessible=active_site_accessible,
        plddt_pass=plddt_pass,
        pae_pass=pae_pass,
        warnings=warnings,
        details=details,
    )


def format_qc_report(qc_result: StructuralQCResult) -> str:
    """Return a Markdown summary of a structural QC result."""
    details = qc_result.details
    disulfide_result = details.get("disulfide", {})
    sasa_result = details.get("active_site_sasa_proxy", {})
    active_site_plddt = details.get("active_site_plddt", {})
    pae_result = details.get("pae", {})

    lines = [
        f"## Structural QC: {qc_result.construct_id}",
        f"**Overall: {'PASS' if qc_result.fold_qc_pass else 'FOLD CORRUPTED'}**",
        "",
        "| Check | Result | Notes |",
        "|-------|--------|-------|",
        _report_row(
            "His triad geometry",
            qc_result.his_triad_pass,
            _his_note(details.get("his_triad", {})),
        ),
        _report_row(
            "Lys-Tyr geometry",
            qc_result.lys_tyr_pass,
            _distance_note(details.get("lys_tyr", {}), "cb_distance"),
        ),
        _report_row(
            "Disulfide geometry",
            None if disulfide_result.get("skipped") else qc_result.disulfide_pass,
            _disulfide_note(disulfide_result),
        ),
        _report_row(
            "Active-site access",
            qc_result.active_site_accessible,
            str(sasa_result.get("accessibility_flag", "unknown")),
        ),
        _report_row(
            "pLDDT (active site)",
            qc_result.plddt_pass,
            _mean_note(active_site_plddt.get("mean_plddt")),
        ),
        _report_row(
            "PAE (inter-domain)",
            qc_result.pae_pass,
            str(pae_result.get("inter_domain_confidence", "N/A")) if pae_result else "N/A",
        ),
        "",
        "### Warnings",
    ]

    if qc_result.warnings:
        lines.extend(f"- {warning}" for warning in qc_result.warnings)
    else:
        lines.append("None.")

    return "\n".join(lines).rstrip() + "\n"


def _run_alphafold_qc(
    *,
    alphafold_pdb_path: str,
    pae_json_path: str | None,
    lox_chain: str,
    active_site_resi: list[int],
    domain_a_resi: list[int] | None,
    domain_b_resi: list[int] | None,
    plddt_threshold: float,
    warnings: list[str],
    details: dict[str, Any],
) -> tuple[bool | None, bool | None]:
    try:
        af_structure = load_structure(alphafold_pdb_path)
    except Exception as exc:
        warnings.append(f"AlphaFold structure load failed: {exc}")
        details["alphafold_load"] = {"pdb_path": alphafold_pdb_path, "warning": str(exc)}
        return False, None

    plddt_result = _run_check("pLDDT", warnings, parse_plddt_from_pdb, af_structure, lox_chain)
    active_site_plddt = _run_check(
        "Active-site pLDDT",
        warnings,
        check_active_site_plddt,
        af_structure,
        lox_chain,
        active_site_resi,
        plddt_threshold,
    )
    details["plddt"] = plddt_result
    details["active_site_plddt"] = active_site_plddt
    plddt_pass = bool(active_site_plddt.get("qc_pass", False))

    pae_result = None
    pae_pass: bool | None = None
    if pae_json_path is not None:
        if domain_a_resi is None or domain_b_resi is None:
            pae_result = {
                "json_path": pae_json_path,
                "format_detected": "unknown",
                "matrix_size": None,
                "domain_a_resi": list(domain_a_resi or []),
                "domain_b_resi": list(domain_b_resi or []),
                "inter_domain_pae_mean": None,
                "inter_domain_pae_max": None,
                "inter_domain_confidence": "unknown",
                "warning": "PAE JSON provided without both domain residue lists.",
            }
            warnings.append("PAE: PAE JSON provided without both domain residue lists.")
        else:
            pae_result = _run_check(
                "PAE",
                warnings,
                parse_pae_from_json,
                pae_json_path,
                domain_a_resi,
                domain_b_resi,
            )
        details["pae"] = pae_result
        pae_pass = pae_result.get("inter_domain_confidence") in {"very_high", "high", "medium"}

    alphafold_summary = summarize_alphafold_qc(plddt_result, active_site_plddt, pae_result)
    details["alphafold_summary"] = alphafold_summary
    if alphafold_summary.get("fold_qc_warning"):
        warnings.append(f"AlphaFold: {alphafold_summary['fold_qc_warning']}")

    return plddt_pass, pae_pass


def _run_check(
    label: str,
    warnings: list[str],
    func: Callable[..., dict],
    *args: Any,
    **kwargs: Any,
) -> dict:
    try:
        result = func(*args, **kwargs)
    except Exception as exc:  # pragma: no cover - called functions are defensive
        exception_warning = f"{label} check failed: {exc}"
        warnings.append(exception_warning)
        return {"warning": exception_warning}

    warning = result.get("warning")
    if warning is not None:
        warnings.append(f"{label}: {str(warning)}")
    return result


def _result(
    *,
    construct_id: str,
    his_triad_pass: bool,
    lys_tyr_pass: bool,
    disulfide_pass: bool,
    active_site_accessible: bool,
    plddt_pass: bool | None,
    pae_pass: bool | None,
    warnings: list[str],
    details: dict,
) -> StructuralQCResult:
    checks = [
        his_triad_pass,
        lys_tyr_pass,
        disulfide_pass,
        active_site_accessible,
        True if plddt_pass is None else plddt_pass,
        True if pae_pass is None else pae_pass,
    ]
    fold_qc_pass = all(checks)
    return StructuralQCResult(
        construct_id=construct_id,
        his_triad_pass=his_triad_pass,
        lys_tyr_pass=lys_tyr_pass,
        disulfide_pass=disulfide_pass,
        active_site_accessible=active_site_accessible,
        plddt_pass=plddt_pass,
        pae_pass=pae_pass,
        fold_qc_pass=fold_qc_pass,
        fold_corrupted=not fold_qc_pass,
        warnings=warnings,
        details=details,
    )


def _append_warning(existing: object, new_warning: str) -> str:
    if existing is None:
        return new_warning
    return f"{existing} {new_warning}"


def _report_row(check: str, passed: bool | None, notes: str) -> str:
    if passed is None:
        result = "SKIPPED" if check == "Disulfide geometry" else "N/A"
    else:
        result = "PASS" if passed else "FAIL"
    return f"| {check} | {result} | {_escape_table_cell(notes)} |"


def _his_note(result: dict) -> str:
    max_distance = result.get("max_ca_distance")
    if max_distance is None:
        return "max CA distance: N/A"
    return f"max CA distance: {float(max_distance):.3f} A"


def _distance_note(result: dict, key: str) -> str:
    distance = result.get(key)
    if distance is None:
        return "distance: N/A"
    return f"distance: {float(distance):.3f} A"


def _disulfide_note(result: dict) -> str:
    if result.get("skipped"):
        return "no disulfide pairs configured"
    return f"{result.get('pairs_intact', 0)} / {result.get('pairs_checked', 0)} intact"


def _mean_note(value: object) -> str:
    if value is None:
        return "mean pLDDT: N/A"
    if isinstance(value, (int, float)):
        return f"mean pLDDT: {float(value):.3f}"
    return "mean pLDDT: N/A"


def _escape_table_cell(value: str) -> str:
    return value.replace("\n", " ").replace("|", "\\|")
