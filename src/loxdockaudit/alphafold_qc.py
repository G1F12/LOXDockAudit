"""AlphaFold confidence quality-control checks for LOXDockAudit."""

from __future__ import annotations

import json
from typing import Any

from loxdockaudit import pdb_parser


def parse_plddt_from_pdb(
    structure: Any,
    chain: str,
    resi_list: list[int] | None = None,
) -> dict:
    """Extract per-residue pLDDT scores from CA atom B-factors."""
    result: dict[str, Any] = {
        "chain": chain,
        "residues_requested": list(resi_list) if resi_list is not None else None,
        "residues_found": [],
        "plddt_per_residue": {},
        "mean_plddt": None,
        "min_plddt": None,
        "max_plddt": None,
        "n_low_confidence": 0,
        "n_medium_confidence": 0,
        "n_high_confidence": 0,
        "n_very_high_confidence": 0,
        "overall_confidence": "unknown",
        "warning": None,
    }

    try:
        ca_atoms = pdb_parser.get_atoms_by_selection(
            structure,
            chains=[chain],
            resi_list=list(resi_list) if resi_list is not None else None,
            atom_names=["CA"],
        )
        requested_order = list(resi_list) if resi_list is not None else None
        plddt_by_resi = {
            _residue_number(atom): _bfactor(atom)
            for atom in ca_atoms
        }
        if requested_order is None:
            residues_found = sorted(plddt_by_resi)
        else:
            residues_found = [resi for resi in requested_order if resi in plddt_by_resi]

        result["residues_found"] = residues_found
        result["plddt_per_residue"] = {
            resi: plddt_by_resi[resi]
            for resi in residues_found
        }

        scores = list(result["plddt_per_residue"].values())
        if not scores:
            result["warning"] = "No CA pLDDT scores found."
            return result

        result["mean_plddt"] = sum(scores) / len(scores)
        result["min_plddt"] = min(scores)
        result["max_plddt"] = max(scores)
        result["n_low_confidence"] = sum(1 for score in scores if score < 50.0)
        result["n_medium_confidence"] = sum(1 for score in scores if 50.0 <= score < 70.0)
        result["n_high_confidence"] = sum(1 for score in scores if 70.0 <= score < 90.0)
        result["n_very_high_confidence"] = sum(1 for score in scores if score >= 90.0)
        result["overall_confidence"] = _plddt_confidence(float(result["mean_plddt"]))

        if resi_list is not None and len(residues_found) != len(resi_list):
            missing = [resi for resi in resi_list if resi not in plddt_by_resi]
            result["warning"] = "Missing CA pLDDT scores for residues: " + _join_ints(missing) + "."

    except Exception as exc:  # pragma: no cover - defensive API boundary
        result["warning"] = f"pLDDT parsing failed: {exc}"

    return result


def parse_pae_from_json(
    json_path: str,
    domain_a_resi: list[int],
    domain_b_resi: list[int],
) -> dict:
    """Parse inter-domain PAE from ColabFold or AlphaFold JSON output."""
    result: dict[str, Any] = {
        "json_path": json_path,
        "format_detected": "unknown",
        "matrix_size": None,
        "domain_a_resi": list(domain_a_resi),
        "domain_b_resi": list(domain_b_resi),
        "inter_domain_pae_mean": None,
        "inter_domain_pae_max": None,
        "inter_domain_confidence": "unknown",
        "warning": None,
    }

    try:
        with open(json_path, encoding="utf-8") as handle:
            payload = json.load(handle)

        matrix, format_detected = _extract_pae_matrix(payload)
        result["format_detected"] = format_detected
        if matrix is None:
            result["warning"] = "No supported PAE matrix found in JSON."
            return result

        matrix_size = len(matrix)
        result["matrix_size"] = matrix_size

        if not domain_a_resi or not domain_b_resi:
            result["warning"] = "Domain residue index lists must not be empty."
            return result

        values = _extract_pae_block(matrix, domain_a_resi, domain_b_resi)
        if not values:
            result["warning"] = "No valid PAE values found for requested domain residues."
            return result

        pae_mean = sum(values) / len(values)
        result["inter_domain_pae_mean"] = pae_mean
        result["inter_domain_pae_max"] = max(values)
        result["inter_domain_confidence"] = _pae_confidence(pae_mean)

        expected_values = len(domain_a_resi) * len(domain_b_resi)
        if len(values) != expected_values:
            result["warning"] = "Some requested PAE residues were outside the matrix."

    except FileNotFoundError:
        result["warning"] = f"PAE JSON file not found: {json_path}"
    except Exception as exc:  # pragma: no cover - defensive API boundary
        result["warning"] = f"PAE JSON parsing failed: {exc}"

    return result


def check_active_site_plddt(
    structure: Any,
    chain: str,
    active_site_resi: list[int],
    threshold: float = 70.0,
) -> dict:
    """Extract pLDDT specifically for active-site residues."""
    result: dict[str, Any] = {
        "active_site_resi": list(active_site_resi),
        "plddt_per_residue": {},
        "mean_plddt": None,
        "all_above_threshold": False,
        "below_threshold_resi": [],
        "threshold_used": threshold,
        "qc_pass": False,
        "warning": None,
    }

    try:
        plddt_result = parse_plddt_from_pdb(structure, chain, active_site_resi)
        plddt_per_residue = dict(plddt_result["plddt_per_residue"])
        result["plddt_per_residue"] = plddt_per_residue

        scores = list(plddt_per_residue.values())
        if not scores:
            result["warning"] = "No active-site pLDDT scores found."
            return result

        mean_plddt = sum(scores) / len(scores)
        below_threshold = [
            resi
            for resi in active_site_resi
            if resi in plddt_per_residue and plddt_per_residue[resi] < threshold
        ]
        missing = [resi for resi in active_site_resi if resi not in plddt_per_residue]

        result["mean_plddt"] = mean_plddt
        result["below_threshold_resi"] = below_threshold
        result["all_above_threshold"] = not below_threshold and not missing
        result["qc_pass"] = mean_plddt >= threshold and not missing

        warnings = []
        if missing:
            warnings.append("Missing active-site CA pLDDT scores for residues: " + _join_ints(missing) + ".")
        if below_threshold:
            warnings.append("Active-site residues below pLDDT threshold: " + _join_ints(below_threshold) + ".")
        result["warning"] = " ".join(warnings) if warnings else None

    except Exception as exc:  # pragma: no cover - defensive API boundary
        result["warning"] = f"Active-site pLDDT check failed: {exc}"

    return result


def summarize_alphafold_qc(
    plddt_result: dict,
    active_site_plddt_result: dict,
    pae_result: dict | None = None,
) -> dict:
    """Combine pLDDT and PAE results into a single QC summary."""
    inter_domain_confidence = None
    if pae_result is not None:
        inter_domain_confidence = pae_result.get("inter_domain_confidence", "unknown")

    active_site_pass = bool(active_site_plddt_result.get("qc_pass", False))
    pae_pass = pae_result is None or inter_domain_confidence in {"very_high", "high", "medium"}
    fold_qc_pass = active_site_pass and pae_pass

    result = {
        "overall_plddt_confidence": plddt_result.get("overall_confidence", "unknown"),
        "active_site_plddt_pass": active_site_pass,
        "active_site_mean_plddt": active_site_plddt_result.get("mean_plddt"),
        "inter_domain_pae_confidence": inter_domain_confidence,
        "fold_qc_pass": fold_qc_pass,
        "fold_qc_warning": None,
    }

    if not fold_qc_pass:
        warnings = []
        if not active_site_pass:
            warnings.append("Active-site pLDDT QC failed.")
        if not pae_pass:
            warnings.append("Inter-domain PAE confidence is low or unknown.")
        result["fold_qc_warning"] = " ".join(warnings)

    return result


def _extract_pae_matrix(payload: Any) -> tuple[list[list[float]] | None, str]:
    if isinstance(payload, dict):
        if "pae" in payload:
            matrix = _coerce_matrix(payload["pae"])
            return matrix, "colabfold" if matrix is not None else "unknown"
        if "predicted_aligned_error" in payload:
            matrix = _coerce_matrix(payload["predicted_aligned_error"])
            return matrix, "alphafold_db" if matrix is not None else "unknown"
    if isinstance(payload, list):
        for item in payload:
            if isinstance(item, dict) and "predicted_aligned_error" in item:
                matrix = _coerce_matrix(item["predicted_aligned_error"])
                return matrix, "af2_scores" if matrix is not None else "unknown"
    return None, "unknown"


def _coerce_matrix(value: Any) -> list[list[float]] | None:
    if not isinstance(value, list) or not value:
        return None

    matrix = []
    for row in value:
        if not isinstance(row, list):
            return None
        try:
            matrix.append([float(cell) for cell in row])
        except (TypeError, ValueError):
            return None
    return matrix


def _extract_pae_block(
    matrix: list[list[float]],
    domain_a_resi: list[int],
    domain_b_resi: list[int],
) -> list[float]:
    values = []
    for row_resi in domain_a_resi:
        row_index = row_resi - 1
        if row_index < 0 or row_index >= len(matrix):
            continue
        row = matrix[row_index]
        for col_resi in domain_b_resi:
            col_index = col_resi - 1
            if col_index < 0 or col_index >= len(row):
                continue
            values.append(row[col_index])
    return values


def _plddt_confidence(mean_plddt: float) -> str:
    if mean_plddt >= 90.0:
        return "very_high"
    if mean_plddt >= 70.0:
        return "high"
    if mean_plddt >= 50.0:
        return "medium"
    return "low"


def _pae_confidence(mean_pae: float) -> str:
    if mean_pae < 5.0:
        return "very_high"
    if mean_pae < 10.0:
        return "high"
    if mean_pae <= 20.0:
        return "medium"
    return "low"


def _join_ints(values: list[int]) -> str:
    return ", ".join(str(value) for value in values)


def _residue_number(atom: Any) -> int:
    return int(atom.get_parent().id[1])


def _bfactor(atom: Any) -> float:
    return float(atom.get_bfactor())
