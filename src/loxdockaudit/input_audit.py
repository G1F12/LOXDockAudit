"""Input validation utilities for LOXDockAudit PDB screens."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any

from loxdockaudit.pdb_parser import (
    compute_sha256,
    get_chain_ids,
    get_sequence_from_chain,
    load_structure,
)
from loxdockaudit.sorting import natural_sort_key


def audit_single_file(
    pdb_path: str,
    expected_chains: list[str],
    expected_n_terminus: str | None = None,
    n_terminus_length: int = 10,
) -> dict[str, Any]:
    """
    Validate one PDB file before docking analysis.
    """
    result: dict[str, Any] = {
        "path": pdb_path,
        "exists": False,
        "parseable": False,
        "sha256": None,
        "chains_found": [],
        "chains_expected": expected_chains,
        "chains_ok": False,
        "n_terminus_detected": None,
        "n_terminus_expected": expected_n_terminus,
        "n_terminus_ok": None,
        "warnings": [],
        "errors": [],
        "valid": False,
    }

    path = Path(pdb_path)
    if not path.is_file():
        result["errors"].append("file not found")
        return result

    result["exists"] = True
    try:
        result["sha256"] = compute_sha256(pdb_path)
    except OSError as exc:
        result["errors"].append(f"failed to compute sha256: {exc}")

    try:
        structure = load_structure(pdb_path)
    except Exception as exc:
        result["errors"].append(f"failed to parse PDB: {exc}")
        return result

    result["parseable"] = True
    chains_found = get_chain_ids(structure)
    result["chains_found"] = chains_found

    missing_chains = [chain_id for chain_id in expected_chains if chain_id not in chains_found]
    result["chains_ok"] = not missing_chains
    if missing_chains:
        result["errors"].append(f"missing expected chains: {', '.join(missing_chains)}")

    unexpected_chains = [chain_id for chain_id in chains_found if chain_id not in expected_chains]
    if unexpected_chains:
        result["warnings"].append(f"unexpected chains found: {', '.join(unexpected_chains)}")

    if expected_n_terminus is not None:
        detected = _detect_n_terminus(
            structure,
            chains_to_check=expected_chains,
            n_terminus_length=n_terminus_length,
        )
        expected = expected_n_terminus[:n_terminus_length]
        result["n_terminus_detected"] = detected
        result["n_terminus_ok"] = detected == expected
        if detected is None:
            result["errors"].append("could not detect N-terminus sequence")
        elif detected != expected:
            result["errors"].append(
                f"N-terminus mismatch: expected {expected}, found {detected}"
            )

    result["valid"] = (
        result["parseable"]
        and result["chains_ok"]
        and not result["errors"]
    )
    return result


def audit_screen(
    construct_configs: list[dict[str, Any]],
) -> dict[str, Any]:
    """
    Audit all construct model PDBs and detect duplicate SHA256 hashes.
    """
    audit_result: dict[str, Any] = {
        "constructs": {},
        "screen_valid": False,
        "critical_errors": [],
    }
    sha_entries: dict[str, list[tuple[str, str]]] = defaultdict(list)

    for config in construct_configs:
        construct_id = config["construct_id"]
        models_dir = Path(config["models_dir"])
        expected_chains = list(config.get("ligand_chains", [])) + list(
            config.get("receptor_chains", [])
        )
        construct_result: dict[str, Any] = {
            "models": [],
            "duplicate_sha256": [],
            "cross_construct_duplicates": [],
            "construct_valid": False,
        }

        pdb_paths = (
            sorted(
                models_dir.glob("*.pdb"),
                key=lambda path: natural_sort_key(str(path)),
            )
            if models_dir.is_dir()
            else []
        )
        if not models_dir.is_dir():
            audit_result["critical_errors"].append(
                f"{construct_id}: models_dir not found: {models_dir}"
            )
        elif not pdb_paths:
            audit_result["critical_errors"].append(
                f"{construct_id}: no PDB files found in {models_dir}"
            )

        for pdb_path in pdb_paths:
            model_audit = audit_single_file(
                str(pdb_path),
                expected_chains=expected_chains,
                expected_n_terminus=config.get("expected_n_terminus"),
            )
            construct_result["models"].append(model_audit)
            if model_audit["sha256"] is not None:
                sha_entries[model_audit["sha256"]].append((construct_id, str(pdb_path)))

        audit_result["constructs"][construct_id] = construct_result

    duplicate_by_construct: dict[str, set[str]] = defaultdict(set)
    cross_duplicates_by_construct: dict[str, list[dict[str, str]]] = defaultdict(list)

    for sha256, entries in sha_entries.items():
        if len(entries) < 2:
            continue

        constructs_with_sha = {construct_id for construct_id, _ in entries}
        for construct_id, _ in entries:
            duplicate_by_construct[construct_id].add(sha256)

        if len(constructs_with_sha) > 1:
            for index, first in enumerate(entries):
                for second in entries[index + 1 :]:
                    first_construct, first_file = first
                    second_construct, second_file = second
                    if first_construct == second_construct:
                        continue

                    duplicate = {
                        "sha256": sha256,
                        "file_a": first_file,
                        "file_b": second_file,
                    }
                    audit_result["critical_errors"].append(
                        "cross-construct duplicate SHA256 "
                        f"{sha256}: {first_file} == {second_file}"
                    )
                    cross_duplicates_by_construct[first_construct].append(duplicate)
                    cross_duplicates_by_construct[second_construct].append(duplicate)

    for construct_id, construct_result in audit_result["constructs"].items():
        duplicate_hashes = sorted(duplicate_by_construct[construct_id])
        construct_result["duplicate_sha256"] = duplicate_hashes
        construct_result["cross_construct_duplicates"] = cross_duplicates_by_construct[
            construct_id
        ]
        models = construct_result["models"]
        construct_result["construct_valid"] = (
            bool(models)
            and all(model["valid"] for model in models)
            and not duplicate_hashes
            and not construct_result["cross_construct_duplicates"]
        )

    audit_result["screen_valid"] = (
        bool(audit_result["constructs"])
        and all(
            construct["construct_valid"]
            for construct in audit_result["constructs"].values()
        )
        and not audit_result["critical_errors"]
    )
    return audit_result


def format_audit_report(audit_result: dict[str, Any]) -> str:
    """
    Return a human-readable Markdown summary of audit_screen output.
    """
    lines = [
        "# LOXDockAudit Input Audit",
        "",
        f"Screen valid: {_format_bool(audit_result.get('screen_valid', False))}",
        "",
    ]

    critical_errors = audit_result.get("critical_errors", [])
    if critical_errors:
        lines.extend(["## Critical Errors", ""])
        lines.extend(f"- {error}" for error in critical_errors)
        lines.append("")

    lines.extend(["## Constructs", ""])
    for construct_id, construct in audit_result.get("constructs", {}).items():
        lines.extend(
            [
                f"### {construct_id}",
                "",
                f"- Construct valid: {_format_bool(construct.get('construct_valid', False))}",
                f"- Models audited: {len(construct.get('models', []))}",
            ]
        )

        duplicate_sha256 = construct.get("duplicate_sha256", [])
        if duplicate_sha256:
            lines.append("- Duplicate SHA256 hashes:")
            lines.extend(f"  - `{sha256}`" for sha256 in duplicate_sha256)

        cross_duplicates = construct.get("cross_construct_duplicates", [])
        if cross_duplicates:
            lines.append("- Cross-construct duplicate files:")
            for duplicate in cross_duplicates:
                lines.append(
                    "  - "
                    f"`{duplicate['sha256']}`: {duplicate['file_a']} == "
                    f"{duplicate['file_b']}"
                )

        invalid_models = [
            model for model in construct.get("models", []) if not model.get("valid", False)
        ]
        if invalid_models:
            lines.append("- Invalid models:")
            for model in invalid_models:
                problems = model.get("errors", []) + model.get("warnings", [])
                problem_text = "; ".join(problems) if problems else "unknown issue"
                lines.append(f"  - {model.get('path')}: {problem_text}")

        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def _detect_n_terminus(
    structure: Any,
    chains_to_check: list[str],
    n_terminus_length: int,
) -> str | None:
    for chain_id in chains_to_check:
        sequence = get_sequence_from_chain(structure, chain_id)
        if sequence:
            return sequence[:n_terminus_length]
    return None


def _format_bool(value: bool) -> str:
    return "PASS" if value else "FAIL"
