"""Configuration loading for LOXDockAudit."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from loxdockaudit.models import (
    ActiveSiteConfig,
    ConstructConfig,
    ControlEntry,
    DomainRange,
    ScreenConfig,
    StrictCriteria,
    TargetResidue,
)


def load_construct_config(yaml_path: str) -> ConstructConfig:
    """
    Load and validate a single construct YAML config file.
    """
    data = _load_yaml_mapping(yaml_path)
    return _construct_from_mapping(data, source=yaml_path)


def load_screen_config(yaml_path: str) -> ScreenConfig:
    """
    Load a screen-level YAML config.
    """
    data = _load_yaml_mapping(yaml_path)
    source = yaml_path
    required_paths = [
        ("screen_id",),
        ("candidate_construct_id",),
        ("constructs",),
    ]
    missing = [".".join(path) for path in required_paths if _get_nested(data, path) is None]
    if missing:
        raise ValueError(f"{source}: missing required field(s): {', '.join(missing)}")

    constructs = data.get("constructs")
    if not isinstance(constructs, list):
        raise ValueError(f"{source}: constructs must be a list")
    if not constructs:
        raise ValueError(f"{source}: constructs must contain at least one entry")

    base_dir = Path(yaml_path).resolve().parent
    control_entries = [
        _control_entry_from_mapping(
            entry,
            source=f"{yaml_path}:constructs[{index}]",
            base_dir=base_dir,
        )
        for index, entry in enumerate(constructs)
    ]

    candidate_construct_id = data["candidate_construct_id"]
    construct_ids = {entry.construct_id for entry in control_entries}
    if candidate_construct_id not in construct_ids:
        raise ValueError(
            f"{source}: candidate_construct_id '{candidate_construct_id}' "
            "does not appear in constructs"
        )
    if not any(entry.control_type == "baseline" for entry in control_entries):
        raise ValueError(f"{source}: at least one baseline control is required")

    return ScreenConfig(
        screen_id=data["screen_id"],
        candidate_construct_id=candidate_construct_id,
        constructs=control_entries,
        strict_criteria=_strict_criteria_from_mapping(data.get("strict_criteria", {})),
    )


def _load_yaml_mapping(yaml_path: str) -> dict[str, Any]:
    path = Path(yaml_path)
    try:
        with path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
    except OSError as exc:
        raise ValueError(f"Could not read YAML file {yaml_path}: {exc}") from exc
    except yaml.YAMLError as exc:
        raise ValueError(f"Could not parse YAML file {yaml_path}: {exc}") from exc

    if not isinstance(data, dict):
        raise ValueError(f"YAML file must contain a mapping: {yaml_path}")
    return data


def _construct_from_mapping(data: Any, source: str) -> ConstructConfig:
    if not isinstance(data, dict):
        raise ValueError(f"{source}: construct entry must be a mapping")

    required_paths = [
        ("construct_id",),
        ("construct_type",),
        ("models_dir",),
        ("receptor", "chains"),
        ("ligand", "chains"),
        ("active_site", "residues"),
        ("receptor", "target_residues"),
    ]
    missing = [".".join(path) for path in required_paths if _get_nested(data, path) is None]
    if missing:
        raise ValueError(f"{source}: missing required field(s): {', '.join(missing)}")

    receptor = data["receptor"]
    ligand = data["ligand"]
    active_site = data["active_site"]

    lox_domain = _domain_from_mapping(
        _first_present(data, ("lox_domain",), ("ligand", "lox_domain")),
        field_name="lox_domain",
        required=True,
        source=source,
    )
    if lox_domain is None:
        raise ValueError(f"{source}: missing required field: lox_domain")
    linker = _domain_from_mapping(
        _first_present(data, ("linker",), ("ligand", "linker")),
        field_name="linker",
        required=False,
        source=source,
    )
    cbd_domain = _domain_from_mapping(
        _first_present(data, ("cbd_domain",), ("ligand", "cbd_domain")),
        field_name="cbd_domain",
        required=False,
        source=source,
    )

    return ConstructConfig(
        construct_id=data["construct_id"],
        construct_type=data["construct_type"],
        models_dir=data["models_dir"],
        receptor_chains=list(receptor["chains"]),
        ligand_chains=list(ligand["chains"]),
        active_site=ActiveSiteConfig(residues=list(active_site["residues"])),
        target_residues=[
            _target_residue_from_mapping(target, source=source)
            for target in receptor["target_residues"]
        ],
        lox_domain=lox_domain,
        linker=linker,
        cbd_domain=cbd_domain,
        baseline_construct=data.get("baseline_construct"),
        productive_distance_threshold=float(
            data.get("productive_distance_threshold", 8.0)
        ),
        contact_distance_threshold=float(data.get("contact_distance_threshold", 4.0)),
        best_productive_rank_max=int(data.get("best_productive_rank_max", 3)),
        is_inactive=bool(data.get("is_inactive", False)),
        control_type=data.get("control_type"),
    )


def _target_residue_from_mapping(data: Any, source: str) -> TargetResidue:
    if not isinstance(data, dict):
        raise ValueError(f"{source}: receptor.target_residues entries must be mappings")

    missing = [field for field in ("chain", "resi", "resn") if field not in data]
    if missing:
        raise ValueError(
            f"{source}: missing required field(s) in receptor.target_residues: "
            f"{', '.join(missing)}"
        )

    return TargetResidue(
        chain=data["chain"],
        resi=int(data["resi"]),
        resn=data["resn"],
        atom=data.get("atom", "NZ"),
    )


def _control_entry_from_mapping(
    data: Any,
    source: str,
    base_dir: Path,
) -> ControlEntry:
    if not isinstance(data, dict):
        raise ValueError(f"{source}: construct entry must be a mapping")

    missing = [
        field
        for field in ("construct_id", "control_type", "config_path")
        if field not in data
    ]
    if missing:
        raise ValueError(f"{source}: missing required field(s): {', '.join(missing)}")

    allowed_control_types = {
        "baseline",
        "candidate",
        "scrambled",
        "polyK",
        "noCBD",
        "inactive",
        "benchmark",
    }
    control_type = str(data["control_type"])
    if control_type not in allowed_control_types:
        raise ValueError(f"{source}: unsupported control_type: {control_type}")

    raw_config_path = str(data["config_path"])
    config_path = Path(raw_config_path)
    resolved_config_path = config_path if config_path.is_absolute() else base_dir / config_path
    if not resolved_config_path.is_file():
        raise ValueError(f"{source}: config_path is not readable: {raw_config_path}")

    return ControlEntry(
        construct_id=data["construct_id"],
        control_type=control_type,
        config_path=str(resolved_config_path),
    )


def _strict_criteria_from_mapping(data: Any) -> StrictCriteria:
    if data is None:
        return StrictCriteria()
    if not isinstance(data, dict):
        raise ValueError("strict_criteria must be a mapping")

    return StrictCriteria(
        pass_count_greater_than_baseline=bool(
            data.get("pass_count_greater_than_baseline", True)
        ),
        best_distance_lte_baseline=bool(data.get("best_distance_lte_baseline", True)),
        best_productive_rank_max=int(data.get("best_productive_rank_max", 3)),
        cbd_must_contribute=bool(data.get("cbd_must_contribute", False)),
        must_beat_scrambled=bool(data.get("must_beat_scrambled", True)),
        must_beat_polyK=bool(data.get("must_beat_polyK", True)),
    )


def _domain_from_mapping(
    data: Any,
    field_name: str,
    required: bool,
    source: str,
) -> DomainRange | None:
    if data is None:
        if required:
            raise ValueError(f"{source}: missing required field: {field_name}")
        return None
    if not isinstance(data, dict):
        raise ValueError(f"{source}: {field_name} must be a mapping")

    missing = [field for field in ("start", "end") if field not in data]
    if missing:
        raise ValueError(
            f"{source}: missing required field(s) in {field_name}: {', '.join(missing)}"
        )

    return DomainRange(
        start=int(data["start"]),
        end=int(data["end"]),
        name=data.get("name", ""),
    )


def _get_nested(data: dict[str, Any], path: tuple[str, ...]) -> Any:
    current: Any = data
    for key in path:
        if not isinstance(current, dict) or key not in current:
            return None
        current = current[key]
    return current


def _first_present(data: dict[str, Any], *paths: tuple[str, ...]) -> Any:
    for path in paths:
        value = _get_nested(data, path)
        if value is not None:
            return value
    return None
