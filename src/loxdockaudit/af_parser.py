"""Parsers for ColabFold and AlphaFold-Multimer output files."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from loxdockaudit.pdb_parser import get_chain_ids, load_structure
from loxdockaudit.sorting import natural_sort_key


@dataclass(frozen=True)
class AFScores:
    """Confidence and ranking scores parsed from AF-Multimer JSON output."""

    iptm: float | None
    ptm: float | None
    ranking_score: float | None
    per_residue_plddt: dict[int, float]


@dataclass(frozen=True)
class AFMultimerModel:
    """One ranked AF-Multimer model and its parsed metadata."""

    pdb_path: Path
    rank: int
    iptm: float | None
    ptm: float | None
    chain_a_residues: frozenset[int]
    chain_b_residues: frozenset[int]


def parse_af_pdb(pdb_path: Path) -> AFMultimerModel:
    """Parse an AF-Multimer PDB file into an ``AFMultimerModel``.

    The parser records rank from common ColabFold/AlphaFold filename patterns
    when available and captures residue numbers for the first two chains in the
    model. Confidence scores are left as ``None`` because they normally live in
    JSON sidecar files.
    """
    structure = load_structure(str(pdb_path))
    chain_ids = get_chain_ids(structure)
    chain_a = chain_ids[0] if chain_ids else ""
    chain_b = chain_ids[1] if len(chain_ids) > 1 else ""
    return AFMultimerModel(
        pdb_path=pdb_path,
        rank=_rank_from_name(pdb_path),
        iptm=None,
        ptm=None,
        chain_a_residues=_residue_numbers(structure, chain_a),
        chain_b_residues=_residue_numbers(structure, chain_b),
    )


def parse_af_json(json_path: Path) -> AFScores:
    """Parse ColabFold/AF-Multimer JSON confidence and ranking scores."""
    with json_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"AF JSON must contain an object: {json_path}")

    return AFScores(
        iptm=_optional_float(_first_present(payload, "iptm", "ipTM")),
        ptm=_optional_float(_first_present(payload, "ptm", "pTM")),
        ranking_score=_optional_float(
            _first_present(
                payload,
                "ranking_score",
                "ranking_confidence",
                "rank_score",
            )
        ),
        per_residue_plddt=_parse_plddt(payload),
    )


def load_af_run(run_dir: Path) -> list[AFMultimerModel]:
    """Load all ranked AF-Multimer models from a run directory."""
    if not run_dir.is_dir():
        raise FileNotFoundError(f"AF run directory not found: {run_dir}")

    pdb_paths = sorted(run_dir.glob("*.pdb"), key=lambda path: natural_sort_key(str(path)))
    models = []
    for pdb_path in pdb_paths:
        model = parse_af_pdb(pdb_path)
        scores = _matching_scores(pdb_path)
        if scores is not None:
            model = AFMultimerModel(
                pdb_path=model.pdb_path,
                rank=model.rank,
                iptm=scores.iptm,
                ptm=scores.ptm,
                chain_a_residues=model.chain_a_residues,
                chain_b_residues=model.chain_b_residues,
            )
        models.append(model)
    return models


def _rank_from_name(path: Path) -> int:
    name = path.stem.lower()
    for pattern in (r"rank[_-]?(\d+)", r"model[_-]?(\d+)", r"_(\d+)$"):
        match = re.search(pattern, name)
        if match is not None:
            return int(match.group(1))
    return 1


def _residue_numbers(structure: Any, chain_id: str) -> frozenset[int]:
    if not chain_id:
        return frozenset()
    model = next(structure.get_models(), None)
    if model is None or chain_id not in model:
        return frozenset()
    return frozenset(
        int(residue.id[1])
        for residue in model[chain_id]
        if residue.id[0] == " "
    )


def _matching_scores(pdb_path: Path) -> AFScores | None:
    candidates = [
        pdb_path.with_suffix(".json"),
        pdb_path.with_name(f"{pdb_path.stem}_scores.json"),
        pdb_path.with_name(f"scores_{pdb_path.stem}.json"),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return parse_af_json(candidate)
    return None


def _first_present(payload: dict[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in payload:
            return payload[key]
    return None


def _optional_float(value: Any) -> float | None:
    if value is None:
        return None
    return float(value)


def _parse_plddt(payload: dict[str, Any]) -> dict[int, float]:
    raw = _first_present(payload, "plddt", "predicted_lddt", "per_residue_plddt")
    if raw is None:
        return {}
    if isinstance(raw, dict):
        return {int(key): float(value) for key, value in raw.items()}
    if isinstance(raw, list):
        return {index: float(value) for index, value in enumerate(raw, start=1)}
    return {}
