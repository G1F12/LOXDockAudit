from __future__ import annotations

import json
from pathlib import Path

import pytest

from loxdockaudit.af_parser import load_af_run, parse_af_json, parse_af_pdb


def test_parse_af_pdb_extracts_rank_and_chain_residues(tmp_path: Path) -> None:
    pdb_path = tmp_path / "rank_002_model.pdb"
    pdb_path.write_text(_af_pdb_text(), encoding="utf-8")

    model = parse_af_pdb(pdb_path)

    assert model.rank == 2
    assert model.chain_a_residues == frozenset({11})
    assert model.chain_b_residues == frozenset({100})
    assert model.iptm is None


def test_parse_af_json_colabfold_scores(tmp_path: Path) -> None:
    json_path = tmp_path / "rank_001_model.json"
    json_path.write_text(
        json.dumps(
            {
                "iptm": 0.82,
                "ptm": 0.71,
                "ranking_confidence": 0.79,
                "plddt": [80.0, 90.0],
            }
        ),
        encoding="utf-8",
    )

    scores = parse_af_json(json_path)

    assert scores.iptm == pytest.approx(0.82)
    assert scores.ptm == pytest.approx(0.71)
    assert scores.ranking_score == pytest.approx(0.79)
    assert scores.per_residue_plddt == {1: 80.0, 2: 90.0}


def test_load_af_run_uses_natural_sort_and_sidecar_scores(tmp_path: Path) -> None:
    for name in ("rank_010_model.pdb", "rank_002_model.pdb"):
        (tmp_path / name).write_text(_af_pdb_text(), encoding="utf-8")
    (tmp_path / "rank_002_model.json").write_text(
        json.dumps({"iptm": 0.9, "ptm": 0.8}),
        encoding="utf-8",
    )

    models = load_af_run(tmp_path)

    assert [model.rank for model in models] == [2, 10]
    assert models[0].iptm == pytest.approx(0.9)


def _af_pdb_text() -> str:
    return "\n".join(
        [
            "ATOM      1  CB  LYS A  11       3.000   0.000   0.000  1.00 80.00           C  ",
            "ATOM      2  NZ  LYS A  11       4.000   0.000   0.000  1.00 80.00           N  ",
            "ATOM      3  CA  HIS D 100       0.000   0.000   0.000  1.00 80.00           C  ",
            "TER",
            "END",
        ]
    )
