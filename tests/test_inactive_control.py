from __future__ import annotations

from pathlib import Path

import yaml

from loxdockaudit.inactive_control import (
    CatalyticResidueExpectation,
    evaluate_catalytic_identity,
    interpret_inactive_control,
    load_inactive_control_config,
    run_inactive_control_analysis,
)
from loxdockaudit.models import ConstructSummary
from tests.conftest import write_pose


def test_evaluate_catalytic_identity_active_his(tmp_path: Path) -> None:
    pdb_path = tmp_path / "active.pdb"
    pdb_path.write_text(_identity_pdb_text("HIS"), encoding="utf-8")

    result = evaluate_catalytic_identity(
        pdb_path,
        "D",
        [
            CatalyticResidueExpectation("His292", 292, "HIS"),
            CatalyticResidueExpectation("His294", 294, "HIS"),
            CatalyticResidueExpectation("His296", 296, "HIS"),
        ],
    )

    assert result["expected_identity_pass"] is True
    assert result["catalytic_status"] == "intact"


def test_evaluate_catalytic_identity_inactive_ala(tmp_path: Path) -> None:
    pdb_path = tmp_path / "inactive.pdb"
    pdb_path.write_text(_identity_pdb_text("ALA"), encoding="utf-8")

    result = evaluate_catalytic_identity(
        pdb_path,
        "D",
        [
            CatalyticResidueExpectation("Ala292", 292, "ALA"),
            CatalyticResidueExpectation("Ala294", 294, "ALA"),
            CatalyticResidueExpectation("Ala296", 296, "ALA"),
        ],
    )

    assert result["expected_identity_pass"] is True
    assert result["catalytic_status"] == "disrupted"


def test_interpret_inactive_control_active_better() -> None:
    active = _summary("active", productive_count=2, best_distance=5.0)
    inactive = _summary("inactive", productive_count=0, best_distance=9.0)

    assert "more favorable productive geometry" in interpret_inactive_control(
        active,
        inactive,
    )


def test_run_inactive_control_analysis_writes_outputs(tmp_path: Path) -> None:
    active_models = tmp_path / "active_models"
    inactive_models = tmp_path / "inactive_models"
    active_models.mkdir()
    inactive_models.mkdir()
    write_pose(active_models / "model_1.pdb", 6.0)
    write_pose(active_models / "model_2.pdb", 12.0)
    _write_inactive_pose(inactive_models / "model_1.pdb", 10.0)
    _write_inactive_pose(inactive_models / "model_2.pdb", 11.0)

    active_config = _write_construct_config(tmp_path, "active", active_models)
    inactive_config = _write_construct_config(
        tmp_path,
        "inactive",
        inactive_models,
        inactive=True,
    )
    comparison_config = tmp_path / "comparison.yaml"
    comparison_config.write_text(
        yaml.safe_dump(
            {
                "comparison_id": "inactive_test",
                "top_n": 2,
                "mutations": ["H124A", "H126A", "H128A"],
                "active": {
                    "construct_id": "active",
                    "config_path": str(active_config),
                    "catalytic_residues": [
                        {"label": "His124", "resi": 124, "expected_resn": "HIS"},
                        {"label": "His126", "resi": 126, "expected_resn": "HIS"},
                        {"label": "His128", "resi": 128, "expected_resn": "HIS"},
                    ],
                },
                "inactive": {
                    "construct_id": "inactive",
                    "config_path": str(inactive_config),
                    "catalytic_residues": [
                        {"label": "Ala124", "resi": 124, "expected_resn": "ALA"},
                        {"label": "Ala126", "resi": 126, "expected_resn": "ALA"},
                        {"label": "Ala128", "resi": 128, "expected_resn": "ALA"},
                    ],
                },
            }
        ),
        encoding="utf-8",
    )

    out_dir = tmp_path / "out"
    result = run_inactive_control_analysis(
        load_inactive_control_config(str(comparison_config)),
        out_dir,
    )

    assert result["active_summary"].productive_count == 1
    assert result["inactive_summary"].productive_count == 0
    assert result["inactive_catalytic"]["catalytic_status"] == "disrupted"
    assert (out_dir / "inactive_test_comparison.csv").exists()
    assert (out_dir / "inactive_test_supplement.md").exists()
    assert (out_dir / "inactive_test_distance_histogram.svg").exists()


def _identity_pdb_text(resname: str) -> str:
    return (
        f"ATOM      1  CA  {resname} D 292       0.000   0.000   0.000  1.00 20.00           C  \n"
        f"ATOM      2  CA  {resname} D 294       1.000   0.000   0.000  1.00 20.00           C  \n"
        f"ATOM      3  CA  {resname} D 296       2.000   0.000   0.000  1.00 20.00           C  \n"
        "TER\n"
        "END\n"
    )


def _write_inactive_pose(path: Path, distance: float) -> None:
    write_pose(path, distance)
    text = path.read_text(encoding="utf-8")
    text = text.replace("HIS D 124", "ALA D 124")
    text = text.replace("HIS D 126", "ALA D 126")
    text = text.replace("HIS D 128", "ALA D 128")
    path.write_text(text, encoding="utf-8")


def _write_construct_config(
    tmp_path: Path,
    construct_id: str,
    models_dir: Path,
    inactive: bool = False,
) -> Path:
    path = tmp_path / f"{construct_id}.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "construct_id": construct_id,
                "construct_type": "lox_inactive_control" if inactive else "lox_baseline",
                "models_dir": str(models_dir),
                "is_inactive": inactive,
                "receptor": {
                    "chains": ["A", "B"],
                    "target_residues": [
                        {"chain": "A", "resi": 11, "resn": "LYS", "atom": "NZ"}
                    ],
                },
                "ligand": {"chains": ["D"]},
                "active_site": {
                    "residues": [
                        {"label": "His124", "resi": 124},
                        {"label": "His126", "resi": 126},
                        {"label": "His128", "resi": 128},
                        {"label": "Lys152", "resi": 152},
                        {"label": "Tyr187", "resi": 187},
                    ]
                },
                "lox_domain": {"start": 1, "end": 200},
                "structural_qc": {
                    "his_resi": [124, 126, 128],
                    "lys_resi": 152,
                    "tyr_resi": 187,
                    "disulfide_pairs": [],
                },
            }
        ),
        encoding="utf-8",
    )
    return path


def _summary(
    construct_id: str,
    productive_count: int,
    best_distance: float,
) -> ConstructSummary:
    return ConstructSummary(
        construct_id=construct_id,
        productive_count=productive_count,
        total_poses=2,
        best_distance=best_distance,
        best_productive_rank=1 if productive_count else None,
        lox_contact_frequency=0.0,
        cbd_contact_frequency=0.0,
        cbd_coupling=0.0,
        strict_pass=None,
    )
