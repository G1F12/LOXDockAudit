from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pytest
import yaml
from click.testing import CliRunner

from loxdockaudit.cli import cli
from loxdockaudit.input_audit import audit_screen


ACTIVE_SITE_RESI = {124, 126, 128, 152, 187}


def test_full_cli_pipeline(tmp_path: Path, fixture_pdb_path: Path) -> None:
    paths = _create_integration_scenario(tmp_path, fixture_pdb_path)

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "--config",
            str(paths["baseline_config"]),
            "--out",
            str(paths["out"]),
            "--top-n",
            "3",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Productive poses: 2/3" in result.output
    assert "Best distance:" in result.output
    assert "Best productive rank: 1" in result.output

    poses_csv = paths["out"] / "test_lox_baseline_poses.csv"
    summary_csv = paths["out"] / "test_lox_baseline_summary.csv"
    report_md = paths["out"] / "test_lox_baseline_report.md"
    assert poses_csv.exists()
    assert summary_csv.exists()
    assert report_md.exists()

    poses = pd.read_csv(poses_csv)
    assert len(poses) == 3
    assert poses["productive"].tolist() == [True, False, True]

    summary = pd.read_csv(summary_csv)
    assert summary.loc[0, "productive_fraction"] == pytest.approx(0.667, abs=0.01)

    report = report_md.read_text(encoding="utf-8")
    assert "## Executive Decision" in report
    assert "## Docking Pose Summary" in report


def test_audit_screen_detects_cross_construct_duplicate(
    tmp_path: Path,
    fixture_pdb_path: Path,
) -> None:
    paths = _create_integration_scenario(tmp_path, fixture_pdb_path)
    result = audit_screen(
        [
            _audit_config("test_lox_baseline", paths["baseline_models"]),
            _audit_config("test_lox_fusion", paths["fusion_models"]),
        ]
    )
    assert result["screen_valid"] is False
    assert result["constructs"]["test_lox_baseline"]["cross_construct_duplicates"]
    assert result["constructs"]["test_lox_fusion"]["cross_construct_duplicates"]


def test_check_config_valid_and_missing_file(tmp_path: Path, fixture_pdb_path: Path) -> None:
    paths = _create_integration_scenario(tmp_path, fixture_pdb_path)
    runner = CliRunner()

    valid = runner.invoke(cli, ["check-config", "--config", str(paths["baseline_config"])])
    assert valid.exit_code == 0
    assert "test_lox_baseline" in valid.output

    missing = runner.invoke(cli, ["check-config", "--config", str(tmp_path / "missing.yaml")])
    assert missing.exit_code == 1
    assert "ERROR" in missing.output


def test_screen_cli_pipeline(tmp_path: Path, fixture_pdb_path: Path) -> None:
    paths = _create_integration_scenario(tmp_path, fixture_pdb_path)
    scrambled_models = tmp_path / "models" / "scrambled"
    scrambled_models.mkdir()
    for index, distance in enumerate([4.0, 5.0, 6.0], start=1):
        _write_shifted_fixture(
            fixture_pdb_path,
            scrambled_models / f"model_{index}.pdb",
            distance,
        )
    scrambled_config = paths["configs"] / "scrambled.yaml"
    scrambled_config.write_text(
        _config_yaml("test_scrambled", "control", scrambled_models),
        encoding="utf-8",
    )
    screen_config = paths["configs"] / "screen.yaml"
    screen_config.write_text(
        yaml.safe_dump(
            {
                "screen_id": "integration_screen",
                "candidate_construct_id": "test_lox_fusion",
                "constructs": [
                    {
                        "construct_id": "test_lox_baseline",
                        "control_type": "baseline",
                        "config_path": str(paths["baseline_config"]),
                    },
                    {
                        "construct_id": "test_lox_fusion",
                        "control_type": "candidate",
                        "config_path": str(paths["fusion_config"]),
                    },
                    {
                        "construct_id": "test_scrambled",
                        "control_type": "scrambled",
                        "config_path": str(scrambled_config),
                    },
                ],
                "strict_criteria": {"must_beat_polyK": False},
            }
        ),
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        cli,
        [
            "screen",
            "--config",
            str(screen_config),
            "--out",
            str(paths["out"]),
            "--top-n",
            "3",
            "--verbose",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Screen: integration_screen" in result.output
    assert "Constructs: 3" in result.output
    assert "Candidate: test_lox_fusion" in result.output
    assert "  Running test_lox_fusion..." in result.output
    assert "Decision: FAIL" in result.output
    assert (paths["out"] / "integration_screen_control_comparison.csv").exists()
    assert (paths["out"] / "integration_screen_summaries.csv").exists()
    assert (paths["out"] / "integration_screen_report.md").exists()
    assert (paths["out"] / "test_lox_baseline_poses.csv").exists()
    assert (paths["out"] / "test_lox_fusion_poses.csv").exists()
    assert (paths["out"] / "test_scrambled_poses.csv").exists()


def test_screen_cli_same_model_dir_controls(
    tmp_path: Path,
    fixture_pdb_path: Path,
) -> None:
    models = tmp_path / "models"
    configs = tmp_path / "configs"
    out = tmp_path / "out"
    models.mkdir()
    configs.mkdir()
    for index in range(1, 4):
        shutil.copyfile(fixture_pdb_path, models / f"model_{index}.pdb")

    candidate_config = configs / "candidate.yaml"
    baseline_config = configs / "baseline.yaml"
    scrambled_config = configs / "scrambled.yaml"
    candidate_config.write_text(
        _config_yaml("test_candidate", "lox_fusion", models)
        + "productive_distance_threshold: 8.0\n",
        encoding="utf-8",
    )
    baseline_config.write_text(
        _config_yaml("test_baseline", "lox_baseline", models),
        encoding="utf-8",
    )
    scrambled_config.write_text(
        _config_yaml("test_scrambled", "control", models),
        encoding="utf-8",
    )
    screen_yaml = configs / "screen.yaml"
    screen_yaml.write_text(
        yaml.safe_dump(
            {
                "screen_id": "same_model_screen",
                "candidate_construct_id": "test_candidate",
                "constructs": [
                    {
                        "construct_id": "test_candidate",
                        "control_type": "candidate",
                        "config_path": str(candidate_config),
                    },
                    {
                        "construct_id": "test_baseline",
                        "control_type": "baseline",
                        "config_path": str(baseline_config),
                    },
                    {
                        "construct_id": "test_scrambled",
                        "control_type": "scrambled",
                        "config_path": str(scrambled_config),
                    },
                ],
                "strict_criteria": {"must_beat_polyK": False},
            }
        ),
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        cli,
        ["screen", "--config", str(screen_yaml), "--out", str(out), "--top-n", "3"],
    )

    assert result.exit_code == 0, result.output
    assert "Screen:" in result.output
    assert "Decision:" in result.output
    report_path = out / "same_model_screen_report.md"
    control_csv = out / "same_model_screen_control_comparison.csv"
    assert report_path.exists()
    assert control_csv.exists()
    report = report_path.read_text(encoding="utf-8")
    assert "## Executive Decision" in report
    assert "## Control Comparison" in report
    control_df = pd.read_csv(control_csv)
    assert len(control_df) >= 2


def test_screen_cli_help_shows_options() -> None:
    result = CliRunner().invoke(cli, ["screen", "--help"])

    assert result.exit_code == 0
    for option in ("--config", "--out", "--top-n", "--verbose"):
        assert option in result.output


def _create_integration_scenario(tmp_path: Path, fixture_pdb_path: Path) -> dict[str, Path]:
    configs = tmp_path / "configs"
    baseline_models = tmp_path / "models" / "lox_baseline"
    fusion_models = tmp_path / "models" / "lox_fusion"
    out = tmp_path / "out"
    configs.mkdir()
    baseline_models.mkdir(parents=True)
    fusion_models.mkdir(parents=True)
    out.mkdir()

    for index, distance in enumerate([6.0, 14.0, 7.0], start=1):
        _write_shifted_fixture(fixture_pdb_path, baseline_models / f"model_{index}.pdb", distance)

    shutil.copyfile(baseline_models / "model_1.pdb", fusion_models / "model_1.pdb")
    _write_shifted_fixture(fixture_pdb_path, fusion_models / "model_2.pdb", 5.0)
    _write_shifted_fixture(fixture_pdb_path, fusion_models / "model_3.pdb", 12.0)

    baseline_config = configs / "lox_baseline.yaml"
    fusion_config = configs / "lox_fusion.yaml"
    baseline_config.write_text(_config_yaml("test_lox_baseline", "lox_baseline", baseline_models), encoding="utf-8")
    fusion_config.write_text(_config_yaml("test_lox_fusion", "lox_fusion", fusion_models), encoding="utf-8")

    return {
        "configs": configs,
        "baseline_config": baseline_config,
        "fusion_config": fusion_config,
        "baseline_models": baseline_models,
        "fusion_models": fusion_models,
        "out": out,
    }


def _write_shifted_fixture(source: Path, destination: Path, distance: float) -> None:
    lines = []
    for line in source.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM") and line[21] == "D" and int(line[22:26]) in ACTIVE_SITE_RESI:
            line = f"{line[:30]}{distance:8.3f}{0.0:8.3f}{0.0:8.3f}{line[54:]}"
        lines.append(line)
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _config_yaml(construct_id: str, construct_type: str, models_dir: Path) -> str:
    return f"""
construct_id: {construct_id}
construct_type: {construct_type}
models_dir: {models_dir.as_posix()}
receptor:
  chains: ["A", "B"]
  target_residues:
    - chain: "A"
      resi: 11
      resn: "LYS"
      atom: "NZ"
ligand:
  chains: ["D"]
lox_domain:
  start: 1
  end: 200
active_site:
  residues:
    - {{label: "His292", resi: 124}}
    - {{label: "His294", resi: 126}}
    - {{label: "His296", resi: 128}}
    - {{label: "Lys320", resi: 152}}
    - {{label: "Tyr355", resi: 187}}
thresholds:
  productive_distance_angstrom: 8.0
  contact_distance_angstrom: 4.0
  best_productive_rank_max: 3
"""


def _audit_config(construct_id: str, models_dir: Path) -> dict:
    return {
        "construct_id": construct_id,
        "models_dir": str(models_dir),
        "ligand_chains": ["D"],
        "receptor_chains": ["A", "B"],
        "expected_n_terminus": None,
    }
