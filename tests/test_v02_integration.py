from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml
from click.testing import CliRunner

from loxdockaudit.cli import cli
from loxdockaudit.config import load_screen_config
from loxdockaudit.screen import run_screen


ACTIVE_SITE_RESI = {124, 126, 128, 152, 187}


def test_v02_full_screen_pipeline(tmp_path: Path, fixture_pdb_path: Path) -> None:
    paths = _build_v02_screen(tmp_path, fixture_pdb_path)
    screen_config = load_screen_config(str(paths["screen_yaml"]))

    screen_result = run_screen(screen_config, top_n=3)

    assert screen_result.overall_pass is True
    assert screen_result.candidate_summary.productive_count == 2
    assert screen_result.baseline_summary is not None
    assert screen_result.baseline_summary.productive_count == 1
    assert screen_result.criteria_results["pass_count_vs_baseline"] is True
    assert screen_result.criteria_results["beats_scrambled"] is True
    assert len(screen_result.control_comparisons) == 3
    assert "advances" in screen_result.decision_text

    result = CliRunner().invoke(
        cli,
        [
            "screen",
            "--config",
            str(paths["screen_yaml"]),
            "--out",
            str(paths["out"]),
            "--top-n",
            "3",
        ],
    )

    assert result.exit_code == 0, result.output
    report_path = paths["out"] / "v02_screen_report.md"
    comparison_path = paths["out"] / "v02_screen_control_comparison.csv"
    assert report_path.exists()
    assert "PASS" in report_path.read_text(encoding="utf-8")
    assert len(pd.read_csv(comparison_path)) == 3


def _build_v02_screen(tmp_path: Path, fixture_pdb_path: Path) -> dict[str, Path]:
    configs = tmp_path / "configs"
    models = tmp_path / "models"
    out = tmp_path / "out"
    configs.mkdir()
    models.mkdir()
    out.mkdir()

    specs = {
        "candidate": ("candidate", [6.0, 14.0, 7.0]),
        "baseline": ("baseline", [8.0, 14.0, 14.0]),
        "scrambled": ("scrambled", [12.0, 15.0, 13.0]),
        "polyK": ("polyK", [8.0, 14.0, 14.0]),
    }
    config_paths: dict[str, Path] = {}
    for construct_id, (_, distances) in specs.items():
        model_dir = models / construct_id
        model_dir.mkdir()
        for index, distance in enumerate(distances, start=1):
            _write_shifted_fixture(
                fixture_pdb_path,
                model_dir / f"model_{index}.pdb",
                distance,
            )
        config_path = configs / f"{construct_id}.yaml"
        config_path.write_text(_construct_yaml(construct_id, model_dir), encoding="utf-8")
        config_paths[construct_id] = config_path

    screen_yaml = configs / "screen.yaml"
    screen_yaml.write_text(
        yaml.safe_dump(
            {
                "screen_id": "v02_screen",
                "candidate_construct_id": "candidate",
                "constructs": [
                    {
                        "construct_id": "candidate",
                        "control_type": "candidate",
                        "config_path": str(config_paths["candidate"]),
                    },
                    {
                        "construct_id": "baseline",
                        "control_type": "baseline",
                        "config_path": str(config_paths["baseline"]),
                    },
                    {
                        "construct_id": "scrambled",
                        "control_type": "scrambled",
                        "config_path": str(config_paths["scrambled"]),
                    },
                    {
                        "construct_id": "polyK",
                        "control_type": "polyK",
                        "config_path": str(config_paths["polyK"]),
                    },
                ],
                "strict_criteria": {
                    "pass_count_greater_than_baseline": True,
                    "best_distance_lte_baseline": True,
                    "best_productive_rank_max": 3,
                    "must_beat_scrambled": True,
                    "must_beat_polyK": True,
                },
            }
        ),
        encoding="utf-8",
    )

    return {"screen_yaml": screen_yaml, "out": out, "configs": configs, "models": models}


def _write_shifted_fixture(source: Path, destination: Path, distance: float) -> None:
    lines = []
    for line in source.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM") and line[21] == "D" and int(line[22:26]) in ACTIVE_SITE_RESI:
            line = f"{line[:30]}{distance:8.3f}{0.0:8.3f}{0.0:8.3f}{line[54:]}"
        lines.append(line)
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _construct_yaml(construct_id: str, models_dir: Path) -> str:
    return f"""
construct_id: {construct_id}
construct_type: lox_baseline
models_dir: {models_dir.as_posix()}
productive_distance_threshold: 8.0
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
"""
