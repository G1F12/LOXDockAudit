from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from loxdockaudit.config import load_construct_config, load_screen_config
from loxdockaudit.models import ScreenResult
from loxdockaudit.screen import pearson_correlation, run_screen, run_single_construct
from tests.conftest import write_pose


def test_pearson_correlation() -> None:
    assert pearson_correlation([1, 2, 3], [1, 2, 3]) == pytest.approx(1.0)
    assert pearson_correlation([1, 2, 3], [3, 2, 1]) == pytest.approx(-1.0)
    assert pearson_correlation([1, 1, 1], [1, 2, 3]) == 0.0
    assert pearson_correlation([1, 2], [1, 2]) == 0.0


def test_run_single_construct(tmp_path: Path) -> None:
    models = tmp_path / "models"
    models.mkdir()
    write_pose(models / "model_1.pdb", 6.0)
    write_pose(models / "model_2.pdb", 14.0)
    write_pose(models / "model_3.pdb", 7.0)
    config_path = _write_construct_config(tmp_path, "candidate", models)

    summary, poses = run_single_construct(load_construct_config(str(config_path)), top_n=3)

    assert summary.construct_id == "candidate"
    assert summary.productive_count == 2
    assert summary.total_poses == 3
    assert summary.best_distance == pytest.approx(6.0)
    assert summary.best_productive_rank == 1
    assert len(poses) == 3


def test_run_single_construct_skips_bad_pose(tmp_path: Path) -> None:
    models = tmp_path / "models"
    models.mkdir()
    write_pose(models / "model_1.pdb", 6.0)
    (models / "model_2.pdb").write_bytes(b"")
    config_path = _write_construct_config(tmp_path, "candidate", models)

    with pytest.warns(RuntimeWarning):
        summary, poses = run_single_construct(load_construct_config(str(config_path)), top_n=2)

    assert summary.total_poses == 1
    assert len(poses) == 1


def test_run_single_construct_raises_when_no_pose_processed(tmp_path: Path) -> None:
    models = tmp_path / "models"
    models.mkdir()
    (models / "model_1.pdb").write_bytes(b"")
    config_path = _write_construct_config(tmp_path, "candidate", models)

    with pytest.warns(RuntimeWarning), pytest.raises(RuntimeError, match="candidate"):
        run_single_construct(load_construct_config(str(config_path)), top_n=1)


def test_run_single_construct_empty_models_dir_raises(tmp_path: Path) -> None:
    models = tmp_path / "models"
    models.mkdir()
    config_path = _write_construct_config(tmp_path, "candidate", models)

    with pytest.raises(RuntimeError, match="candidate"):
        run_single_construct(load_construct_config(str(config_path)), top_n=10)


def test_run_screen(tmp_path: Path) -> None:
    baseline_models = tmp_path / "baseline_models"
    candidate_models = tmp_path / "candidate_models"
    scrambled_models = tmp_path / "scrambled_models"
    for path in (baseline_models, candidate_models, scrambled_models):
        path.mkdir()

    for index, distance in enumerate([6.0, 14.0, 12.0], start=1):
        write_pose(baseline_models / f"model_{index}.pdb", distance)
    for index, distance in enumerate([5.0, 7.0, 9.0], start=1):
        write_pose(candidate_models / f"model_{index}.pdb", distance)
    for index, distance in enumerate([4.0, 5.0, 6.0], start=1):
        write_pose(scrambled_models / f"model_{index}.pdb", distance)

    baseline_config = _write_construct_config(tmp_path, "baseline", baseline_models)
    candidate_config = _write_construct_config(tmp_path, "candidate", candidate_models)
    scrambled_config = _write_construct_config(tmp_path, "scrambled", scrambled_models)
    screen_config_path = tmp_path / "screen.yaml"
    screen_config_path.write_text(
        yaml.safe_dump(
            {
                "screen_id": "screen",
                "candidate_construct_id": "candidate",
                "constructs": [
                    {
                        "construct_id": "baseline",
                        "control_type": "baseline",
                        "config_path": str(baseline_config),
                    },
                    {
                        "construct_id": "candidate",
                        "control_type": "candidate",
                        "config_path": str(candidate_config),
                    },
                    {
                        "construct_id": "scrambled",
                        "control_type": "scrambled",
                        "config_path": str(scrambled_config),
                    },
                ],
                "strict_criteria": {
                    "must_beat_polyK": False,
                },
            }
        ),
        encoding="utf-8",
    )

    result = run_screen(load_screen_config(str(screen_config_path)), top_n=3)

    assert isinstance(result, ScreenResult)
    assert result.screen_id == "screen"
    assert result.candidate_summary.construct_id == "candidate"
    assert result.baseline_summary is not None
    assert result.candidate_summary.strict_pass is result.overall_pass
    assert result.criteria_results["beats_scrambled"] is False
    assert result.overall_pass is False
    assert any(comparison.control_type == "baseline" for comparison in result.control_comparisons)
    assert "does not advance" in result.decision_text


def test_run_screen_adds_screen_id_to_runtime_error(tmp_path: Path) -> None:
    baseline_models = tmp_path / "baseline_models"
    candidate_models = tmp_path / "candidate_models"
    baseline_models.mkdir()
    candidate_models.mkdir()
    write_pose(baseline_models / "model_1.pdb", 6.0)

    baseline_config = _write_construct_config(tmp_path, "baseline", baseline_models)
    candidate_config = _write_construct_config(tmp_path, "candidate", candidate_models)
    screen_config_path = tmp_path / "screen.yaml"
    screen_config_path.write_text(
        yaml.safe_dump(
            {
                "screen_id": "screen_error",
                "candidate_construct_id": "candidate",
                "constructs": [
                    {
                        "construct_id": "baseline",
                        "control_type": "baseline",
                        "config_path": str(baseline_config),
                    },
                    {
                        "construct_id": "candidate",
                        "control_type": "candidate",
                        "config_path": str(candidate_config),
                    },
                ],
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="screen_error"):
        run_screen(load_screen_config(str(screen_config_path)), top_n=3)


def _write_construct_config(tmp_path: Path, construct_id: str, models_dir: Path) -> Path:
    path = tmp_path / f"{construct_id}.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "construct_id": construct_id,
                "construct_type": "lox_baseline",
                "models_dir": str(models_dir),
                "receptor": {
                    "chains": ["A", "B"],
                    "target_residues": [
                        {"chain": "A", "resi": 11, "resn": "LYS", "atom": "NZ"}
                    ],
                },
                "ligand": {"chains": ["D"]},
                "active_site": {
                    "residues": [
                        {"label": "His292", "resi": 124},
                        {"label": "His294", "resi": 126},
                        {"label": "His296", "resi": 128},
                        {"label": "Lys320", "resi": 152},
                        {"label": "Tyr355", "resi": 187},
                    ]
                },
                "lox_domain": {"start": 1, "end": 200},
            }
        ),
        encoding="utf-8",
    )
    return path
