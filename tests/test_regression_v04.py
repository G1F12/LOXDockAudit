from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml

from loxdockaudit.config import load_construct_config
from loxdockaudit.inactive_control import (
    load_inactive_control_config,
    run_inactive_control_analysis,
)
from loxdockaudit.screen import run_single_construct


ROOT = Path(__file__).resolve().parents[1]
ACTIVE_ID = "r5_real_lox169"
INACTIVE_ID = "r5_inactive_h292a_h294a_h296a"
COMPARISON_ID = "r5_lox169_active_vs_inactive"


def test_v03_metrics_unchanged_with_orientation_disabled(tmp_path: Path) -> None:
    result = _run_inactive_control(
        tmp_path,
        config_path=ROOT / "configs" / "r5_inactive_control.yaml",
    )

    active = result["active_summary"]
    inactive = result["inactive_summary"]

    assert active.productive_count == 2
    assert active.best_distance == pytest.approx(6.159, abs=0.001)
    assert active.best_productive_rank == 3
    assert inactive.productive_count == 1
    assert inactive.best_distance == pytest.approx(6.421, abs=0.001)
    assert inactive.best_productive_rank == 7


def test_orientation_angle_in_valid_range(tmp_path: Path) -> None:
    config_paths = _write_round5_configs(tmp_path, orientation_enabled=True)
    active_config = load_construct_config(str(config_paths["active"]))
    _, poses = run_single_construct(active_config, top_n=10)

    angles = [pose.orientation_angle_deg for pose in poses]

    assert len(angles) == 10
    for angle in angles:
        assert angle is not None
        assert pd.notna(angle)
        assert 0.0 <= angle <= 180.0


def test_fully_productive_subset_of_productive(tmp_path: Path) -> None:
    config_paths = _write_round5_configs(tmp_path, orientation_enabled=True)
    result = _run_inactive_control(tmp_path, config_path=config_paths["comparison"])

    for pose_df_name in (
        f"{ACTIVE_ID}_poses.csv",
        f"{INACTIVE_ID}_poses.csv",
    ):
        pose_df = pd.read_csv(result["out_dir"] / pose_df_name)
        for _, row in pose_df.iterrows():
            if bool(row["fully_productive"]):
                assert bool(row["productive"]) is True


def test_orientation_columns_in_output_csv(tmp_path: Path) -> None:
    config_paths = _write_round5_configs(tmp_path, orientation_enabled=True)
    result = _run_inactive_control(tmp_path, config_path=config_paths["comparison"])

    active_pose_df = pd.read_csv(result["out_dir"] / f"{ACTIVE_ID}_poses.csv")

    for column in (
        "orientation_angle_deg",
        "orientation_productive",
        "fully_productive",
    ):
        assert column in active_pose_df.columns


def test_scatter_svg_generated(tmp_path: Path) -> None:
    config_paths = _write_round5_configs(tmp_path, orientation_enabled=True)
    result = _run_inactive_control(tmp_path, config_path=config_paths["comparison"])

    for construct_id in (ACTIVE_ID, INACTIVE_ID):
        svg_path = result["out_dir"] / f"{construct_id}_geometry_scatter.svg"
        assert svg_path.exists()
        assert svg_path.stat().st_size > 0


def test_natural_sort_still_correct(tmp_path: Path) -> None:
    result = _run_inactive_control(
        tmp_path,
        config_path=ROOT / "configs" / "r5_inactive_control.yaml",
    )

    assert result["active_summary"].best_productive_rank == 3


def _run_inactive_control(
    tmp_path: Path,
    config_path: Path,
) -> dict:
    out_dir = tmp_path / "out"
    result = run_inactive_control_analysis(
        load_inactive_control_config(str(config_path)),
        out_dir,
        top_n=10,
    )
    result["out_dir"] = out_dir
    return result


def _write_round5_configs(
    tmp_path: Path,
    orientation_enabled: bool,
) -> dict[str, Path]:
    active_config = _load_yaml(ROOT / "examples" / "real_hdock_round5" / "config.yaml")
    inactive_config = _load_yaml(ROOT / "configs" / "r5_inactive.yaml")
    comparison_config = _load_yaml(ROOT / "configs" / "r5_inactive_control.yaml")

    orientation = {
        "enabled": orientation_enabled,
        "threshold_deg": 90.0,
        "sidechain_atoms": ["CB", "NZ"],
        "activesite_centroid_atoms": ["CA"],
    }
    active_config["orientation"] = orientation
    inactive_config["orientation"] = orientation

    active_path = tmp_path / "active.yaml"
    inactive_path = tmp_path / "inactive.yaml"
    comparison_path = tmp_path / "comparison.yaml"

    _write_yaml(active_path, active_config)
    _write_yaml(inactive_path, inactive_config)

    comparison_config["active"]["config_path"] = str(active_path)
    comparison_config["inactive"]["config_path"] = str(inactive_path)
    _write_yaml(comparison_path, comparison_config)

    return {
        "active": active_path,
        "inactive": inactive_path,
        "comparison": comparison_path,
    }


def _load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    assert isinstance(data, dict)
    return data


def _write_yaml(path: Path, data: dict) -> None:
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
