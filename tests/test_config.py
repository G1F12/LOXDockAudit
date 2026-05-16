from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from loxdockaudit.config import load_construct_config, load_screen_config
from loxdockaudit.models import ConstructConfig, ScreenConfig


VALID_YAML = """
construct_id: baseline
construct_type: lox_baseline
models_dir: models/baseline
receptor:
  chains: [A, B]
  target_residues:
    - chain: A
      resi: 11
      resn: LYS
ligand:
  chains: [D]
active_site:
  residues:
    - label: His292
      resi: 124
lox_domain:
  start: 1
  end: 200
"""


def test_load_construct_config_valid(tmp_path: Path) -> None:
    path = tmp_path / "valid.yaml"
    path.write_text(VALID_YAML, encoding="utf-8")
    config = load_construct_config(str(path))
    assert isinstance(config, ConstructConfig)
    assert config.construct_id == "baseline"
    assert config.receptor_chains == ["A", "B"]
    assert config.ligand_chains == ["D"]


def test_load_construct_config_missing_construct_id(tmp_path: Path) -> None:
    path = tmp_path / "missing.yaml"
    path.write_text(VALID_YAML.replace("construct_id: baseline\n", ""), encoding="utf-8")
    with pytest.raises(ValueError, match="construct_id"):
        load_construct_config(str(path))


def test_load_construct_config_missing_active_site(tmp_path: Path) -> None:
    path = tmp_path / "missing.yaml"
    path.write_text(VALID_YAML.replace("active_site:\n  residues:\n    - label: His292\n      resi: 124\n", ""), encoding="utf-8")
    with pytest.raises(ValueError, match="active_site"):
        load_construct_config(str(path))


def test_load_construct_config_no_cbd_domain(tmp_path: Path) -> None:
    path = tmp_path / "valid.yaml"
    path.write_text(VALID_YAML, encoding="utf-8")
    assert load_construct_config(str(path)).cbd_domain is None


def test_load_construct_config_is_inactive_true(tmp_path: Path) -> None:
    path = tmp_path / "inactive.yaml"
    path.write_text(VALID_YAML + "is_inactive: true\n", encoding="utf-8")
    assert load_construct_config(str(path)).is_inactive is True


def test_load_screen_config_three_entries(tmp_path: Path) -> None:
    path = tmp_path / "screen.yaml"
    construct_config = tmp_path / "construct.yaml"
    construct_config.write_text(VALID_YAML, encoding="utf-8")
    constructs = []
    for construct_id, control_type in (
        ("baseline", "baseline"),
        ("candidate", "candidate"),
        ("scrambled", "scrambled"),
    ):
        constructs.append(
            {
                "construct_id": construct_id,
                "control_type": control_type,
                "config_path": "construct.yaml",
            }
        )
    path.write_text(
        yaml.safe_dump(
            {
                "screen_id": "screen",
                "candidate_construct_id": "candidate",
                "constructs": constructs,
            }
        ),
        encoding="utf-8",
    )
    screen = load_screen_config(str(path))
    assert isinstance(screen, ScreenConfig)
    assert screen.screen_id == "screen"
    assert len(screen.constructs) == 3
    assert screen.strict_criteria.best_productive_rank_max == 3


def test_load_screen_config_candidate_must_appear(tmp_path: Path) -> None:
    construct_config = tmp_path / "construct.yaml"
    construct_config.write_text(VALID_YAML, encoding="utf-8")
    path = tmp_path / "screen.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "screen_id": "screen",
                "candidate_construct_id": "missing",
                "constructs": [
                    {
                        "construct_id": "baseline",
                        "control_type": "baseline",
                        "config_path": "construct.yaml",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="candidate_construct_id"):
        load_screen_config(str(path))


def test_load_screen_config_requires_baseline(tmp_path: Path) -> None:
    construct_config = tmp_path / "construct.yaml"
    construct_config.write_text(VALID_YAML, encoding="utf-8")
    path = tmp_path / "screen.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "screen_id": "screen",
                "candidate_construct_id": "candidate",
                "constructs": [
                    {
                        "construct_id": "candidate",
                        "control_type": "candidate",
                        "config_path": "construct.yaml",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="baseline"):
        load_screen_config(str(path))


def test_load_screen_config_requires_readable_config_path(tmp_path: Path) -> None:
    path = tmp_path / "screen.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "screen_id": "screen",
                "candidate_construct_id": "candidate",
                "constructs": [
                    {
                        "construct_id": "baseline",
                        "control_type": "baseline",
                        "config_path": "missing.yaml",
                    },
                    {
                        "construct_id": "candidate",
                        "control_type": "candidate",
                        "config_path": "missing.yaml",
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="config_path"):
        load_screen_config(str(path))


def test_productive_distance_threshold_default(tmp_path: Path) -> None:
    path = tmp_path / "valid.yaml"
    path.write_text(VALID_YAML, encoding="utf-8")
    assert load_construct_config(str(path)).productive_distance_threshold == 8.0
