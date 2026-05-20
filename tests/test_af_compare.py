from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml
from click.testing import CliRunner

from loxdockaudit.af_compare import load_af_compare_config, run_af_compare
from loxdockaudit.cli import cli


def test_af_compare_orientation_and_outputs(tmp_path: Path) -> None:
    config_path, out_dir = _write_af_compare_fixture(tmp_path)

    result = run_af_compare(load_af_compare_config(config_path), out_dir)

    assert len(result.hdock_poses) == 2
    assert result.hdock_poses[0].productive is True
    assert result.hdock_poses[0].orientation_angle_deg == pytest.approx(0.0)
    assert result.af_geometry.is_productive is True
    assert result.af_geometry.is_fully_productive is True
    assert result.convergence_category == "strong convergence"

    summary = pd.read_csv(out_dir / "synthetic_af_compare_af_compare.csv")
    assert summary.loc[0, "hdock_productive_count"] == 1
    assert summary.loc[0, "af_fully_productive"] is True or bool(
        summary.loc[0, "af_fully_productive"]
    )
    assert (out_dir / "synthetic_af_compare_interface_overlap.csv").exists()
    assert (out_dir / "synthetic_af_compare_contact_persistence.csv").exists()
    report = (out_dir / "synthetic_af_compare_report.md").read_text(
        encoding="utf-8"
    )
    assert "computational agreement" in report.lower()
    assert "does not establish enzymatic activity" in report
    for forbidden in ("validated", "proved", "confirmed mechanism", "demonstrated catalysis"):
        assert forbidden not in report.lower()


def test_af_compare_cli(tmp_path: Path) -> None:
    config_path, out_dir = _write_af_compare_fixture(tmp_path)

    result = CliRunner().invoke(
        cli,
        ["af-compare", "--config", str(config_path), "--out", str(out_dir)],
    )

    assert result.exit_code == 0, result.output
    assert "Comparison: synthetic_af_compare" in result.output
    assert "Convergence category: strong convergence" in result.output
    assert "does not establish enzymatic activity" in result.output


def _write_af_compare_fixture(tmp_path: Path) -> tuple[Path, Path]:
    hdock_dir = tmp_path / "hdock"
    hdock_dir.mkdir()
    out_dir = tmp_path / "out"
    (hdock_dir / "model_1.pdb").write_text(
        _pose_text(active_site_x=0.0),
        encoding="utf-8",
    )
    (hdock_dir / "model_2.pdb").write_text(
        _pose_text(active_site_x=20.0),
        encoding="utf-8",
    )
    af_path = tmp_path / "rank_001_af.pdb"
    af_path.write_text(_pose_text(active_site_x=0.0), encoding="utf-8")
    config_path = tmp_path / "af_compare.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "comparison_id": "synthetic_af_compare",
                "hdock_models_dir": str(hdock_dir),
                "alphafold_model_path": str(af_path),
                "target_residue": {
                    "chain": "A",
                    "resi": 11,
                    "resn": "LYS",
                    "atom": "NZ",
                },
                "active_site_chain": "D",
                "active_site_residues": [100],
                "productive_distance_threshold": 8.0,
                "contact_distance_threshold": 5.0,
                "orientation_enabled": True,
                "orientation_threshold_deg": 90.0,
            }
        ),
        encoding="utf-8",
    )
    return config_path, out_dir


def _pose_text(active_site_x: float) -> str:
    return "\n".join(
        [
            "ATOM      1  CB  LYS A  11       3.000   0.000   0.000  1.00 80.00           C  ",
            "ATOM      2  NZ  LYS A  11       4.000   0.000   0.000  1.00 80.00           N  ",
            f"ATOM      3  CA  HIS D 100      {active_site_x:6.3f}   0.000   0.000  1.00 80.00           C  ",
            "TER",
            "END",
        ]
    )
