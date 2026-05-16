from __future__ import annotations

import shutil
from pathlib import Path

from click.testing import CliRunner

from loxdockaudit.cli import cli


def test_v03_run_with_structural_qc_configured(
    tmp_path: Path,
    fixture_pdb_path: Path,
) -> None:
    paths = _write_cli_fixture(tmp_path, fixture_pdb_path, "v03_qc", structural_qc=True)

    result = CliRunner().invoke(
        cli,
        ["run", "--config", str(paths["config"]), "--out", str(paths["out"])],
    )

    assert result.exit_code == 0, result.output
    assert "Structural QC:" in result.output
    assert (paths["out"] / "v03_qc_structural_qc.md").exists()
    assert "## Structural QC" in (paths["out"] / "v03_qc_report.md").read_text(encoding="utf-8")


def test_v03_run_without_structural_qc(
    tmp_path: Path,
    fixture_pdb_path: Path,
) -> None:
    paths = _write_cli_fixture(tmp_path, fixture_pdb_path, "v03_no_qc", structural_qc=False)

    result = CliRunner().invoke(
        cli,
        ["run", "--config", str(paths["config"]), "--out", str(paths["out"])],
    )

    assert result.exit_code == 0, result.output
    assert "Structural QC:" not in result.output
    report = (paths["out"] / "v03_no_qc_report.md").read_text(encoding="utf-8")
    assert "Structural QC not configured." in report


def test_v03_strict_qc_with_failing_geometry_exits(
    tmp_path: Path,
    fixture_pdb_path: Path,
) -> None:
    paths = _write_cli_fixture(
        tmp_path,
        fixture_pdb_path,
        "v03_bad_qc",
        structural_qc=True,
        bad_his=True,
    )

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "--config",
            str(paths["config"]),
            "--out",
            str(paths["out"]),
            "--strict-qc",
        ],
    )

    assert result.exit_code == 1
    assert "ERROR" in result.output


def _write_cli_fixture(
    tmp_path: Path,
    fixture_pdb_path: Path,
    construct_id: str,
    *,
    structural_qc: bool,
    bad_his: bool = False,
) -> dict[str, Path]:
    models = tmp_path / f"{construct_id}_models"
    out = tmp_path / f"{construct_id}_out"
    config = tmp_path / f"{construct_id}.yaml"
    models.mkdir()
    shutil.copyfile(fixture_pdb_path, models / "model_1.pdb")

    qc_block = ""
    if structural_qc:
        his_resi = "[999, 998, 997]" if bad_his else "[124, 126, 128]"
        qc_block = f"""
structural_qc:
  his_resi: {his_resi}
  lys_resi: 152
  tyr_resi: 187
  disulfide_pairs: []
"""

    config.write_text(
        f"""
construct_id: {construct_id}
construct_type: lox_baseline
models_dir: {models.as_posix()}
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
{qc_block}
""",
        encoding="utf-8",
    )
    return {"models": models, "out": out, "config": config}
