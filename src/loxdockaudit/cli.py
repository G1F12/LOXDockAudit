"""Command line interface for LOXDockAudit."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import click
import numpy as np

from loxdockaudit.config import load_construct_config, load_screen_config
from loxdockaudit.contacts import analyze_pose_contacts, compute_contact_frequencies
from loxdockaudit.distances import active_site_to_target_distance
from loxdockaudit.input_audit import audit_screen, format_audit_report
from loxdockaudit.inactive_control import (
    load_inactive_control_config,
    run_inactive_control_analysis,
)
from loxdockaudit.fold_qc import format_qc_report, run_structural_qc
from loxdockaudit.models import ConstructConfig, ConstructSummary, PoseMetrics, StructuralQCResult
from loxdockaudit.pdb_parser import load_structure
from loxdockaudit.reporting import (
    control_comparisons_to_dataframe,
    generate_markdown_report,
    generate_screen_markdown_report,
    poses_to_dataframe,
    save_outputs,
    save_screen_outputs,
    summaries_to_dataframe,
)
from loxdockaudit.screen import run_screen, run_single_construct
from loxdockaudit.sorting import natural_sort_key


@click.group()
def cli() -> None:
    """LOXDockAudit command line tools."""
    _configure_stdio()
    pass


@cli.command()
@click.option("--config", required=True, help="Path to construct YAML config")
@click.option("--out", required=True, help="Output directory")
@click.option("--top-n", default=10, help="Number of top poses to analyze")
@click.option("--strict-qc", is_flag=True, help="Exit if structural QC fails")
@click.option("--verbose", is_flag=True)
def run(config: str, out: str, top_n: int, strict_qc: bool, verbose: bool) -> None:
    """Run a LOXDockAudit analysis."""
    try:
        construct_config = load_construct_config(config)
        pdb_paths = sorted(
            Path(construct_config.models_dir).glob("*.pdb"),
            key=lambda path: natural_sort_key(str(path)),
        )
        qc_result = _run_configured_structural_qc(construct_config, out)
        if qc_result is not None:
            _print_qc_summary(qc_result)
            if qc_result.fold_corrupted and strict_qc:
                click.echo("ERROR: Structural QC failed. Use --no-strict-qc to override.")
                sys.exit(1)
            if qc_result.fold_corrupted:
                click.echo("WARNING: Structural QC failed. Docking results may be unreliable.")

        audit_config = _audit_config_from_construct(construct_config)
        audit_result = audit_screen([audit_config])
        audit_report = format_audit_report(audit_result)

        if verbose:
            click.echo(audit_report.rstrip())

        if not audit_result["screen_valid"]:
            if not verbose:
                click.echo(audit_report.rstrip())
            click.echo("ERROR: input audit failed")
            sys.exit(1)

        pose_metrics: list[PoseMetrics] = []
        contact_results: list[dict[str, Any]] = []
        target_residue = construct_config.target_residues[0]
        active_site_resi = [
            int(residue["resi"])
            for residue in construct_config.active_site.residues
            if "resi" in residue
        ]

        for rank, pdb_path in enumerate(pdb_paths[:top_n], start=1):
            structure = load_structure(str(pdb_path), structure_id=f"pose_{rank}")
            distance_result = active_site_to_target_distance(
                structure=structure,
                active_site_resi=active_site_resi,
                active_site_chains=construct_config.ligand_chains,
                target_resi=target_residue.resi,
                target_chain=target_residue.chain,
                threshold=construct_config.productive_distance_threshold,
                target_atom_names=[target_residue.atom],
            )
            contact_result = analyze_pose_contacts(
                structure=structure,
                receptor_chains=construct_config.receptor_chains,
                ligand_chains=construct_config.ligand_chains,
                lox_domain_resi_range=(
                    construct_config.lox_domain.start,
                    construct_config.lox_domain.end,
                ),
                cbd_domain_resi_range=_domain_tuple(construct_config.cbd_domain),
                contact_cutoff=construct_config.contact_distance_threshold,
            )
            contact_results.append(contact_result)

            pose_metrics.append(
                PoseMetrics(
                    rank=rank,
                    model_path=str(pdb_path),
                    active_site_to_target_distance=distance_result[
                        "distance_angstrom"
                    ],
                    productive=distance_result["productive"],
                    lox_contacts=contact_result["lox_contacts"],
                    cbd_contacts=contact_result["cbd_contacts"],
                    closest_active_site_resi=distance_result[
                        "closest_active_site_resi"
                    ],
                    warning=distance_result["warning"],
                )
            )

        summary = _build_summary(construct_config, pose_metrics, contact_results)
        pose_df = poses_to_dataframe(pose_metrics)
        summary_df = summaries_to_dataframe([summary])
        markdown_report = generate_markdown_report(
            construct_config=construct_config,
            pose_df=pose_df,
            summary=summary,
            audit_report=audit_report,
            qc_result=qc_result,
        )
        save_outputs(
            out_dir=out,
            construct_id=construct_config.construct_id,
            pose_df=pose_df,
            summary_df=summary_df,
            markdown_report=markdown_report,
        )

        click.echo(f"Construct: {construct_config.construct_id}")
        click.echo(f"Productive poses: {summary.productive_count}/{summary.total_poses}")
        click.echo(f"Best distance: {summary.best_distance:.2f} Å")
        click.echo(f"Best productive rank: {summary.best_productive_rank}")
        click.echo(f"Output saved to: {out}")
    except Exception as exc:
        click.echo(f"ERROR: {exc}")
        sys.exit(1)


@cli.command("check-config")
@click.option("--config", required=True, help="Path to construct YAML config")
def check_config(config: str) -> None:
    """Load and validate a LOXDockAudit construct YAML config."""
    try:
        construct_config = load_construct_config(config)
        click.echo("Config valid")
        click.echo(f"Construct: {construct_config.construct_id}")
        click.echo(f"Type: {construct_config.construct_type}")
        click.echo(f"Models dir: {construct_config.models_dir}")
        click.echo(f"Receptor chains: {', '.join(construct_config.receptor_chains)}")
        click.echo(f"Ligand chains: {', '.join(construct_config.ligand_chains)}")
        click.echo(f"Active-site residues: {len(construct_config.active_site.residues)}")
        click.echo(f"Target residues: {len(construct_config.target_residues)}")
    except Exception as exc:
        click.echo(f"ERROR: {exc}")
        sys.exit(1)


@click.command("screen")
@click.option("--config", required=True, help="Path to screen-level YAML config")
@click.option("--out", required=True, help="Output directory")
@click.option("--top-n", default=10, help="Number of top poses per construct")
@click.option("--verbose", is_flag=True)
def screen(config: str, out: str, top_n: int, verbose: bool) -> None:
    """Run a full control-aware screen comparison."""
    try:
        screen_config = load_screen_config(config)
        click.echo(f"Screen: {screen_config.screen_id}")
        click.echo(f"Constructs: {len(screen_config.constructs)}")
        click.echo(f"Candidate: {screen_config.candidate_construct_id}")

        if verbose:
            for entry in screen_config.constructs:
                click.echo(f"  Running {entry.construct_id}...")

        screen_result = run_screen(screen_config, top_n=top_n, verbose=verbose)

        pose_dfs = {}
        summaries_by_construct = {
            screen_result.candidate_summary.construct_id: screen_result.candidate_summary
        }
        if screen_result.baseline_summary is not None:
            summaries_by_construct[
                screen_result.baseline_summary.construct_id
            ] = screen_result.baseline_summary

        for entry in screen_config.constructs:
            construct_config = load_construct_config(entry.config_path)
            qc_result = _run_configured_structural_qc(construct_config, out)
            if qc_result is not None:
                click.echo(
                    f"Structural QC {entry.construct_id}: "
                    f"{'PASS' if qc_result.fold_qc_pass else 'FOLD CORRUPTED'}"
                )
            summary, pose_metrics = run_single_construct(
                construct_config,
                top_n=top_n,
                verbose=verbose,
            )
            if entry.construct_id in summaries_by_construct:
                summary = summaries_by_construct[entry.construct_id]
            else:
                summaries_by_construct[entry.construct_id] = summary
            pose_dfs[entry.construct_id] = poses_to_dataframe(pose_metrics)

        control_df = control_comparisons_to_dataframe(
            screen_result.control_comparisons
        )
        all_summaries_df = summaries_to_dataframe(
            [
                summaries_by_construct[entry.construct_id]
                for entry in screen_config.constructs
                if entry.construct_id in summaries_by_construct
            ]
        )
        markdown_report = generate_screen_markdown_report(screen_result, pose_dfs)

        save_screen_outputs(
            out,
            screen_config.screen_id,
            screen_result,
            control_df,
            all_summaries_df,
            markdown_report,
        )

        output_dir = Path(out)
        output_dir.mkdir(parents=True, exist_ok=True)
        for construct_id, pose_df in pose_dfs.items():
            pose_df.to_csv(output_dir / f"{construct_id}_poses.csv", index=False)

        _print_screen_summary(screen_result, out)
    except Exception as exc:
        click.echo(f"ERROR: {exc}")
        sys.exit(1)


@click.command("inactive-control")
@click.option("--config", required=True, help="Path to inactive-control YAML config")
@click.option("--out", required=True, help="Output directory")
@click.option("--top-n", default=None, type=int, help="Override top-N poses to analyze")
def inactive_control(config: str, out: str, top_n: int | None) -> None:
    """Compare active LOX against a catalytically inactive control."""
    try:
        comparison_config = load_inactive_control_config(config)
        result = run_inactive_control_analysis(
            comparison_config=comparison_config,
            out_dir=out,
            top_n=top_n,
        )
        active_summary = result["active_summary"]
        inactive_summary = result["inactive_summary"]

        click.echo(f"Comparison: {comparison_config.comparison_id}")
        click.echo(
            "Active productive poses: "
            f"{active_summary.productive_count}/{active_summary.total_poses}"
        )
        click.echo(
            "Inactive productive poses: "
            f"{inactive_summary.productive_count}/{inactive_summary.total_poses}"
        )
        click.echo(f"Active best distance: {active_summary.best_distance:.2f} A")
        click.echo(f"Inactive best distance: {inactive_summary.best_distance:.2f} A")
        click.echo(f"Interpretation: {result['interpretation']}")
        click.echo(f"Output saved to: {out}")
    except Exception as exc:
        click.echo(f"ERROR: {exc}")
        sys.exit(1)


def _audit_config_from_construct(construct_config: ConstructConfig) -> dict[str, Any]:
    return {
        "construct_id": construct_config.construct_id,
        "models_dir": construct_config.models_dir,
        "ligand_chains": construct_config.ligand_chains,
        "receptor_chains": construct_config.receptor_chains,
        "expected_n_terminus": None,
    }


def _run_configured_structural_qc(
    construct_config: ConstructConfig,
    out_dir: str,
) -> StructuralQCResult | None:
    qc_config = construct_config.structural_qc
    if qc_config is None:
        return None

    model_path = _first_model_path(construct_config.models_dir)
    active_site_resi = [
        int(residue["resi"])
        for residue in construct_config.active_site.residues
        if "resi" in residue
    ]
    qc_result = run_structural_qc(
        pdb_path=str(model_path),
        construct_id=construct_config.construct_id,
        lox_chain=construct_config.ligand_chains[0],
        his_resi=qc_config.his_resi,
        lys_resi=qc_config.lys_resi,
        tyr_resi=qc_config.tyr_resi,
        disulfide_pairs=qc_config.disulfide_pairs,
        active_site_resi=active_site_resi,
        alphafold_pdb_path=qc_config.alphafold_pdb_path,
        pae_json_path=qc_config.pae_json_path,
        domain_a_resi=qc_config.domain_a_resi,
        domain_b_resi=qc_config.domain_b_resi,
        plddt_threshold=qc_config.plddt_threshold,
        his_max_ca_distance=qc_config.his_max_ca_distance,
        lys_tyr_max_cb_distance=qc_config.lys_tyr_max_cb_distance,
        disulfide_max_sg_distance=qc_config.disulfide_max_sg_distance,
    )
    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / f"{construct_config.construct_id}_structural_qc.md").write_text(
        format_qc_report(qc_result),
        encoding="utf-8",
    )
    return qc_result


def _first_model_path(models_dir: str) -> Path:
    pdb_paths = sorted(
        Path(models_dir).glob("*.pdb"),
        key=lambda path: natural_sort_key(str(path)),
    )
    if not pdb_paths:
        raise ValueError(f"No PDB files found in models_dir: {models_dir}")
    return pdb_paths[0]


def _print_qc_summary(qc_result: StructuralQCResult) -> None:
    click.echo(f"Structural QC: {'PASS' if qc_result.fold_qc_pass else 'FOLD CORRUPTED'}")
    click.echo(f"His triad: {_status(qc_result.his_triad_pass)}")
    click.echo(f"Lys-Tyr: {_status(qc_result.lys_tyr_pass)}")
    disulfide_result = qc_result.details.get("disulfide", {})
    disulfide_status = "SKIPPED" if disulfide_result.get("skipped") else _status(qc_result.disulfide_pass)
    click.echo(f"Disulfide: {disulfide_status}")
    click.echo(f"Active-site access: {_status(qc_result.active_site_accessible)}")
    click.echo(f"pLDDT: {_status(qc_result.plddt_pass)}")
    click.echo(f"PAE: {_status(qc_result.pae_pass)}")


def _status(value: bool | None) -> str:
    if value is None:
        return "N/A"
    return "PASS" if value else "FAIL"


def _configure_stdio() -> None:
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            reconfigure(encoding="utf-8")


def _domain_tuple(domain: Any) -> tuple[int, int] | None:
    if domain is None:
        return None
    return domain.start, domain.end


def _build_summary(
    construct_config: ConstructConfig,
    pose_metrics: list[PoseMetrics],
    contact_results: list[dict[str, Any]],
) -> ConstructSummary:
    productive_count = sum(1 for pose in pose_metrics if pose.productive)
    total_poses = len(pose_metrics)
    best_distance = (
        min(pose.active_site_to_target_distance for pose in pose_metrics)
        if pose_metrics
        else float(np.inf)
    )
    best_productive_rank = next(
        (pose.rank for pose in pose_metrics if pose.productive),
        None,
    )
    frequencies = compute_contact_frequencies(contact_results)

    return ConstructSummary(
        construct_id=construct_config.construct_id,
        productive_count=productive_count,
        total_poses=total_poses,
        best_distance=best_distance,
        best_productive_rank=best_productive_rank,
        lox_contact_frequency=frequencies["lox_contact_frequency"],
        cbd_contact_frequency=frequencies["cbd_contact_frequency"],
        cbd_coupling=0.0,
        strict_pass=None,
    )


def _print_screen_summary(screen_result: Any, out: str) -> None:
    candidate = screen_result.candidate_summary
    divider = "══════════════════════════════════════"
    click.echo(divider)
    click.echo(f"Screen: {screen_result.screen_id}")
    click.echo(divider)
    click.echo(f"Candidate: {candidate.construct_id}")
    click.echo(f"Productive poses: {candidate.productive_count}/{candidate.total_poses}")
    click.echo(f"Best distance: {candidate.best_distance:.2f} Å")
    click.echo(f"Best productive rank: {candidate.best_productive_rank}")
    click.echo("")
    click.echo("Criteria:")
    for criterion in [
        "pass_count_vs_baseline",
        "best_distance_vs_baseline",
        "best_productive_rank",
        "cbd_contribution",
        "beats_scrambled",
        "beats_polyK",
    ]:
        click.echo(f"  {criterion + ':':<27} {_criteria_status(screen_result.criteria_results.get(criterion))}")
    click.echo("")
    click.echo(f"Decision: {'PASS' if screen_result.overall_pass else 'FAIL'}")
    click.echo(screen_result.decision_text)
    click.echo("")
    click.echo(f"Output saved to: {out}")
    click.echo(divider)


def _criteria_status(value: object) -> str:
    if value is True:
        return "PASS"
    if value is False:
        return "FAIL"
    return "SKIPPED"


cli.add_command(run)
cli.add_command(screen)
cli.add_command(check_config)
cli.add_command(inactive_control)


if __name__ == "__main__":
    cli()
