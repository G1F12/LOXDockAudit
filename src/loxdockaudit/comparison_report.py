"""Final HDOCK versus AF-Multimer comparison reporting."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from loxdockaudit.convergence import ConvergenceReport


@dataclass(frozen=True)
class ComparisonReport:
    """Combined HDOCK, AF-Multimer, and convergence interpretation."""

    hdock_summary: Any
    af_summary: Any
    convergence: ConvergenceReport
    scientific_narrative: str
    version: str = "0.5.1"


def generate_comparison_report(
    hdock_results: Any,
    af_results: Any,
    convergence: ConvergenceReport,
) -> ComparisonReport:
    """Build a structured comparison report from HDOCK and AF-Multimer results."""
    return ComparisonReport(
        hdock_summary=hdock_results,
        af_summary=af_results,
        convergence=convergence,
        scientific_narrative=(
            "The comparison reports computational agreement between docking and "
            "AF-Multimer under the current geometric metric. It does not "
            "establish enzymatic activity."
        ),
    )


def render_comparison_markdown(report: ComparisonReport) -> str:
    """Render a comparison report as Markdown."""
    lines = [
        "# HDOCK vs AF-Multimer Comparison",
        "",
        f"Version: {report.version}",
        "",
        "## Summary",
        "",
        report.scientific_narrative,
        "",
        "## Interface Convergence",
        "",
        f"- Shared residues: {report.convergence.overlap_count}",
        f"- Jaccard index: {report.convergence.jaccard_index:.3f}",
        f"- Interpretation: {report.convergence.convergence_narrative}",
        "",
        "This computational agreement does not establish enzymatic activity, "
        "copper loading, LTQ chemistry, collagen oxidation, or biological efficacy.",
    ]
    return "\n".join(lines).rstrip() + "\n"


def render_comparison_svg(report: ComparisonReport) -> str:
    """Render a comparison report summary as SVG."""
    width = 640
    height = 220
    shared = report.convergence.overlap_count
    jaccard = report.convergence.jaccard_index
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        '<rect width="100%" height="100%" fill="white"/>\n'
        '<text x="32" y="40" font-family="Arial" font-size="18" '
        'font-weight="700">HDOCK vs AF-Multimer interface comparison</text>\n'
        f'<text x="32" y="82" font-family="Arial" font-size="14">Shared residues: {shared}</text>\n'
        f'<text x="32" y="112" font-family="Arial" font-size="14">Jaccard index: {jaccard:.3f}</text>\n'
        '<text x="32" y="150" font-family="Arial" font-size="12">'
        'Computational agreement only; not enzymatic validation.</text>\n'
        "</svg>\n"
    )
