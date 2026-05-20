"""HDOCK versus AF-Multimer interface convergence analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from loxdockaudit.interface_extractor import InterfaceResidues


@dataclass(frozen=True)
class ConvergenceReport:
    """Summary of residue-level overlap between HDOCK and AF-Multimer interfaces."""

    converged_residues: frozenset[int]
    hdock_only_residues: frozenset[int]
    af_only_residues: frozenset[int]
    jaccard_index: float
    overlap_count: int
    convergence_narrative: str


@dataclass(frozen=True)
class ResidueClassification:
    """Partition of interface residues by method support."""

    converged: frozenset[int]
    hdock_only: frozenset[int]
    af_only: frozenset[int]
    total_unique: int


@dataclass(frozen=True)
class CorrelationReport:
    """High-level correlation between productive HDOCK poses and AF models."""

    productive_in_both: int
    productive_hdock_only: int
    productive_af_only: int
    summary: str


def compute_interface_overlap(
    hdock_interface: InterfaceResidues,
    af_interface: InterfaceResidues,
) -> ConvergenceReport:
    """Compute residue overlap and narrative for HDOCK and AF-Multimer interfaces."""
    classification = classify_residues(hdock_interface, af_interface)
    union_count = classification.total_unique
    overlap_count = len(classification.converged)
    jaccard_index = overlap_count / union_count if union_count else 0.0
    return ConvergenceReport(
        converged_residues=classification.converged,
        hdock_only_residues=classification.hdock_only,
        af_only_residues=classification.af_only,
        jaccard_index=jaccard_index,
        overlap_count=overlap_count,
        convergence_narrative=_convergence_narrative(jaccard_index),
    )


def classify_residues(
    hdock: InterfaceResidues,
    af: InterfaceResidues,
) -> ResidueClassification:
    """Classify residues as converged, HDOCK-only, or AF-only."""
    hdock_residues = _combined_residues(hdock)
    af_residues = _combined_residues(af)
    converged = hdock_residues & af_residues
    hdock_only = hdock_residues - af_residues
    af_only = af_residues - hdock_residues
    return ResidueClassification(
        converged=frozenset(sorted(converged)),
        hdock_only=frozenset(sorted(hdock_only)),
        af_only=frozenset(sorted(af_only)),
        total_unique=len(hdock_residues | af_residues),
    )


def productive_geometry_correlation(
    hdock_poses: list[Any],
    af_models: list[Any],
) -> CorrelationReport:
    """Compare productive geometry calls between HDOCK poses and AF models."""
    hdock_productive = any(_productive(item) for item in hdock_poses)
    af_productive = any(_productive(item) for item in af_models)
    return CorrelationReport(
        productive_in_both=int(hdock_productive and af_productive),
        productive_hdock_only=int(hdock_productive and not af_productive),
        productive_af_only=int(af_productive and not hdock_productive),
        summary=_correlation_summary(hdock_productive, af_productive),
    )


def _combined_residues(interface: InterfaceResidues) -> frozenset[int]:
    return frozenset(interface.chain_a_residues | interface.chain_b_residues)


def _convergence_narrative(jaccard_index: float) -> str:
    if jaccard_index >= 0.5:
        return (
            "strong convergence: the methods share a substantial interface "
            "residue set under the current contact metric"
        )
    if jaccard_index > 0.0:
        return (
            "partial convergence: the methods share some interface residues "
            "under the current contact metric"
        )
    return (
        "divergent interfaces: no shared interface residues were detected under "
        "the current contact metric"
    )


def _productive(item: Any) -> bool:
    for attribute in ("fully_productive", "is_fully_productive", "productive", "is_productive"):
        if hasattr(item, attribute):
            return bool(getattr(item, attribute))
    if isinstance(item, dict):
        return bool(
            item.get("fully_productive")
            or item.get("is_fully_productive")
            or item.get("productive")
            or item.get("is_productive")
        )
    return False


def _correlation_summary(hdock_productive: bool, af_productive: bool) -> str:
    if hdock_productive and af_productive:
        return "productive geometry was observed in both workflows under configured heuristics"
    if hdock_productive:
        return "productive geometry was observed in HDOCK only under configured heuristics"
    if af_productive:
        return "productive geometry was observed in AF-Multimer only under configured heuristics"
    return "productive geometry was not observed in either workflow under configured heuristics"
