from __future__ import annotations

from loxdockaudit.af_compare import convergence_category
from loxdockaudit.convergence import (
    classify_residues,
    compute_interface_overlap,
    productive_geometry_correlation,
)
from loxdockaudit.interface_extractor import InterfaceResidues


def test_interface_overlap_logic_partial() -> None:
    hdock = InterfaceResidues(frozenset({10, 11}), frozenset({100}))
    af = InterfaceResidues(frozenset({11, 12}), frozenset({101}))

    report = compute_interface_overlap(hdock, af)

    assert report.converged_residues == frozenset({11})
    assert report.hdock_only_residues == frozenset({10, 100})
    assert report.af_only_residues == frozenset({12, 101})
    assert report.jaccard_index == 0.2
    assert convergence_category(report) == "partial convergence"


def test_convergence_categorization_strong_and_divergent() -> None:
    strong = compute_interface_overlap(
        InterfaceResidues(frozenset({1, 2}), frozenset({3})),
        InterfaceResidues(frozenset({1, 2}), frozenset({4})),
    )
    divergent = compute_interface_overlap(
        InterfaceResidues(frozenset({1}), frozenset({2})),
        InterfaceResidues(frozenset({3}), frozenset({4})),
    )

    assert convergence_category(strong) == "strong convergence"
    assert convergence_category(divergent) == "divergent interfaces"


def test_classify_residues_total_unique() -> None:
    result = classify_residues(
        InterfaceResidues(frozenset({1}), frozenset({2})),
        InterfaceResidues(frozenset({2}), frozenset({3})),
    )

    assert result.converged == frozenset({2})
    assert result.total_unique == 3


def test_productive_geometry_correlation_reads_dataclass_like_attrs() -> None:
    class Pose:
        productive = True

    class Model:
        is_productive = False

    report = productive_geometry_correlation([Pose()], [Model()])

    assert report.productive_hdock_only == 1
    assert "HDOCK only" in report.summary
