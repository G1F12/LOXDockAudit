# LOXDockAudit Project Summary

## Problem

Docking score alone is not enough for LOX-collagen engineering. A top-ranked
docked pose can still be biologically unhelpful if the LOX active-site region is
not positioned near a collagen Lys/Hyl target residue. LOXDockAudit addresses
this gap by auditing whether docked poses satisfy a defined productive-geometry
criterion.

## Tool

LOXDockAudit is a Python command-line tool for productive geometry analysis in
LOX-collagen docking screens. It parses docking poses, preserves natural numeric
model ordering, validates inputs, measures active-site-to-target distances,
counts contacts, and reports productive pose summaries for candidate and control
constructs.

## Scientific Motivation

For enzyme-matrix systems, generic binding is not the same as plausible catalytic
positioning. LOXDockAudit is designed to separate nonspecific docking success
from geometry that places a LOX active-site region near a collagen substrate
residue. Negative-control-aware reporting helps compare candidate constructs
against baseline, scrambled, polyK, inactive, or benchmark controls.

## Real Round 5 Reproduction

The repository includes a reproducible HDOCK Round 5 example for LOX169-417
against collagen-I 6VZX. The expected result is:

- Productive poses: 2 / 10
- Best distance: 6.159 A
- Best productive rank: 3
- Productive distance threshold: 8.000 A

The README includes a real output excerpt and a plot of active-site-to-target
distance across the ten HDOCK-ranked models.

## Validation Status

The current release has automated tests and continuous integration on GitHub
Actions for Python 3.11 and 3.12. Local and CI validation cover unit tests,
integration tests, parser behavior, reporting, `ruff`, and `mypy`.

## Limitations

LOXDockAudit does not prove enzymatic activity, collagen oxidation, crosslink
formation, tendon strengthening, or in vivo safety. It is a reproducibility and
triage tool for docking outputs. Results depend on input structure quality,
docking protocol, residue mapping, and the declared productive-geometry
threshold.

## Future Work

Future work includes broader benchmark datasets, richer negative-control
templates, optional structural visualizations, more detailed contact context,
and export formats for manuscript supplements or screening reports.
