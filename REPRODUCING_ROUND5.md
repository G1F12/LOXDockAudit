# Reproducing the Round 5 HDOCK Validation

This document reproduces the real-data LOXDockAudit v0.2.1 validation example:
R5 LOX169-417 vs collagen-I 6VZX.

## Dataset

Packaged example:

```text
examples/real_hdock_round5/
├── README.md
├── config.yaml
├── cleaned_models/
└── expected_output/
```

The cleaned model directory contains HDOCK-ranked `model_1.pdb` through
`model_10.pdb`. LOXDockAudit uses natural numeric sorting, so rank 10 is
processed after rank 9.

## Commands

Run from the repository root:

```bash
pip install -e .
loxdockaudit check-config examples/real_hdock_round5/config.yaml
loxdockaudit run --config examples/real_hdock_round5/config.yaml --out examples/real_hdock_round5/expected_output --top-n 10 --verbose
```

If the console script is not on `PATH` on Windows, run it through the user
scripts directory or use your active virtual environment's `Scripts` directory.

## Expected CLI Summary

```text
Construct: r5_real_lox169
Productive poses: 2/10
Best distance: 6.16 A
Best productive rank: 3
Output saved to: examples/real_hdock_round5/expected_output
```

## Expected Output Files

```text
examples/real_hdock_round5/expected_output/r5_real_lox169_poses.csv
examples/real_hdock_round5/expected_output/r5_real_lox169_summary.csv
examples/real_hdock_round5/expected_output/r5_real_lox169_report.md
examples/real_hdock_round5/expected_output/r5_real_lox169_structural_qc.md
```

## Expected Metrics

| Metric | Expected value |
|--------|----------------|
| Productive poses | 2/10 |
| Best distance | 6.16 A |
| Best productive rank | 3 |
| Productive fraction | 0.20 |

Expected productive ranks are 3 and 5. `model_10.pdb` must appear as rank 10 in
the poses CSV.

## Structural QC Notes

Correct active-site residue numbers for this construct use full-length LOX
numbering on chain D:

- His triad: His292, His294, His296
- LTQ precursor residues: Lys320 and Tyr355
- Candidate disulfide pairs: Cys238-Cys291, Cys244-Cys361, Cys324-Cys340,
  Cys330-Cys351, Cys398-Cys412

QC result: PASS.

Observed structural QC distances in `model_1.pdb`:

- His292-His294 CA distance: 5.714 A
- His292-His296 CA distance: 9.158 A
- His294-His296 CA distance: 6.555 A
- Lys320-Tyr355 CB distance: 6.775 A
- Disulfides intact by the configured 2.5 A tolerance: 5/5

Interpretation: the real Round 5 LOX169-417 model retains a plausible active-site
arrangement before docking interpretation. The His triad is compact enough for
the current v0.3 check (`all pairs < 10 A`), Lys320 and Tyr355 are spatially
proximal, and the configured disulfides are present in the model.

## Interpretation

This validation confirms that the real HDOCK dataset is parsed successfully and
that productive-rank calculation follows HDOCK's numeric model ordering. It also
demonstrates why raw docking rank alone is insufficient: rank 1 is not
productive, while productive geometry first appears at rank 3.
