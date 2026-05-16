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

## Interpretation

This validation confirms that the real HDOCK dataset is parsed successfully and
that productive-rank calculation follows HDOCK's numeric model ordering. It also
demonstrates why raw docking rank alone is insufficient: rank 1 is not
productive, while productive geometry first appears at rank 3.
