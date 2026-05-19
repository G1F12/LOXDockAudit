# LOX169-417 Inactive-Control Validation

This workflow compares active LOX169-417 with the catalytically inactive
H292A/H294A/H296A control. It is a computational structural-control experiment,
not an enzymatic-activity assay.

Run from the repository root:

```bash
python -m loxdockaudit.cli inactive-control --config configs/r5_inactive_control.yaml --out examples/real_hdock_inactive_control/expected_output --top-n 10
```

Expected current result:

| Metric | Active | Inactive |
| --- | --- | --- |
| Productive poses | 2/10 | 1/10 |
| Best distance | 6.159 A | 6.421 A |
| Best productive rank | 3 | 7 |
| Productive fraction | 0.200 | 0.100 |
| QC status | PASS | PASS |
| Active-site geometry | intact | disrupted |

Interpretation:

Active LOX169-417 retained more favorable productive geometry under the current
docking metric. This is not proof of enzymatic activity.

Primary outputs:

```text
examples/real_hdock_inactive_control/expected_output/r5_lox169_active_vs_inactive_supplement.md
examples/real_hdock_inactive_control/expected_output/r5_lox169_active_vs_inactive_comparison.csv
examples/real_hdock_inactive_control/expected_output/r5_lox169_active_vs_inactive_distance_histogram.svg
examples/real_hdock_inactive_control/expected_output/r5_lox169_active_vs_inactive_distance_distribution.csv
examples/real_hdock_inactive_control/expected_output/r5_lox169_active_vs_inactive_contact_persistence.csv
examples/real_hdock_inactive_control/expected_output/r5_lox169_active_vs_inactive_pose_clusters.csv
```

The inactive ColabFold provenance and paired HDOCKlite output are stored under
`examples/real_hdock_inactive_control/provenance`.
