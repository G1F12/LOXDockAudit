# Real HDOCK Inactive-Control Example

This example compares active LOX169-417 against the catalytically inactive
H292A/H294A/H296A control using paired Round 5 HDOCK runs against the same
6VZX lysine-site receptor.

Run from the repository root:

```bash
loxdockaudit inactive-control --config configs/r5_inactive_control.yaml --out examples/real_hdock_inactive_control/expected_output --top-n 10
```

The inactive model was generated as a ColabFold/AlphaFold2-ptm single-chain
prediction, then docked with HDOCKlite using the same receptor and top-10 model
pipeline as the active control.

This is a computational structural-control experiment. It does not test or prove
enzymatic activity.
