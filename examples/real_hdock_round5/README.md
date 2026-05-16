# Real HDOCK Round 5 Example

This example packages the R5 LOX169-417 vs collagen-I 6VZX validation dataset
used for the LOXDockAudit v0.2.1 release checks.

Run from the repository root:

```bash
loxdockaudit check-config examples/real_hdock_round5/config.yaml
loxdockaudit run --config examples/real_hdock_round5/config.yaml --out examples/real_hdock_round5/expected_output --top-n 10
```

Expected metrics:

- Productive poses: 2/10
- Best distance: 6.16 A
- Best productive rank: 3

The cleaned models are ordered by numeric model rank, so `model_10.pdb` is
processed after `model_9.pdb`, not between `model_1.pdb` and `model_2.pdb`.
