# LOXDockAudit v0.2.1 Release Notes

LOXDockAudit v0.2.1 is a reproducibility and parser-hardening release for the
Round 5 real HDOCK validation dataset.

## Highlights

- Natural numeric sorting for model files: `model_1.pdb`, `model_2.pdb`, ...,
  `model_10.pdb`.
- Real HDOCK compatibility validated on R5 LOX169-417 vs collagen-I 6VZX.
- Official PDB preprocessing for HDOCK wrapper and malformed records.
- Strict productive-rank correction for v0.2 criteria.
- Full typing and lint checks remain clean.
- 105-test baseline passed before this release pass; preprocessing tests added.

## Real HDOCK Result

Dataset: R5 LOX169-417 vs collagen-I 6VZX

- Productive poses: 2/10
- Best distance: 6.16 A
- Best productive rank: 3

## Reproduce

```bash
pip install -e .
loxdockaudit check-config examples/real_hdock_round5/config.yaml
loxdockaudit run --config examples/real_hdock_round5/config.yaml --out examples/real_hdock_round5/expected_output --top-n 10
```

Expected output files:

- `examples/real_hdock_round5/expected_output/r5_real_lox169_poses.csv`
- `examples/real_hdock_round5/expected_output/r5_real_lox169_summary.csv`
- `examples/real_hdock_round5/expected_output/r5_real_lox169_report.md`

## Git Tag

Recommended tag command from a real git checkout:

```bash
git tag -a v0.2.1 -m "LOXDockAudit v0.2.1"
git push origin v0.2.1
```
