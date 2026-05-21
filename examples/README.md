# Examples

This directory contains lightweight inputs and reproducible workflow examples.

## `sample.pdb`

Minimal PDB file for smoke-testing the parser through:

```bash
python main.py --input examples/sample.pdb
```

## `real_hdock_round5/`

Cleaned HDOCK Round 5 LOX169-417 versus collagen-I 6VZX pose files used for productive-geometry regression checks.

## `real_af_multimer_example/`

Derived minimal AF-Multimer-style comparison example using the same LOX169-417 and collagen-I target context.
It demonstrates `loxdockaudit af-compare` without redistributing a large raw ColabFold output directory.
