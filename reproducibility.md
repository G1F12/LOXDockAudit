# Reproducibility

## Environment

- Python version: 3.10 or 3.11.
- Key dependencies with pinned versions:
  - biopython==1.87
  - numpy==2.2.6
  - pandas==2.3.1
  - pyyaml==6.0.2
  - click==8.1.8
  - pytest==9.0.2
  - matplotlib==3.10.5
  - seaborn==0.13.2

## Data

- HDOCK Round 5 PDB files are included in `examples/real_hdock_round5/cleaned_models/`.
- The derived AF-Multimer-style example input is included as `examples/real_af_multimer_example/ranked_0.pdb`.
- The minimal sample input for smoke testing is `examples/sample.pdb`.
- Source biological context: LOX169-417 and collagen-I 6VZX lysine-site target region.

## Steps to reproduce

1. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

2. Run the smoke-test parser:

   ```bash
   python main.py --input examples/sample.pdb
   ```

3. Run the AF comparison example:

   ```bash
   loxdockaudit af-compare \
     --config examples/real_af_multimer_example/r5_af_compare_real.yaml \
     --out examples/real_af_multimer_example/expected_output
   ```

4. Regenerate figures:

   ```bash
   python docs/figures/generate_figures.py
   ```

## Expected output

- `examples/real_af_multimer_example/expected_output/real_af_compare_af_compare.csv`
- `examples/real_af_multimer_example/expected_output/real_af_compare_interface_overlap.csv`
- `examples/real_af_multimer_example/expected_output/real_af_compare_contact_persistence.csv`
- `examples/real_af_multimer_example/expected_output/real_af_compare_report.md`
- `docs/figures/figure1_pipeline.png`
- `docs/figures/figure2_active_inactive.png`
- `docs/figures/figure3_af_hdock_overlap.png`

## Validation

Verify that the AF comparison CSV reports 2/10 HDOCK productive poses, AF productive status `True`, interface overlap count 26, Jaccard index 0.667, and convergence label `strong convergence`.
Then run:

```bash
python -m pytest -q
python -m ruff check .
python -m mypy src/loxdockaudit
```
