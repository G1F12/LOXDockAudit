# LOXDockAudit Usage Guide

LOXDockAudit is a command-line tool for auditing LOX-collagen docking screens. It helps separate generic binding from potentially productive catalytic geometry by checking whether the LOX active site is positioned near a target collagen Lys/Hyl residue.

## Install

From the repository root:

```bash
pip install -e .
```

Check that the CLI is available:

```bash
loxdockaudit --help
```

Available commands:

```text
loxdockaudit run
loxdockaudit screen
loxdockaudit check-config
loxdockaudit inactive-control
loxdockaudit af-compare
```

## Single Construct Analysis

Use `run` when you want to analyze one construct.

```bash
loxdockaudit check-config configs/example.yaml
loxdockaudit run --config configs/example.yaml --out examples/reports --top-n 3
```

Options:

```text
--config   Path to construct YAML config
--out      Output directory
--top-n    Number of top-ranked PDB poses to analyze
--verbose  Print additional audit output
```

Single-construct outputs:

```text
{construct_id}_poses.csv
{construct_id}_summary.csv
{construct_id}_report.md
```

## Screen-Level Analysis

Use `screen` when you want to compare a candidate construct against baseline and control constructs.

```bash
loxdockaudit screen --config configs/example_screen.yaml --out examples/reports --top-n 3
```

Screen outputs:

```text
{screen_id}_control_comparison.csv
{screen_id}_summaries.csv
{screen_id}_report.md
{construct_id}_poses.csv
```

## Inactive-Control Analysis

Use `inactive-control` to compare active LOX against a catalytically inactive
control with the same docking and productive-geometry rules.

```bash
loxdockaudit inactive-control --config configs/r5_inactive_control.yaml --out examples/real_hdock_inactive_control/expected_output --top-n 10
```

Inactive-control outputs:

```text
{comparison_id}_comparison.csv
{comparison_id}_distance_distribution.csv
{comparison_id}_distance_histogram.svg
{comparison_id}_contact_persistence.csv
{comparison_id}_pose_clusters.csv
{comparison_id}_supplement.md
```

## Independent Structural Prediction Comparison

Use `af-compare` to compare productive HDOCK poses against an
AlphaFold-Multimer or ColabFold predicted complex for the same LOX construct
and collagen target region.

```bash
loxdockaudit af-compare \
  --config examples/real_af_multimer_example/r5_af_compare_real.yaml \
  --out examples/real_af_multimer_example/expected_output
```

AF comparison outputs:

```text
{comparison_id}_af_compare.csv
{comparison_id}_interface_overlap.csv
{comparison_id}_contact_persistence.csv
{comparison_id}_report.md
```

Minimal AF comparison config:

```yaml
comparison_id: example_hdock_vs_af
hdock_models_dir: examples/real_hdock_round5/cleaned_models
alphafold_model_path: examples/af_multimer/ranked_0.pdb

target_residue:
  chain: "A"
  resi: 11
  resn: "LYS"
  atom: "NZ"

active_site_chain: "D"
active_site_residues: [292, 294, 296, 320, 355]

productive_distance_threshold: 8.0
contact_distance_threshold: 5.0
orientation_enabled: true
orientation_threshold_deg: 90.0

compare_contact_fingerprints: true
compare_interface_overlap: true
compare_productive_geometry: true
```

The convergence categories are heuristic:

- `strong convergence`: substantial residue overlap under the configured contact metric.
- `partial convergence`: some shared interface residues.
- `divergent interfaces`: no shared interface residues detected.

These labels describe computational agreement only. They do not establish
enzymatic activity, catalytic mechanism, collagen oxidation, crosslink
formation, or biological efficacy.

## Construct Config Format

Minimal construct config:

```yaml
construct_id: test_lox_baseline
construct_type: lox_baseline
models_dir: examples/data/lox_baseline
productive_distance_threshold: 8.0

receptor:
  chains: ["A", "B"]
  target_residues:
    - chain: "A"
      resi: 11
      resn: "LYS"
      atom: "NZ"

ligand:
  chains: ["D"]

lox_domain:
  start: 1
  end: 200

active_site:
  residues:
    - {label: "His292", resi: 124}
    - {label: "His294", resi: 126}
    - {label: "His296", resi: 128}
    - {label: "Lys320", resi: 152}
    - {label: "Tyr355", resi: 187}
```

Important fields:

- `models_dir`: folder containing `*.pdb` docking poses.
- `receptor.chains`: collagen/receptor chains.
- `ligand.chains`: LOX or fusion construct chains.
- `receptor.target_residues`: substrate residue or residues to check.
- `active_site.residues`: LOX active-site residue numbers in the model numbering.
- `productive_distance_threshold`: distance cutoff in angstroms. Default is `8.0`.

## Screen Config Format

Minimal screen config:

```yaml
screen_id: example_screen
candidate_construct_id: example_candidate

constructs:
  - construct_id: example_baseline
    control_type: baseline
    config_path: example.yaml

  - construct_id: example_candidate
    control_type: candidate
    config_path: example.yaml

  - construct_id: example_scrambled
    control_type: scrambled
    config_path: example.yaml

  - construct_id: example_polyK
    control_type: polyK
    config_path: example.yaml

strict_criteria:
  pass_count_greater_than_baseline: true
  best_distance_lte_baseline: true
  best_productive_rank_max: 3
  cbd_must_contribute: false
  must_beat_scrambled: true
  must_beat_polyK: true
```

`config_path` values are resolved relative to the screen YAML file if they are not absolute paths.

Supported `control_type` values:

```text
candidate
baseline
scrambled
polyK
noCBD
inactive
benchmark
```

## What The Metrics Mean

`productive`:
The pose is considered productive if the minimum heavy-atom distance between the configured active-site residues and the target substrate residue is less than or equal to the configured threshold.

`best_productive_rank`:
The first ranked pose that is productive. `None` means no productive pose was found.

`lox_contacts`:
Number of ligand LOX-domain atom pairs contacting receptor atoms within the contact cutoff.

`cbd_contacts`:
Number of ligand CBD/FMOD-domain atom pairs contacting receptor atoms. This is `0` if no CBD domain is configured.

`cbd_coupling`:
Pearson correlation between CBD contacts and active-site distance across poses. This is a triage metric, not a mechanistic proof.

`interface_jaccard_index`:
Residue-overlap score between AF-Multimer and productive HDOCK interfaces. It is
a reproducibility and cross-method comparison metric, not validation of
catalysis.

## Input Audit

LOXDockAudit audits PDB inputs before analysis:

- Checks file existence.
- Checks PDB parseability.
- Computes SHA-256 hashes.
- Checks expected chains.
- Detects duplicate SHA-256 files across constructs.

Cross-construct duplicate files are treated as critical audit errors.

## Typical Workflow

1. Prepare one folder of PDB files per construct.
2. Create one construct YAML per construct.
3. Run `check-config` on each construct config.
4. Create a screen YAML containing candidate, baseline, and controls.
5. Run `loxdockaudit screen`.
6. Review the Markdown report and CSV outputs.

Example:

```bash
loxdockaudit check-config configs/example.yaml
loxdockaudit screen --config configs/example_screen.yaml --out examples/reports --top-n 3 --verbose
```

## Limitations

LOXDockAudit does not prove enzymatic activity, Cu loading, LTQ/topaquinone maturation, collagen oxidation, crosslink formation, tissue strengthening, or in vivo safety.

It is a computational triage and reproducibility tool for docking-screen interpretation.

## Troubleshooting

If `loxdockaudit` is not found after installation, use:

```bash
python -m pip install -e .
```

On Windows, the script may be installed under:

```text
%APPDATA%\Python\Python311\Scripts
```

If a run fails with an input audit error, check:

- `models_dir` exists.
- PDB files are present.
- Chain IDs match the YAML config.
- The PDB files are not empty.
- Candidate and controls are not accidental duplicate files unless intended.
