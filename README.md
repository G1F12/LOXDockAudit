# LOXDockAudit

[![CI](https://github.com/G1F12/LOXDockAudit/actions/workflows/ci.yml/badge.svg)](https://github.com/G1F12/LOXDockAudit/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20241797.svg)](https://doi.org/10.5281/zenodo.20241797)

LOXDockAudit is a reproducible Python framework for auditing LOX-collagen docking poses by productive geometry rather than docking score alone. It is intended for computational docking analysis, reproducibility checks, and control-aware interpretation of lysyl oxidase (LOX) docking screens.

This project is computational only. It does not demonstrate enzymatic activity, collagen oxidation, crosslink formation, biomechanical improvement, therapeutic utility, or in vivo safety.

## Why This Exists

Generic docking engines rank poses by binding-oriented objectives. For LOX-collagen systems, a high-ranking pose is not necessarily catalytically interpretable unless the LOX active-site region is positioned near the target collagen Lys/Hyl residue and the substrate side chain has a plausible approach orientation.

LOXDockAudit adds post-docking checks that make these assumptions explicit:

- active-site-to-substrate distance scoring
- optional Lys/Hyl side-chain orientation scoring
- fully productive pose scoring when distance and orientation both pass
- control-aware comparisons across baseline, scrambled, polyK, noCBD, inactive, and benchmark constructs
- structural QC for active-site geometry, disulfides, and accessibility proxies
- deterministic natural sorting of docking models
- SHA256 input auditing and duplicate detection
- CSV, Markdown, and SVG reporting

## Installation

LOXDockAudit requires Python 3.10 or newer.

```bash
git clone https://github.com/G1F12/LOXDockAudit.git
cd LOXDockAudit
python -m pip install -e ".[dev]"
```

Run checks:

```bash
python -m pytest -q
python -m ruff check .
python -m mypy src
```

Current v0.4.0 QA status: `215 passed`, `ruff` passed, `mypy` passed.

## Quickstart

Validate a construct config:

```bash
loxdockaudit check-config configs/example.yaml
```

Run a single construct:

```bash
loxdockaudit run \
  --config configs/example.yaml \
  --out examples/reports \
  --top-n 10
```

Run the Round 5 inactive-control comparison:

```bash
loxdockaudit inactive-control \
  --config configs/r5_inactive_control.yaml \
  --out examples/real_hdock_inactive_control/expected_output \
  --top-n 10
```

## Configuration

Constructs are defined in YAML. The core fields identify docking models, receptor and ligand chains, target substrate residue, and active-site residues.

```yaml
construct_id: r5_real_lox169
construct_type: lox_baseline
models_dir: examples/real_hdock_round5/cleaned_models
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

active_site:
  residues:
    - {label: "His292", resi: 292}
    - {label: "His294", resi: 294}
    - {label: "His296", resi: 296}
```

## Orientation Scoring

v0.4.0 adds optional substrate orientation scoring:

```yaml
orientation:
  enabled: true
  threshold_deg: 90.0
  sidechain_atoms: ["CB", "NZ"]
  activesite_centroid_atoms: ["CA"]
```

The orientation score is the angle between:

- the vector from the active-site centroid to the target Lys/Hyl NZ atom
- the vector from target Lys/Hyl C-beta to NZ

Small angles indicate that the side chain points toward the active-site region. When orientation is enabled, a fully productive pose must satisfy both the distance threshold and the orientation threshold. When orientation is disabled, LOXDockAudit preserves the v0.3 distance-only behavior.

## Real Round 5 Reproduction

The repository includes a reproducible HDOCK Round 5 example:

- construct: LOX169-417
- receptor: collagen-I 6VZX lysine-site model
- poses analyzed: 10
- docking engine: HDOCK

Distance-only v0.3-compatible results:

| Construct | Productive poses | Best distance | Best productive rank |
|---|---:|---:|---:|
| Active LOX169-417 | 2/10 | 6.159 A | 3 |
| Inactive H292A/H294A/H296A | 1/10 | 6.421 A | 7 |

With v0.4 orientation scoring enabled, the active construct has one fully productive pose and the inactive catalytic-control construct has zero fully productive poses under the current geometric heuristic.

Interpretation: the active construct retained more favorable productive geometry under this computational metric. This is not evidence of enzymatic activity or biological efficacy.

## Outputs

LOXDockAudit writes:

- per-pose CSV files with distance, productivity, optional orientation, and contacts
- summary CSV files with productive and fully productive counts
- screen and inactive-control comparison CSV files
- Markdown reports and supplements
- SVG distance histograms and geometry scatter plots generated with pure Python

## Reproducibility Focus

LOXDockAudit treats docking poses as auditable evidence:

- numeric natural sorting prevents `model_10.pdb` from being ranked before `model_2.pdb`
- SHA256 input auditing detects duplicate or reused structures
- expected outputs are included for the real Round 5 examples
- regression tests lock down v0.3 distance metrics while adding v0.4 orientation outputs

## Limitations

LOXDockAudit is a post-docking analysis framework. It does not prove catalysis.

Important limitations:

- Distance and orientation thresholds are geometric heuristics.
- Docking poses are static and do not model conformational dynamics.
- Structural QC uses practical geometry proxies, not full physical validation.
- HDOCK and similar docking engines do not model copper loading, LTQ/topaquinone chemistry, or catalytic turnover.
- Productive geometry is not equivalent to enzymatic activity, collagen oxidation, crosslink formation, tissue strengthening, medical efficacy, or in vivo safety.

## Documentation

- [Technical report](docs/technical_report.md)
- [Usage guide](docs/USAGE.md)
- [Round 5 reproduction notes](REPRODUCING_ROUND5.md)
- [Inactive-control workflow](docs/INACTIVE_CONTROL.md)
- [PDB preprocessing notes](docs/PDB_PREPROCESSING.md)

## Citation

If you use LOXDockAudit in research or teaching, cite the software metadata in [CITATION.cff](CITATION.cff).

```bibtex
@software{loxdockaudit_2026,
  title = {LOXDockAudit: Control-aware productive-geometry analysis for LOX-collagen docking screens},
  author = {Karatseyeu, Aliaksandr},
  year = {2026},
  version = {0.4.0},
  url = {https://github.com/G1F12/LOXDockAudit}
}
```

## Roadmap

Near-term priorities:

- add compact PyMOL/ChimeraX visualization exports for selected poses
- improve documentation around orientation-score interpretation
- add benchmark examples beyond Round 5
- separate lightweight example outputs from large local validation artifacts
- add optional JSON outputs for downstream analysis

## License

MIT License. See [LICENSE](LICENSE).

## Disclaimer

LOXDockAudit is research software for computational docking analysis. It is not a diagnostic, therapeutic, clinical, or wet-lab validation tool.
