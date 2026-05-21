# LOXDockAudit

LOXDockAudit audits LOX-collagen docking and AF-Multimer comparison outputs with reproducible productive-geometry metrics for cautious structural-bioinformatics interpretation.

![CI](https://github.com/G1F12/LOXDockAudit/actions/workflows/ci.yml/badge.svg)
![Python](https://img.shields.io/badge/python-3.10%20|%203.11-blue)
![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)
![mypy](https://img.shields.io/badge/mypy-checked-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Release](https://img.shields.io/badge/release-v0.5.1-orange)

## Overview

Docking rank alone does not show whether LOX active-site residues are positioned near a target collagen Lys/Hyl side chain.
LOXDockAudit measures active-site distance, optional side-chain orientation, interface overlap, control behavior, and AF-Multimer agreement using deterministic Python workflows.
It writes CSV and Markdown outputs that make productive-geometry assumptions explicit and reproducible.
The project is computational only and does not claim enzymatic activity or biological efficacy.

## Quick Start

```bash
pip install -r requirements.txt
python main.py --input examples/sample.pdb
```

## Repository Structure

```text
project/
├── README.md
├── docs/
│   └── visual_summary.md
├── technical_report.md
├── examples/
│   └── sample.pdb
└── reproducibility.md
```

## Results Summary

Round 5 LOX169-417 analysis finds 2/10 productive HDOCK poses, while the inactive H292A/H294A/H296A control finds 1/10.
The derived AF-Multimer-style example reports strong computational interface convergence with HDOCK under the configured contact metric.
See [docs/visual_summary.md](docs/visual_summary.md) for the figures and one-page interpretation.

## Citation

See `CITATION.cff`.

## License

MIT
