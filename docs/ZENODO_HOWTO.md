# How to publish LOXDockAudit to Zenodo

## Prerequisites
- GitHub account (G1F12)
- Zenodo account (create at https://zenodo.org using GitHub login)

## Steps

### Step 1 - Create Zenodo account
Go to https://zenodo.org
Click "Log in with GitHub"
Authorize Zenodo to access your GitHub account

### Step 2 - Enable the repository
Go to https://zenodo.org/account/settings/github/
Find "LOXDockAudit" in the repository list
Toggle it ON (flip the switch to enabled)

### Step 3 - Create GitHub Release v0.3.0
Go to https://github.com/G1F12/LOXDockAudit/releases
Click "Draft a new release"
Tag: v0.3.0
Title: LOXDockAudit v0.3.0 - Structural QC and Round 5 Reproducibility Release
Description (copy this exactly):

---
## What's new in v0.3.0

### Structural QC pipeline
- His triad CA geometry check (His292/His294/His296)
- Lys320-Tyr355 CB distance check
- Disulfide bond geometry check (5/5 intact in Round 5 dataset)
- Active-site burial proxy (SASA approximation)
- pLDDT/PAE parsing from AlphaFold2/ColabFold outputs
- Fold-corrupted flag with --strict-qc CLI option

### Real data validation
- LOX169-417 vs collagen-I 6VZX (HDOCK Round 5)
- Structural QC: PASS (all 4 checks)
- Productive poses: 2/10
- Best distance: 6.16 A
- Best productive rank: 3

### Bug fix: lexicographic file ordering
- model_10.pdb was sorted before model_2.pdb
- This changed best_productive_rank from 3 to 4
- Fixed with natural sort key
- Regression test added

### Test coverage
- 182 tests passing
- mypy: 0 errors
- ruff: 0 violations

### Documentation
- Technical report: docs/technical_report.md
- Reproducibility guide: REPRODUCING_ROUND5.md
---

Click "Publish release"

### Step 4 - Zenodo picks it up automatically
Wait 2-5 minutes after publishing the GitHub release.
Go to https://zenodo.org/account/settings/github/
You should see LOXDockAudit with a DOI assigned.
Click on the record to see your DOI (format: 10.5281/zenodo.XXXXXXX)

### Step 5 - Update README badge
Replace the placeholder DOI badge in README.md with the real one:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Replace XXXXXXX with your actual Zenodo record number.

Commit:
  git add README.md
  git commit -m "docs: add Zenodo DOI badge"
  git push

### Step 6 - Verify
Go to https://doi.org/10.5281/zenodo.XXXXXXX
Confirm the record shows:
- Correct title
- Correct version v0.3.0
- Correct author name
- MIT license
- Link to GitHub repository

## Result
You now have a citable DOI for LOXDockAudit v0.3.0.
This can be included in:
- Application materials
- Portfolio
- Any future publication or poster
