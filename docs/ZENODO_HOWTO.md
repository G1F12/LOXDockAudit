# How to publish LOXDockAudit to Zenodo

## Prerequisites

- GitHub repository with a clean v0.4.0 release tag
- Zenodo account connected to GitHub
- Updated `CITATION.cff` and `.zenodo.json`

## Steps

### 1. Enable the repository

Go to <https://zenodo.org/account/settings/github/> and enable the `LOXDockAudit` repository.

### 2. Create the GitHub release

Go to <https://github.com/G1F12/LOXDockAudit/releases> and draft a release.

- Tag: `v0.4.0`
- Title: `LOXDockAudit v0.4.0 - Orientation and inactive-control release`

Use release notes that describe the software changes without implying experimental validation. The release should state that LOXDockAudit is a computational docking-analysis and reproducibility framework.

### 3. Confirm Zenodo metadata

After publishing the GitHub release, wait a few minutes for Zenodo to archive it. Confirm that the record shows:

- version `v0.4.0`
- MIT license
- GitHub repository link
- software title from `CITATION.cff`
- no claims of enzymatic activity, therapeutic efficacy, or biomechanical validation

### 4. Update DOI badge if needed

If Zenodo creates a new version-specific DOI, update the README badge only after confirming the DOI record.

## Release caution

LOXDockAudit reports computational geometry metrics from docking poses. It does not validate LOX activity, collagen oxidation, crosslink formation, tissue strengthening, or in vivo behavior.
