# Changelog

## v0.2.1 - Round 5 Reproducibility Release

- Fixed numeric/natural sorting for HDOCK model filenames, so `model_10.pdb` no longer ranks before `model_2.pdb`.
- Added official HDOCK-oriented PDB preprocessing for `MODEL`, `ENDMDL`, `REMARK`, `TER`, malformed coordinate rows, and blank chain IDs.
- Added real HDOCK Round 5 example package with cleaned models, config, and expected outputs.
- Added Round 5 reproducibility documentation with exact commands and expected metrics.
- Confirmed real HDOCK compatibility on R5 LOX169-417 vs collagen-I 6VZX.
- Corrected strict productive-rank handling in v0.2 criteria.
- Kept full typing and lint checks clean.
- Verified 105-test baseline before release hardening; added preprocessing coverage.

## v0.1.0 - Validated MVP

- Added CLI commands: run, check-config
- Added productive pose analysis
- Added distance/contact reporting
- Added CSV and Markdown outputs
- Added duplicate SHA audit detection
- Added Windows UTF-8 CLI handling
- Added pytest suite
- Added integration tests
- Passed mypy and ruff
