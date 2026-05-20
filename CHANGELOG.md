# Changelog

## v0.5.1 - Realistic AF Example and Visual Summary

- Added a lightweight derived AF-Multimer-style comparison example for LOX169-417 versus the collagen-I 6VZX lysine-site target region.
- Added expected `af-compare` outputs for the example, including summary CSV, interface-overlap CSV, contact-persistence CSV, Markdown report, and a PNG visual summary.
- Added `docs/visual_summary.md` and a GitHub-readable workflow overview figure for portfolio and mentor review.
- Standardized AF comparison Markdown output as `{comparison_id}_report.md`.
- Updated documentation with the runnable real-example command and conservative interpretation language.
- Added `loxdockaudit af-compare` for HDOCK versus AlphaFold-Multimer/ColabFold cross-method comparison.
- Added AF/ColabFold PDB and JSON parsing helpers.
- Added interface residue extraction and residue-pair contact fingerprinting for HDOCK and AF structures.
- Added productive-geometry scoring for AF-Multimer models using the existing distance and orientation metrics.
- Added interface overlap, contact persistence, and conservative convergence categories.
- Added CSV and Markdown outputs for AF comparison workflows.
- Added synthetic unit and integration tests for AF parsing, interface overlap, orientation comparison, convergence categorization, and CLI output.

## v0.4.0 - Orientation and Inactive-Control Release

- Added optional substrate orientation scoring for target Lys/Hyl side-chain geometry.
- Added fully productive pose scoring when both distance and orientation criteria pass.
- Added catalytically inactive H292A/H294A/H296A Round 5 control workflow.
- Added pure-Python SVG geometry scatter plots and two-panel distance/orientation histograms.
- Added orientation-aware CSV and Markdown reporting fields.
- Added v0.4 regression tests on the real Round 5 active/inactive dataset.
- Preserved v0.3 distance-only metrics when orientation is disabled.
- Verified release QA: 215 tests passing, ruff passing, mypy passing.

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
