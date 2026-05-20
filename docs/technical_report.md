---
title: "LOXDockAudit: A Reproducible Framework for Control-Aware 
        Productive Geometry Analysis in LOX–Collagen Docking Screens"
author: "Sasha"
date: "May 2026"
version: "v0.5.1"
repository: "https://github.com/G1F12/LOXDockAudit"
---

# LOXDockAudit: A Reproducible Framework for Control-Aware
# Productive Geometry Analysis in LOX–Collagen Docking Screens

**Author:** Sasha  
**Date:** May 2026  
**Version:** v0.5.1  
**Repository:** https://github.com/G1F12/LOXDockAudit  

---

## Abstract

LOXDockAudit is a reproducible analysis framework for evaluating LOX-collagen docking screens through control-aware productive geometry metrics. It addresses the central limitation that docking score does not equal productive catalytic geometry: a pose that binds well does not necessarily place the LOX active site near the target collagen Lys or Hyl residue. The tool computes active-site-to-substrate distance metrics, compares results against controls, performs structural quality control, and records auditable inputs and deterministic analysis behavior. In the real Round 5 reproduction of LOX169-417 docked against collagen-I 6VZX, 2 of 10 poses satisfied the productive geometry criterion, with the best distance measuring 6.16 Å at docking rank 3. During development, lexicographic file ordering changed the pass/fail outcome; LOXDockAudit fixed this with natural sorting and a regression test.

## 1. Background

LOXDockAudit sits at the intersection of enzyme mechanism, extracellular matrix biology, and docking-screen reproducibility. Its central premise is that LOX engineering requires evaluating whether docking poses support catalysis, not only whether protein surfaces contact one another.

### 1.1 Lysyl oxidase and collagen crosslinking

Lysyl oxidase oxidizes lysine and hydroxylysine residues in collagen and elastin, initiating covalent crosslink formation in the extracellular matrix. These reactions produce hydroxylysyl pyridinoline and lysyl pyridinoline crosslinks that are critical for tissue tensile strength, matrix stability, and ECM mechanics. Productive LOX chemistry requires the copper cofactor, the lysyl tyrosylquinone cofactor, and correct substrate positioning. Engineering LOX for enhanced collagen crosslinking therefore requires designs that place the active site near target Lys or Hyl residues rather than merely increasing surface binding.

### 1.2 Why docking score alone is insufficient

Docking programs optimize and rank poses primarily by binding-energy proxies, shape complementarity, and contact patterns, not by enzyme-specific catalytic geometry. A top-ranked pose can bury a large interface while leaving the active site far from the reactive collagen side chain. In a LOX fusion construct, a collagen-binding domain can bind collagen effectively while trapping the catalytic domain in the wrong orientation. Contact count has the same limitation: it measures binding extent, not whether the docked complex supports oxidation of a target Lys or Hyl residue.

### 1.3 The productive geometry problem

For LOX catalysis, the active-site region must be positioned within approximately 8 Å of the target Lys or Hyl NZ atom. His292, His294, and His296 coordinate the catalytic copper, so substrate proximity to this region is required for a productive pose. This requirement combines distance and orientation constraints; it is not satisfied by binding affinity or interface size alone. Generic docking tools do not report this LOX-specific productive geometry criterion, which creates the need for a dedicated post-docking audit.

## 2. Problem Statement

LOX engineering screens require evidence that a docked complex supports catalysis. Generic docking outputs provide binding-oriented rankings, while LOXDockAudit evaluates whether those rankings also satisfy active-site geometry, control performance, structural plausibility, and input integrity.

### 2.1 What generic docking tools report

HDOCK and similar docking tools output binding scores and ranked pose models. Pose rank reflects predicted binding affinity, but it provides no information about active-site geometry relative to the target substrate residue. These tools also do not include built-in comparison to negative controls. File ordering is often assumed to be sorted by rank, but filesystem or lexical ordering does not guarantee that model files are processed in numeric rank order.

### 2.2 What is actually needed for LOX engineering

LOX engineering requires pre-specified criteria that connect docking output to catalytic interpretation:

- Active-site-to-Lys/Hyl distance <= 8 Å (productive geometry)
- Productive pose count > baseline (not just absolute count)
- First productive pose within top 3 ranks
- CBD/FMOD contact must not anti-correlate with active-site proximity
- Candidate must outperform scrambled and polyK controls
- Input files must be verified unique (SHA256)
- Structural fold must be plausible before docking is interpreted

### 2.3 Specific failure modes observed

Four concrete failure modes were observed during development and Round 5 review. CBD trapping occurred when the FMOD fusion produced contact_freq=1.0 but 0/10 productive poses, demonstrating binding without catalytic placement. Nonspecific binding occurred when scrambledSILY outperformed the SILY candidate, invalidating a sequence-specific interpretation. Lexicographic file ordering placed model_10 before model_2, changing best_productive_rank from 3 to 4. Identical input files were detected when a candidate and control had the same SHA256 hash, making the comparison scientifically invalid.

## 3. LOXDockAudit: Method and Design

LOXDockAudit implements a post-docking analysis workflow that treats docking poses as raw evidence rather than final conclusions. The design combines productive geometry measurement, explicit controls, structural quality checks, and input auditing into a reproducible decision framework.

### 3.1 Core concept: productive geometry vs binding score

LOXDockAudit separates binding signal from productive catalytic geometry. It evaluates docking poses not by score but by active-site proximity to the target substrate residue. A high-ranking pose therefore contributes to a positive result only when the LOX active-site region is positioned close enough to the collagen Lys or Hyl side chain to support catalytic interpretation.

### 3.2 Active-site-to-substrate distance metric

The primary metric is the minimum heavy-atom distance between atoms from His292, His294, His296, Lys320, and Tyr355 and the target Lys/Hyl NZ atom. The default threshold is 8.0 Å and is configurable for sensitivity analysis or alternative mechanistic assumptions. A productive pose is defined as any pose with distance <= threshold. LOXDockAudit reports the productive pose count, the best observed distance, and the first productive docking rank; in Round 5, the real LOX169-417 versus collagen-I 6VZX screen produced 2/10 productive poses, a best distance of 6.16 Å, and first productive rank 3.

### 3.3 Control-aware screen comparison

The screen configuration includes candidate, baseline, scrambled, polyK, noCBD, and inactive constructs. A strict criteria engine evaluates each pre-specified rule and reports it as PASS, FAIL, or SKIPPED. The workflow prevents post-hoc selection by defining all criteria before analysis and applying them consistently across constructs. The output includes a per-criterion verdict, an overall pass/fail result, and decision text that explains the conclusion in terms of productive geometry and control behavior.

### 3.4 Structural QC pipeline

The structural QC pipeline flags fold-corrupted models before docking interpretation. In the real Round 5 QC result, His triad CA geometry had a maximum distance of 9.158 Å and passed the 10 Å threshold. The Lys320-Tyr355 CB distance was 6.775 Å and passed the 12 Å threshold. Disulfide geometry passed with 5/5 intact disulfides, and the active-site burial proxy reported the site as accessible, yielding an overall PASS for interpreting the docking screen.

### 3.5 Input audit and reproducibility guarantees

The input audit records a SHA256 hash for every file used in the analysis. It performs cross-construct duplicate detection, N-terminus sequence checks, and chain ID verification. These checks prevent scientifically invalid comparisons caused by reused structures, mislabeled controls, truncated constructs, or chain mismatches. Together with deterministic sorting and recorded configuration, the audit makes the analysis reproducible and reviewable.

### 3.6 Substrate Orientation Score

LOXDockAudit supports an optional substrate orientation score for target Lys/Hyl residues. The score measures the angle between the approach vector from the active-site centroid to the target NZ atom and the side-chain vector from Lys C-beta to NZ. Smaller angles indicate that the side chain points toward the active site. When orientation scoring is enabled, a fully productive pose must satisfy both the distance threshold and the orientation threshold; when it is disabled, distance-only behavior is preserved.

### 3.7 Independent structural prediction comparison

Version v0.5.1 includes an AlphaFold-Multimer/ColabFold comparison workflow through `loxdockaudit af-compare`. The purpose is to compare productive HDOCK docking poses against an independently generated predicted complex for the same LOX construct and collagen target region. The workflow measures AF active-site distance, optional substrate orientation, HDOCK productive pose count, interface residue overlap, contact-pair persistence across productive HDOCK poses, and an interface convergence category.

The convergence categories are intentionally conservative. Strong convergence means substantial interface overlap under the configured contact metric, partial convergence means some shared interface residues, and divergent interfaces means no shared residues were detected. These are heuristic computational labels. They do not establish enzymatic activity, catalytic mechanism, collagen oxidation, crosslink formation, biomechanical improvement, therapeutic utility, or in vivo safety.

## 4. Software Implementation

LOXDockAudit is implemented as a tested Python package with a command-line interface, configuration-driven analysis, and report generation. The software design keeps scientific assumptions explicit in configuration files and separates parsing, measurement, QC, criteria evaluation, and reporting into focused modules.

### 4.1 Architecture overview

LOXDockAudit is a Python package using a `src` layout and `pyproject.toml`, with focused modules, a Click CLI, and a YAML-based configuration system. All analysis is driven by config files; no hardcoded paths or residue numbers appear in the core logic.

### 4.2 Key modules

| Module           | Responsibility                              |
|------------------|---------------------------------------------|
| pdb_parser.py    | PDB loading, atom selection, SHA256         |
| distances.py     | Active-site-to-substrate distance and pose geometry |
| orientation.py   | Target Lys/Hyl side-chain orientation score |
| contacts.py      | Ligand-receptor contact counting            |
| structural_qc.py | His triad, Lys-Tyr, disulfide, SASA proxy  |
| alphafold_qc.py  | pLDDT/PAE parsing from AF2/ColabFold output |
| fold_qc.py       | QC aggregation and fold-corrupted flag      |
| criteria.py      | Pass/fail criteria engine                   |
| screen.py        | Multi-construct screen orchestration        |
| reporting.py     | CSV, Markdown report generation             |
| af_parser.py     | AF-Multimer/ColabFold PDB and JSON parsing  |
| interface_extractor.py | Interface residue and contact fingerprint extraction |
| convergence.py   | HDOCK/AF interface overlap and convergence labels |
| af_geometry.py   | AF-Multimer productive-geometry scoring     |
| af_compare.py    | End-to-end HDOCK versus AF-Multimer workflow |

### 4.3 CLI interface

`loxdockaudit run` analyzes a single docking result and writes per-pose geometry outputs.

```bash
loxdockaudit run --config configs/round5_baseline.yaml --output results/round5_baseline
```

`loxdockaudit screen` evaluates a multi-construct screen with controls and produces criteria-level verdicts.

```bash
loxdockaudit screen --config configs/round5_screen.yaml --output results/round5_screen
```

`loxdockaudit check-config` validates a YAML configuration before analysis.

```bash
loxdockaudit check-config --config configs/round5_screen.yaml
```

`loxdockaudit af-compare` compares productive HDOCK poses against an AF-Multimer or ColabFold predicted complex.

```bash
loxdockaudit af-compare --config examples/real_af_multimer_example/r5_af_compare_real.yaml --out examples/real_af_multimer_example/expected_output
```

### 4.4 Test coverage

LOXDockAudit v0.5.1 extends the test suite with synthetic AF-Multimer fixtures. Unit tests cover `pdb_parser`, `distances`, `orientation`, `contacts`, `input_audit`, `config`, `reporting`, criteria/screen behavior, structural QC, inactive-control analysis, AF parsing, interface extraction, convergence logic, and Round 5 regression behavior. Integration tests cover the CLI workflows, control-aware screen analysis, inactive catalytic controls, AF comparison outputs, and backward compatibility with the v0.3 distance-only Round 5 metrics. Static checks are expected to report mypy: 0 type errors and ruff: 0 linting violations, and CI runs through a GitHub Actions workflow on push.

## 5. Results: Real Round 5 Reproduction

The Round 5 reproduction tested the full LOXDockAudit workflow on a real LOX-collagen docking result. The analysis combined structural QC, productive geometry scoring, domain-contact accounting, and deterministic rank handling.

### 5.1 Dataset

The construct was LOX169-417, human lysyl oxidase residues 169-417. The receptor was the collagen-I triple helix from PDB 6VZX, and docking was performed with HDOCK. The analysis included 10 docked models. The target residue was chain A, Lys11, atom NZ. The active-site residues were His292, His294, His296, Lys320, and Tyr355. The confirmed disulfide pairs were (238,291), (244,361), (324,340), (330,351), and (398,412).

### 5.2 Structural QC result

| Check              | Result | Value      |
|--------------------|--------|------------|
| His triad (max CA) | PASS   | 9.158 Å    |
| Lys320-Tyr355 CB   | PASS   | 6.775 Å    |
| Disulfides         | PASS   | 5/5 intact |
| Active-site access | PASS   | accessible |

His triad CA distance of 9.158 Å is near the upper boundary of the 10 Å threshold, consistent with computational model flexibility relative to crystal structures. All four checks passing indicates the LOX169-417 model is structurally adequate for docking interpretation.

### 5.3 Docking analysis result

| Metric                    | Value  |
|---------------------------|--------|
| Total poses analyzed      | 10     |
| Productive poses          | 2 / 10 |
| Best active-site distance | 6.16 Å |
| Best productive rank      | 3      |
| LOX contact frequency     | 1.000  |
| CBD contact frequency     | 0.000  |

The screen produced 2 productive poses out of 10. The best active-site distance was 6.16 Å, well within the 8 Å productive-geometry threshold, and the first productive pose occurred at rank 3, satisfying the pre-specified maximum rank criterion of 3. The LOX domain contacted collagen in all 10 poses, while CBD contact frequency was 0.000 because LOX169-417 is the baseline construct and contains no fusion CBD. This result established LOX169-417 alone as the baseline for future fusion construct comparisons.

### 5.4 Interpretation

The Round 5 reproduction showed that baseline LOX169-417 produced a measurable productive-geometry signal without a fusion domain. This baseline set the minimum comparison standard for engineered constructs: future candidates must improve productive pose count, preserve early-rank productive geometry, and outperform negative controls rather than only increasing collagen contact.

### 5.5 Active vs catalytically inactive control

An H292A/H294A/H296A inactive control was added for LOX169-417 to test whether
productive geometry metrics are sensitive to catalytic-site disruption. The
control was generated as a ColabFold/AlphaFold2-ptm single-chain mutant model and
docked against the same 6VZX lysine-site receptor using the paired Round 5
HDOCKlite setup. This is a computational structural-control experiment, not an
enzymatic-activity assay.

| Metric | Active | Inactive |
| --- | --- | --- |
| Productive poses | 2/10 | 1/10 |
| Best active-site distance | 6.159 A | 6.421 A |
| Best productive rank | 3 | 7 |
| fully_productive_count | 1 | 0 |
| fully_productive_fraction | 0.100 | 0.000 |
| best_orientation_angle_deg | 68.698 | 29.197 |
| best_fully_productive_rank | 5 | N/A |
| Productive fraction | 0.200 | 0.100 |
| QC status | PASS | PASS |
| Active-site geometry | intact | disrupted |

The active construct retained more favorable productive geometry under the
current docking metric. This supports the value of the inactive structural
control, while avoiding an overclaim: the result does not prove enzymatic
activity, copper loading, LTQ chemistry, or collagen oxidation.

## 6. Bug Discovery: Lexicographic Ordering

The strongest reproducibility finding during development was a rank-ordering bug caused by filename sorting. The bug changed a scientific verdict, not only a display order.

### 6.1 The bug

During development, PDB files were sorted using Python's default string sort. For files named model_1.pdb through model_10.pdb, this produced the order:

model_1, model_10, model_2, model_3, ..., model_9

The file model_10.pdb was assigned rank 2, and model_2.pdb was assigned rank 3.

### 6.2 Scientific impact

The best productive pose in the Round 5 LOX169-417 dataset is model_3.pdb, with active-site-to-target distance 6.159 Å.

Under lexicographic ordering, model_3.pdb received rank 4. The pre-specified criterion required first productive rank <= 3. Under lexicographic ordering, this construct would have received a FAIL verdict for the rank criterion.

Under correct numeric ordering, model_3.pdb receives rank 3, and the rank criterion passes.

A file-naming convention caused a pass/fail outcome to change. This was not detectable by inspecting docking scores alone.

### 6.3 Fix and regression test

The fix replaced `sorted()` with `sorted(..., key=natural_sort_key)`, where `natural_sort_key` splits filenames into numeric and text chunks. The change was applied in `screen.py`, `cli.py`, and `input_audit.py`.

The regression test `test_natural_sort_order()` verifies that `model_1, model_10, model_2, model_3, model_9` sorts to `model_1, model_2, model_3, model_9, model_10`. The fix changed the reported `best_productive_rank` from 4 to 3 on the real Round 5 dataset.

## 7. Limitations

1. Distance threshold of 8.0 Å is a geometric heuristic, not derived from LOX crystal structure kinetics data.
2. Orientation threshold of 90° is a geometric heuristic.
3. Active-site-to-substrate distance uses minimum heavy-atom distance and does not model docking pocket shape.
4. Structural QC uses CA/CB distances and a burial proxy, not true SASA or normal-mode analysis.
5. pLDDT/PAE parsing requires AlphaFold/ColabFold output; homology models are not supported.
6. CBD coupling metric (Pearson r over 10 poses) is statistically underpowered; interpret with caution.
7. AF-Multimer interface convergence is a computational agreement metric, not wet-lab validation and not proof of mechanism.
8. LOXDockAudit does not prove: enzymatic activity, Cu loading, LTQ/topaquinone maturation, collagen oxidation, crosslink formation, mechanical strengthening, or in vivo safety.
9. All analysis is based on static docking poses or static predicted complexes; molecular dynamics and conformational flexibility are not modeled.

## 8. Future Work

1. Expand example datasets: include a failed FMOD trapping case and a CNA-like positive benchmark.
2. PyMOL script generation for visualizing productive poses and active-site geometry.
3. Statistical comparison across controls using bootstrap resampling rather than deterministic criteria.
4. Support for AlphaFold3 output format (updated JSON schema).
5. Add curated real AF-Multimer example outputs once lightweight redistribution constraints are resolved.

## References

No external references are included in this technical portfolio report.
