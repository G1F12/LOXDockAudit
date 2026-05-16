# LOXDockAudit: A Reproducible Framework for Productive Geometry Analysis in LOX-Collagen Docking Screens

## 1. Problem: Docking Score Alone Is Not Enough

Protein-protein docking screens commonly rank candidate poses by a docking score or
energy-like value. That ranking is useful for triage, but it does not answer the
core biochemical question for lysyl oxidase (LOX)-collagen docking: whether a
pose places the LOX catalytic site in a geometry that could plausibly reach a
collagen lysine or hydroxylysine target.

For LOX-collagen systems, a high-scoring docked complex can still be
non-productive if the target Lys/Hyl side chains are far from the active-site
region or oriented away from the catalytic surface. Conversely, a lower-ranked
pose may be more mechanistically relevant if it places the active-site region
near a chemically plausible substrate residue. LOXDockAudit therefore treats
docking score as an input to audit, not as the final decision criterion.

The goal of the framework is reproducibility: given a set of docked model files,
a declared active-site definition, and target residue rules, the same productive
geometry summary should be reproduced by another user without hidden manual
steps.

## 2. Method: Active-Site-to-Lys/Hyl Productive Geometry

LOXDockAudit evaluates each docked pose by measuring the distance from a declared
active-site anchor on LOX to candidate collagen Lys/Hyl residues. A pose is
classified as productive when it satisfies the configured productive geometry
criteria, including the relevant chain, residue identity, and distance threshold.

This design makes the audit explicit. The configuration file records which chain
is treated as enzyme, which chain or chains are treated as substrate, and which
residue types are eligible substrate targets. The output reports per-pose
measurements as well as aggregate summaries such as the number of productive
poses, the best productive rank, and the closest observed active-site-to-target
distance.

The framework is intentionally conservative. It does not attempt to infer
enzymatic activity directly from a docked structure. It only reports whether a
docked pose satisfies a reproducible productive-geometry rule.

## 3. Controls: Baseline, Scrambled, PolyK, Inactive, and Benchmark Logic

Docking audits are most useful when candidate results are compared against
controls. LOXDockAudit supports screen-style workflows that can compare candidate
complexes with baseline or negative-control sets. Typical control logic includes:

- Baseline docking sets that establish the normal productive-pose frequency for a
  chosen receptor and collagen input.
- Scrambled or sequence-altered controls that test whether productive geometry is
  enriched beyond a nonspecific substrate-like surface.
- PolyK controls that stress-test whether the method overcalls productivity when
  many lysines are available.
- Inactive or active-site-disrupted controls that separate docking proximity from
  catalytic plausibility.
- Benchmark logic that keeps the same scoring, sorting, and reporting rules
  across all candidate and control groups.

The important engineering requirement is that candidate and control sets pass
through the same parser, sorting, filtering, and reporting path. That prevents
manual ranking differences or file-order artifacts from changing the conclusion.

## 4. Engineering: Preprocessing, Natural Numeric Sorting, and Audit Checks

The v0.2.1 release hardens the reproducibility layer around real HDOCK outputs.
Real docking files often contain formatting details that are valid enough for a
viewer but brittle for downstream parsers. LOXDockAudit includes preprocessing
utilities for HDOCK-style PDB files so malformed or wrapper records can be cleaned
before measurement.

Natural numeric sorting is also required for reproducibility. Lexicographic file
ordering places `model_10.pdb` before `model_2.pdb`, which silently corrupts rank
interpretation. LOXDockAudit sorts model filenames numerically where possible, so
`model_1.pdb`, `model_2.pdb`, and `model_10.pdb` retain the intended docking
rank order.

Input audit checks are part of the workflow. Configuration validation,
chain/residue presence checks, parser tests, and reporting tests reduce the
chance that a run succeeds while measuring the wrong target. The release also
keeps expected output files for the Round 5 example, allowing reproduced output
to be compared against a known result.

## 5. Real Example: Round 5 Reproduction

The v0.2.1 reproducibility example is the Round 5 LOX169-417 versus collagen-I
6VZX HDOCK dataset. The release notes define the expected reproduction command:

```bash
pip install -e .
loxdockaudit check-config examples/real_hdock_round5/config.yaml
loxdockaudit run --config examples/real_hdock_round5/config.yaml --out examples/real_hdock_round5/expected_output --top-n 10
```

For this dataset, the expected summary is:

- Productive poses: 2/10
- Best distance: 6.16 A
- Best productive rank: 3

These numbers should be interpreted as a reproducible geometry audit of the
docking output, not as a claim that the docked complex has experimentally
validated catalytic activity.

## 6. Limitations

LOXDockAudit does not prove enzymatic activity, substrate turnover, collagen
crosslink formation, or tendon strengthening. It only reports whether docked
poses satisfy declared proximity-based productive-geometry criteria.

The method also depends on the quality of the input structures, docking protocol,
active-site definition, and residue mapping. A productive pose can be a false
positive if the docking model is physically unrealistic. A non-productive pose can
be a false negative if the docking search missed a relevant conformation or if
the active-site and substrate definitions are too restrictive.

For publication or decision-making use, LOXDockAudit results should be paired
with orthogonal evidence: biochemical assays, mutational controls, molecular
dynamics or refinement, collagen-context validation, and clear reporting of all
configuration files and input structures.
