# Technical Report

## Background

LOXDockAudit evaluates whether LOX-collagen docking poses satisfy productive geometric constraints rather than relying on docking rank alone.
For lysyl oxidase systems, a binding-like interface is not sufficient: the LOX active-site region must be positioned near the target collagen Lys/Hyl residue with a plausible approach orientation.
The software therefore treats docking and AF-Multimer models as computational evidence that must be audited with explicit thresholds, controls, and reproducibility checks.

## Data

The primary repository example uses LOX169-417 and the collagen-I 6VZX lysine-site target region.
HDOCK pose inputs are stored under `examples/real_hdock_round5/cleaned_models/`.
The derived AF-Multimer-style example is stored under `examples/real_af_multimer_example/ranked_0.pdb`.
It is a lightweight redistributed example for demonstrating the comparison workflow, not a full raw ColabFold run directory.

Input target definition:

- Target residue: chain A, Lys11, atom NZ.
- LOX active-site chain: D.
- Active-site residues: His292, His294, His296, Lys320, Tyr355.
- Productive distance threshold: 8.0 angstrom.
- Orientation threshold: 90 degrees when enabled.
- Interface contact threshold: 5.0 angstrom for AF/HDOCK comparison.

## Methods

Structures are parsed with Biopython after lightweight PDB cleanup for HDOCK-style files.
Productive geometry is measured as the minimum heavy-atom distance between configured LOX active-site residues and the target substrate atom.
Optional orientation scoring measures the angle between the active-site-to-target approach vector and the Lys/Hyl CB-to-NZ side-chain vector.
When orientation is enabled, a fully productive pose must satisfy both distance and orientation thresholds.

Control-aware analysis compares active LOX169-417 against the H292A/H294A/H296A inactive control using the same distance and orientation criteria.
AF-Multimer comparison extracts interface residues and residue-pair contacts from productive HDOCK poses and from the AF-style ranked model.
Interface convergence is summarized with overlap count, Jaccard index, and conservative labels: strong convergence, partial convergence, or divergent interfaces.

## Results

Round 5 active LOX169-417 has 2 productive poses out of 10 and best productive rank 3.
The inactive H292A/H294A/H296A control has 1 productive pose out of 10 and best productive rank 7.
Orientation-aware scoring gives the active construct one fully productive pose and the inactive control zero fully productive poses.
Both active and inactive examples pass the configured structural QC checks.

The derived AF-Multimer-style comparison reports 26 shared interface residues between productive HDOCK poses and the AF-style model.
The Jaccard index is 0.667, producing a strong convergence label under the configured contact metric.
This supports computational interface consistency only and does not establish enzymatic activity.

## Limitations

Distance and orientation thresholds are geometric heuristics.
Docking and AF-Multimer models are static structures and do not model conformational dynamics.
The AF-Multimer-style example is lightweight and redistributed for workflow demonstration, not a complete prediction archive.
LOXDockAudit does not prove copper loading, LTQ chemistry, catalytic turnover, collagen oxidation, crosslink formation, therapeutic utility, or in vivo safety.
