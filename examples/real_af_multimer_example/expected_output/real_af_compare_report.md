# HDOCK vs AF-Multimer Comparison: real_af_compare

## Executive Summary

The comparison shows strong convergence between productive HDOCK interfaces and the AF-Multimer interface under the current geometric metric.

This is a computational agreement check. It is consistent with a cross-method structural comparison, but it does not establish enzymatic activity, catalytic turnover, collagen oxidation, or biological efficacy.

## Productive Geometry

- HDOCK productive poses: 2 / 10
- HDOCK best productive rank: 3
- AF active-site distance: 6.159 A
- AF productive: True
- AF fully productive: False

## Interface Convergence

- Shared interface residues: 26
- Jaccard index: 0.667
- Category: strong convergence
- Narrative: strong convergence: the methods share a substantial interface residue set under the current contact metric

## Contact Persistence

| receptor_resnum | ligand_resnum | productive_pose_count | productive_pose_frequency | present_in_af | shared_contact |
| --- | --- | --- | --- | --- | --- |
| 4 | 191 | 1 | 0.5 | True | True |
| 6 | 335 | 1 | 0.5 | True | True |
| 7 | 334 | 1 | 0.5 | True | True |
| 7 | 335 | 1 | 0.5 | True | True |
| 7 | 336 | 1 | 0.5 | True | True |
| 8 | 328 | 1 | 0.5 | True | True |
| 8 | 329 | 1 | 0.5 | True | True |
| 8 | 330 | 1 | 0.5 | True | True |
| 8 | 336 | 1 | 0.5 | True | True |
| 9 | 328 | 1 | 0.5 | True | True |
| 9 | 329 | 1 | 0.5 | True | True |
| 9 | 336 | 1 | 0.5 | True | True |
| 9 | 341 | 1 | 0.5 | False | False |
| 10 | 268 | 1 | 0.5 | True | True |
| 10 | 329 | 1 | 0.5 | True | True |
| 10 | 330 | 1 | 0.5 | True | True |
| 10 | 331 | 1 | 0.5 | True | True |
| 10 | 336 | 1 | 0.5 | True | True |
| 10 | 338 | 1 | 0.5 | False | False |
| 10 | 339 | 1 | 0.5 | False | False |
| 10 | 340 | 1 | 0.5 | False | False |
| 10 | 341 | 1 | 0.5 | False | False |
| 10 | 342 | 1 | 0.5 | False | False |
| 11 | 221 | 1 | 0.5 | True | True |
| 11 | 268 | 1 | 0.5 | True | True |
| 11 | 326 | 1 | 0.5 | False | False |
| 11 | 329 | 1 | 0.5 | True | True |
| 11 | 341 | 1 | 0.5 | False | False |
| 11 | 354 | 1 | 0.5 | True | True |
| 12 | 221 | 1 | 0.5 | True | True |
| 12 | 326 | 1 | 0.5 | False | False |
| 12 | 336 | 1 | 0.5 | False | False |
| 13 | 220 | 1 | 0.5 | True | True |
| 13 | 328 | 1 | 0.5 | False | False |
| 13 | 329 | 1 | 0.5 | False | False |
| 13 | 336 | 1 | 0.5 | False | False |
| 13 | 354 | 1 | 0.5 | False | False |
| 14 | 220 | 1 | 0.5 | True | True |
| 14 | 221 | 2 | 1.0 | True | True |
| 14 | 224 | 1 | 0.5 | True | True |
| 14 | 225 | 1 | 0.5 | True | True |
| 14 | 266 | 1 | 0.5 | True | True |
| 14 | 326 | 1 | 0.5 | False | False |
| 14 | 328 | 1 | 0.5 | False | False |
| 14 | 354 | 1 | 0.5 | False | False |
| 14 | 355 | 1 | 0.5 | False | False |
| 14 | 356 | 1 | 0.5 | False | False |
| 14 | 359 | 1 | 0.5 | False | False |
| 15 | 221 | 1 | 0.5 | False | False |
| 15 | 224 | 1 | 0.5 | True | True |
| 16 | 221 | 1 | 0.5 | False | False |
| 17 | 221 | 1 | 0.5 | False | False |
| 17 | 225 | 1 | 0.5 | False | False |
| 17 | 266 | 1 | 0.5 | False | False |
| 18 | 225 | 1 | 0.5 | False | False |
| 20 | 224 | 1 | 0.5 | False | False |

## Interpretation Limits

AlphaFold-Multimer and docking have different objective functions and failure modes. Partial convergence is useful as independent computational support for an interface hypothesis, but divergent interfaces do not by themselves disprove docking and convergent interfaces do not prove mechanism.
