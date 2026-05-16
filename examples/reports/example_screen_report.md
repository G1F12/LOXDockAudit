# LOXDockAudit Screen Report: example_screen

## Executive Decision

**FAIL**

example_candidate does not advance. The candidate did not exceed the baseline productive pose count. Candidate does not advance.

## Criteria Evaluation

| Criterion | Result | Notes |
| --- | --- | --- |
| pass_count_vs_baseline | FAIL |  |
| best_distance_vs_baseline | PASS |  |
| best_productive_rank | PASS |  |
| cbd_contribution | SKIPPED | criterion disabled |
| beats_scrambled | FAIL |  |
| beats_polyK | FAIL |  |

## Candidate Summary

| construct_id | productive_count | total_poses | best_distance | best_productive_rank | lox_contact_frequency | cbd_contact_frequency | cbd_coupling | strict_pass |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| test_lox_baseline | 2 | 3 | 6.000 | 1 | 1.000 | 0.000 | 0.000 | False |

## Control Comparison

| construct_id | control_type | candidate_productive_count | control_productive_count | candidate_best_distance | control_best_distance | candidate_beats_control_count | candidate_beats_control_distance | conclusion |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| test_lox_baseline | baseline | 2 | 2 | 6.0 | 6.0 | False | False | Candidate matches baseline. |
| test_lox_baseline | scrambled | 2 | 2 | 6.0 | 6.0 | False | False | Candidate matches scrambled control. Binding may be nonspecific. |
| test_lox_baseline | polyK | 2 | 2 | 6.0 | 6.0 | False | False | Candidate matches polyK control. Binding may be nonspecific. |

## Per-Construct Pose Details

### example_baseline

| rank | active_site_to_target_distance | productive | lox_contacts | cbd_contacts |
| --- | --- | --- | --- | --- |
| 1 | 6.000 | True | 136 | 0 |
| 2 | 14.000 | False | 16 | 0 |
| 3 | 7.000 | True | 128 | 0 |

### example_candidate

| rank | active_site_to_target_distance | productive | lox_contacts | cbd_contacts |
| --- | --- | --- | --- | --- |
| 1 | 6.000 | True | 136 | 0 |
| 2 | 14.000 | False | 16 | 0 |
| 3 | 7.000 | True | 128 | 0 |

### example_scrambled

| rank | active_site_to_target_distance | productive | lox_contacts | cbd_contacts |
| --- | --- | --- | --- | --- |
| 1 | 6.000 | True | 136 | 0 |
| 2 | 14.000 | False | 16 | 0 |
| 3 | 7.000 | True | 128 | 0 |

### example_polyK

| rank | active_site_to_target_distance | productive | lox_contacts | cbd_contacts |
| --- | --- | --- | --- | --- |
| 1 | 6.000 | True | 136 | 0 |
| 2 | 14.000 | False | 16 | 0 |
| 3 | 7.000 | True | 128 | 0 |

## Warnings

No warnings.
