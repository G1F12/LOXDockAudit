# Active vs Catalytically Inactive Control Comparison

## Design

The inactive control was introduced to test whether productive geometry metrics were sensitive to catalytic-site disruption.

- Active construct: r5_real_lox169
- Inactive construct: r5_inactive_h292a_h294a_h296a
- Mutations: H292A, H294A, H296A
- Interpretation scope: computational structural-control experiment, not enzymatic activity validation.

## Results

| Metric | Active | Inactive |
| --- | --- | --- |
| Productive poses | 2/10 | 1/10 |
| Best distance | 6.159 | 6.421 |
| Best productive rank | 3 | 7 |
| Productive fraction | 0.200 | 0.100 |
| QC status | PASS | PASS |
| Active-site access | PASS | PASS |
| Active-site geometry | intact | disrupted |

## Interpretation

Active LOX169-417 retained more favorable productive geometry under the current docking metric. This is not proof of enzymatic activity.

## Limitations

The comparison uses docking-derived productive geometry and simple structural QC. HDOCK scoring does not include copper coordination, LTQ chemistry, catalytic turnover, or wet-lab activity.
