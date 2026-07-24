Median expression values for the same 22 chemokines:

| Dataset        | Median control log₂ expression | SD across chemokines |
| -------------- | -----------------------------: | -------------------: |
| mouse_gse85    |                           4.88 |                 2.18 |
| mouse_gse14365 |                           6.99 |                 2.29 |
| mouse_gse33878 |                           7.82 |                 2.26 |
| mouse_gse2585  |                           9.21 |                 2.45 |
| rat_mtec_high  |                          10.46 |                 3.03 |


The rat dataset is:
- shifted toward higher raw values;
- somewhat more dispersed;
- much wider in total range.

Percentiles:
The four mouse datasets agree reasonably well:
- mouse–mouse Spearman correlations: 0.75–0.88
- average mouse–mouse correlation: approximately 0.80

The rat dataset agrees less strongly with each mouse dataset:
- rat–mouse correlations: 0.44–0.58
- average rat–mouse correlation: approximately 0.51


The higher absolute rat expression baseline is probably substantially technical because the raw mouse datasets also have
very different median expression levels. Different platforms and preprocessing pipelines produce different log₂ scales.



DEG:
mouse_gse85 up: 57
mouse_gse85 down: 71
mouse_gse85 total: 128
mouse_gse2585 up: 10
mouse_gse2585 down: 101
mouse_gse2585 total: 111
mouse_gse14365 up: 1
mouse_gse14365 down: 3
mouse_gse14365 total: 4
mouse_gse33878 up: 72
mouse_gse33878 down: 356
mouse_gse33878 total: 428
rat_mtec_high up: 172
rat_mtec_high down: 308
rat_mtec_high total: 480