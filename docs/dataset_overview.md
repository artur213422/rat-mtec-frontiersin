# Dataset Overview

## mouse_gse14365 is outlier
mouse_gse14365 is a clear outlier. It has only four DEG, and its chemokine and 
kallikrein expression patterns do not match the other mouse datasets or the rat dataset.

The study used an Aire mutation designed to mimic a human mutation. This model produced
a weak Aire phenotype and is not equivalent to a conventional Aire knockout.

Because the biological model and expression response differ substantially from the other
datasets, mouse_gse14365 is excluded from further comparative analysis.



## Expression level and dataset stability

`mouse_gse14365` was excluded.

### Median expression
Rat expression is not consistently higher than all mouse datasets:

| Dataset        | Combined median expression |
| -------------- | -------------------------: |
| mouse_gse85    |                      3.859 |
| mouse_gse2585  |                      8.327 |
| mouse_gse33878 |                      7.022 |
| rat_mtec_high  |                      7.668 |


The rat dataset also has a wider expression range:

| Dataset        | Expression Q25–Q75 |   IQR |
| -------------- | -----------------: | ----: |
| mouse_gse85    |        2.818–5.139 | 2.321 |
| mouse_gse2585  |        6.985–9.616 | 2.632 |
| mouse_gse33878 |        5.725–8.279 | 2.554 |
| rat_mtec_high  |       5.119–10.567 | 5.447 |

The rat expression IQR is approximately twice as wide as in the mouse datasets.

### Gene-level variation

| Dataset        | Combined median gene SD | Combined median gene CV |
| -------------- | ----------------------: | ----------------------: |
| mouse_gse85    |                   0.173 |                  0.0470 |
| mouse_gse2585  |                   0.188 |                  0.0232 |
| mouse_gse33878 |                   0.139 |                  0.0210 |
| rat_mtec_high  |                   0.366 |                  0.0530 |

Rat has the highest gene-level SD:

* Approximately 2.1 times higher than `mouse_gse85`.
* Approximately 1.9 times higher than `mouse_gse2585`.
* Approximately 2.6 times higher than `mouse_gse33878`.

Only about 2.26% of rat genes are classified as DEG, so DEG alone cannot strongly raise the median SD across all genes.

Rat also has the highest median CV, although it is only slightly higher than `mouse_gse85`.

Within the rat dataset, control and KO variability are nearly identical:

| Statistic      | Control |     KO |
| -------------- | ------: | -----: |
| Mean gene SD   |  0.3790 | 0.3795 |
| Median gene SD |  0.3182 | 0.3196 |
| Mean gene CV   |  0.0551 | 0.0552 |
| Median gene CV |  0.0438 | 0.0443 |

There is no evidence that one rat group is substantially noisier than the other.

### Replicate correlations

| Dataset        | Mean control correlation | Mean KO correlation |
| -------------- | -----------------------: | ------------------: |
| mouse_gse85    |                    0.970 |               0.984 |
| mouse_gse2585  |                    0.992 |               0.977 |
| mouse_gse33878 |                    0.994 |               0.996 |
| rat_mtec_high  |                    0.980 |               0.980 |

All datasets have high within-group replicate correlations.

The average mouse correlation is approximately 0.986 for both control and KO samples. 
Rat correlations are approximately 0.980, slightly lower than the mouse average but still high.

Rat minimum correlations are also high:

* Control minimum correlation: 0.978.
* KO minimum correlation: 0.979.

Therefore, the rat dataset has a wider expression range and higher gene-level SD, but its samples remain strongly correlated. 
The higher SD does not indicate failed or unstable replicates.

### Conclusion

Rat does not have uniformly higher median expression than all mouse datasets. Its median lies within the range of the mouse datasets.

The rat platform produces a substantially wider expression distribution and higher gene-level SD. This is likely at least partly technical and related to the different platform and summarization method.

Despite the higher SD, rat control and KO samples are internally balanced and have correlations around 0.98. The rat dataset is therefore reasonably stable, although less tightly correlated than `mouse_gse33878`.
