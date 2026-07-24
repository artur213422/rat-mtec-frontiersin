

```
article/
├── data/
│   ├── raw/
│   │   ├── mouse/
│   │   └── rat/
│   │
│   ├── intermediate/
│   ├── processed/
│   └── exploration/
│       ├── chemokines/
│       ├── expression_patterns/
│       └── mtec_species_comparison/
├── src/
│   ├── data_preprocessing/
│   │    ├── pipelines
│   │    │    ├── mouse_dataset_pipeline/
│   │    │    └── rat_dataset_pipeline/
│   │    └── internal/
│   ├── common/
│   │    ├── types/
│   │    └── models/
│   ├── exploration/
│   └── manuscript/
│
├── figures/
│   ├── exploratory/
│   └── manuscript/
│
└── config/
```

rat annotation tables.
https://download.rgd.mcw.edu/data_release

https://www.ensembl.org/biomart/martview/4b392f59783586f3d35d436d27b689ff

- [ ] All down regulated genes mouse vs rat 
  - How much changed
  - Top genes

- [x] Chemokines mouse vs rat. Highlight thymus important genes (ccl19, ccl21, ccl25, xcl1, cxcl12, ccl17, ccl22).

- [ ] mTEC (and also skin) differentiation related genes (like epidermal differentiation complexi (EDC) genes).

- [ ] T1 and T2 ISG rat vs rot. Total ISG comparison is already done.

- [x] Kallikreins cluster expression pattern

- [x] Major Histocompatibility Complex (MHC) & Antigen Presentation comparison

- [x] TSA comparison

- [ ] Imprinted gene clusters

- [ ] Thymic mimetic cell marker genes

- [x] Peptide hormones

- [ ] Põhimõtteliselt võiks veel vaadata ka nt. kõrge vs madalaekspressioonitasemega geenid (st et basaalekspressioon oleks madal)
või ka nt suured vs väikesed geenid. aga see eeldab muidugi jälle kvaliteetse nimekirja olemasolu.

  interfeone



GSE85     = MAS5 / Norm signal; not explicitly log2; treat as linear signal
GSE2585   = MAS5 signal; not log2; treat as linear signal
GSE14365  = RMA in R; treat as log2 RMA expression
GSE33878  = RMA signal; scale not explicitly stated in GEO; treat as ambiguous