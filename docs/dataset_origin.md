# Dataset Origin

## mouse_gse85

* **GEO:** [GSE85](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85)
* **Publication:** *Projection of an immunological self shadow within the thymus by the Aire protein*. Science, 2002. PMID: [12376594](https://pubmed.ncbi.nlm.nih.gov/12376594/)
* **Authors:** Mark S. Anderson, Emily S. Venanzi, Ludger Klein, Zhibin Chen, Stuart P. Berzins, Shannon J. Turley, Harald von Boehmer, Roderick Bronson, Andrée Dierich, Christophe Benoist, Diane Mathis.

### Samples

* **Comparison:** 3 control arrays vs 3 Aire-KO arrays.
* **Control:** GSM2444, GSM2445, GSM2446.
* **Aire-KO:** GSM2447, GSM2448, GSM2449.
* Each array was generated from thymic tissue pooled from five mice.
* Mice were 3–5 weeks old.
* Both male and female mice were used.
* Genetic background: B6/129 F2.
* Therefore, the array-level sample size is 3 vs 3, but each group represents pooled material from approximately 15 mice.

### Cell preparation

* Thymi were digested using collagenase and dispase.
* Stromal cells were enriched by density-gradient separation.
* mTECs were isolated by flow cytometry.
* Reported phenotype included:

  * CD45-negative
  * G8.8-positive
  * CDR1-intermediate
  * B7.1/CD80-high
* The samples therefore represent relatively mature mTECs rather than all thymic epithelial cells.

### RNA and microarray

* Total RNA was extracted from sorted cells.
* RNA underwent two rounds of T7 polymerase-based amplification.
* **Platform:** Affymetrix Murine Genome U74A Version 2 Array.
* **GEO platform:** GPL81.
* **Bioconductor CDF:** `mgu74av2cdf`.
* Array type: in situ synthesized oligonucleotide, classical Affymetrix expression array.
* GPL81 contains approximately 12,488 probe sets and is based on an old mouse transcript annotation, including UniGene Build 74.

### Data processing

* Raw CEL files are available.
* Original signal values were generated using Affymetrix Microarray Suite 5.0.
* GEO contains:

  * raw signal values
  * present/absent/marginal calls
  * normalized values stored as `NORM_VALUE`
* The original processed values are not equivalent to RMA-normalized values used in some of the newer datasets.

### Potential technical or biological issues

* Each array is a pool of five mice. Pooling reduces the effect of individual outliers but prevents estimation of mouse-to-mouse biological variability.
* Sex-specific variation cannot be assessed because male and female tissue was pooled or analyzed together.
* The B6/129 F2 genetic background differs from the backgrounds used in the other datasets.
* The old U74Av2 array has substantially lower transcript coverage than the Mouse 430 2.0 and Mouse Gene 1.0 ST arrays.
* Probe annotations are old. Some probes may now map to different genes, multiple genes, obsolete identifiers, or no current gene.
* Two rounds of T7 amplification can increase amplification variability, particularly for low-abundance transcripts.
* The U74Av2 platform is a 3′-biased expression-array design. RNA degradation and amplification can therefore strongly affect measured expression.
* MAS5-derived values should not be directly compared with RMA-normalized values from other datasets.
* Small array-level sample size: 3 vs 3.

---

## mouse_gse2585

* **GEO:** [GSE2585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2585)
* **Publication:** *Promiscuous gene expression in thymic epithelial cells is regulated at multiple levels*. Journal of Experimental Medicine, 2005.
* **Authors:** Jens Derbinski, Jana Gäbler, Benedikt Brors, Sascha Tierling, Sunitha Jonnakuty, Manfred Hergenhahn, Leena Peltonen, Jörn Walter, Bruno Kyewski.

### Samples

* **Selected comparison:** 2 control arrays vs 2 Aire-KO arrays.
* **Control:** GSM49735, GSM49736.
* **Aire-KO:** GSM49731, GSM49732.
* These samples represent the CD80-high mTEC population.
* The full GEO series contains 16 samples, including additional epithelial-cell subsets.
* The selected comparison is therefore **2 vs 2**, not 3 vs 3.
* Mice were 6–12 weeks old.
* Wild-type mice were C57BL/6.
* Aire-KO mice were reported as having a mixed genetic background.
* Samples were obtained from two independent experiments.

### Cell preparation

* Thymic stromal cells were isolated using sequential:

  * collagenase digestion
  * collagenase/dispase digestion
  * trypsin digestion
  * density-gradient centrifugation
* mTECs were identified approximately as:

  * CD45-negative or CD45-low
  * CDR1/Ly51-negative
  * EpCAM-positive
* mTECs were further separated according to CD80 expression.
* The selected samples are CD80-high mature mTECs.

### RNA and microarray

* RNA was pre-amplified and biotin-labelled.
* The protocol used two rounds of cDNA synthesis and in vitro transcription.
* **Platform:** Affymetrix Murine Genome U74A Version 2 Array.
* **GEO platform:** GPL81.
* **Bioconductor CDF:** `mgu74av2cdf`.
* This is the same microarray platform as GSE85.

### Data processing

* Raw CEL files are available.
* Fluorescence intensity was evaluated using Affymetrix Microarray Suite 5.0.
* GEO processed data include:

  * MAS5 expression values
  * absolute detection calls
  * detection p-values

### Potential technical or biological issues

* Very small selected sample size: 2 vs 2.
* Wild-type and Aire-KO mice do not have the same reported genetic background. Genotype effects may therefore be partly confounded by strain effects.
* The selected samples are only the CD80-high subset and should not be interpreted as total mTEC expression.
* The full series contains multiple cell types and maturation states. Sample selection must be checked carefully.
* The U74Av2 platform has limited transcript coverage and old probe annotations.
* Two rounds of amplification may increase noise or bias in low-expression genes.
* The platform has a 3′-biased probe design and may be sensitive to RNA quality.
* Original MAS5 values should not be compared directly with RMA values from GSE14365 or GSE33878.
* Although GSE85 and GSE2585 use the same platform, they differ in:

  * mouse age
  * genetic background
  * cell-isolation protocol
  * mTEC gating
  * pooling strategy
  * laboratory and experiment batch

---

## mouse_gse14365

* **GEO:** [GSE14365](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14365)
* **Publication:** *Aire-deficient C57BL/6 mice mimicking the common human 13-base pair deletion mutation present with only a mild autoimmune phenotype*. Journal of Immunology, 2009. PMID: [19265170](https://pubmed.ncbi.nlm.nih.gov/19265170/)
* **Authors:** François-Xavier Hubert, Sarah A. Kinkel, Pauline E. Crewther, Ping Z.F. Cannon, Kylie E. Webster, Maire Link, Raivo Uibo, Moira K. O’Bryan, Anthony Meager, Simon P. Forehan, Gordon K. Smyth, Lauréane Mittaz, Stylianos E. Antonarakis, Pärt Peterson, William R. Heath, Hamish S. Scott.

### Samples

* **Selected comparison:** 3 control arrays vs 2 Aire-mutant arrays.
* **Control:** GSM359245, GSM359246, GSM359247.
* **Aire mutant:** GSM359248, GSM359249.
* The cell population was reported as approximately:

  * CD45-negative
  * MHC class II-high
  * Ly51-low
* These cells represent a mature mTEC-enriched population.

### Aire model

* The mice carried an Aire exon 8 mutation designed to mimic a common human 13-base-pair deletion.
* The model is on a C57BL/6 background.
* The mutation produced a truncated Aire–β-galactosidase fusion product.
* This model may not be biologically identical to a complete Aire-null allele used in other studies.

### RNA and microarray

* Total RNA was isolated using the RNeasy Micro kit.
* Approximately 100 ng total RNA was used.
* RNA amplification used:

  * T7-oligo(dT) priming
  * MEGAscript T7 amplification
  * second-strand cDNA synthesis
  * Affymetrix GeneChip in vitro transcription and biotin-labelling
* Approximately 15 µg fragmented labelled material was hybridized to each array.
* **Platform:** Affymetrix Mouse Genome 430 2.0 Array.
* **GEO platform:** GPL1261.
* **Bioconductor CDF:** `mouse4302cdf`.
* Arrays were hybridized overnight and scanned using an Affymetrix GeneChip Scanner 3000.

### Platform annotation

* GPL1261 is an in situ oligonucleotide array.
* Probe design used sequences from GenBank, dbEST, RefSeq, UniGene Build 107 and the mouse-genome resources available at the time.
* The platform provides broader gene coverage than GPL81 but still uses older probe-set definitions and a classical 3′-expression design.

### Data processing

* Expression values were normalized and summarized using RMA.
* Differential expression was analyzed using `limma`.
* Linear models and empirical Bayes statistics were used.
* Multiple-testing correction used the Benjamini–Hochberg method.
* The publication used an FDR threshold of 0.05.
* Three additional closely related samples were included when fitting the original linear model, although they were not part of the reported wild-type-versus-mutant t-statistic comparison.

### Potential technical or biological issues

* Small and unbalanced selected sample size: 3 vs 2.
* Reanalysis using only the five selected arrays may not exactly reproduce the published statistics because the original model included three additional related samples.
* The mutation is not necessarily equivalent to the complete Aire-knockout alleles used in GSE85, GSE2585 or GSE33878.
* The sorting definition differs from the CD80-high sorting used in GSE85 and GSE2585.
* Starting with only 100 ng RNA and applying T7-based amplification may introduce amplification variability.
* Mouse Genome 430 2.0 is a 3′-expression array and may be affected by RNA degradation or differences in transcript 3′ ends.
* GPL1261 probe-to-gene mappings are based on old transcript and genome annotations.
* Probe-level values cannot be directly combined with GPL81 or GPL6246 data.
* RMA-normalized values cannot be directly compared with the original MAS5 values from GSE85 or GSE2585.

---

## mouse_gse33878

* **GEO:** [GSE33878](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33878)
* **Publication:** *Aire unleashes stalled RNA polymerase to induce ectopic gene expression in thymic epithelial cells*. Proceedings of the National Academy of Sciences, 2012.
* **Authors:** Matthieu Giraud, Hideyuki Yoshida, Jakub Abramson, Peter B. Rahl, Richard A. Young, Diane Mathis, Christophe Benoist.

### Samples

* **Comparison:** 3 control arrays vs 3 Aire-KO arrays.
* **Control:** GSM838272, GSM838273, GSM838274.
* **Aire-KO:** GSM838275, GSM838276, GSM838277.
* Samples were obtained from individual mice rather than pooled groups.
* Mice were 4–6 weeks old.
* Genetic background: `(B6 × NOD)F1`.
* Cells were described as MHC class II-high medullary thymic epithelial cells.

### Cell preparation

* Thymic epithelial cells were isolated and sorted using an Aria flow cytometer.
* The reported phenotype included:

  * CD45-negative
  * Ly51-intermediate
  * MHC class II-high
* Cells were sorted directly into TRIzol.
* This gating strategy differs from the CD80-high and Ly51-low definitions used in the other studies.

### RNA and microarray

* RNA was isolated using TRIzol.
* RNA underwent two amplification rounds using the MessageAmp aRNA kit.
* Biotin-labelled cRNA was generated using a standard Affymetrix protocol.
* **Platform:** Affymetrix Mouse Gene 1.0 ST Array.
* **GEO platform:** GPL6246.
* **Bioconductor platform package:** `pd.mogene.1.0.st.v1`.
* This is a whole-transcript array rather than the classical 3′-biased array design used for GPL81 and GPL1261.

### Hybridization and scanning

* Hybridization was performed for 16 hours at 45°C and 60 rpm.
* Instrument: GeneChip Hybridization Oven 640.
* Washing and staining were performed using an Affymetrix Fluidics Station 400.
* Staining included streptavidin R-phycoerythrin.
* Protocol: EukGE-WS2v4.
* Arrays were scanned with a GeneChip Scanner 3000.
* GeneChip Operating Software was used during scanning.

### Data processing

* Raw CEL files are available.
* Expression data were preprocessed using RMA.
* The study was designed partly to investigate variation in exon usage.
* The GEO gene-level table contains approximately 35,512 transcript-level rows.

### Potential technical or biological issues

* The `(B6 × NOD)F1` genetic background differs strongly from the backgrounds used in the other three studies.
* Baseline expression differences may therefore reflect strain as well as experimental or platform differences.
* The whole-transcript probe architecture is fundamentally different from the GPL81 and GPL1261 3′-expression arrays.
* Whole-transcript arrays can detect transcript-wide signal, but their expression scale and probe summarization are not directly comparable with older Affymetrix arrays.
* Two rounds of MessageAmp amplification may alter relative transcript-region representation.
* The sorting gate differs from the other datasets and may capture a somewhat different mTEC population.
* The publication investigated exon usage. Gene-level processing with `pd.mogene.1.0.st.v1` and a gene-level GPL6246 annotation table may discard exon- or probeset-level information relevant to that analysis.
* Exon-level analysis requires an appropriate probeset-level annotation and summarization workflow.
* Sample size remains small: 3 vs 3.

---

## rat_mtec_high

* **Dataset:** Internal rat mTEC-high dataset.

### Samples

* **Comparison:** 3 control arrays vs 3 Aire-KO arrays.
* **Control:** ctrl1, ctrl2, ctrl3.
* **Aire-KO:** ko1, ko2, ko3.
* **Age:** 4 weeks.
* **Species:** Rat.
* **Cell population:** mTEC-high.
* **Analysis comparison label:** HE vs KO.

### RNA and microarray

* **Platform:** Applied Biosystems/Thermo Fisher Clariom S Rat HT Array.
* **Array type:** `Clariom_S_Rat_HT`.
* **Bioconductor platform package:** `pd.clariom.s.rat.ht`.
* **Analysis type:** Gene-level expression.
* Clariom S is a whole-transcript expression array with probes distributed across transcripts rather than restricted mainly to the 3′ end.

### Data processing

* **Analysis version:** version 2.
* **Summarization method:** Gene Level – SST-RMA.
* **Positive-vs-negative AUC threshold:** 0.7.
* **Genome version:** rn6, *Rattus norvegicus*.
* **Annotation file:** `Clariom_S_Rat_HT.r1.na36.rn6.a1.transcript.csv`.
* Data processing pipeline: `rat`.

### Important technical points

* The rat dataset uses a newer array and annotation system than all four mouse datasets.
* SST-RMA processing is not identical to the RMA or MAS5 processing used for the mouse datasets.
* Rat and mouse expression values should not be compared directly on the raw or normalized intensity scale.
* Clariom S Rat HT and Mouse Gene 1.0 ST are both whole-transcript arrays, but they use different probe designs, annotation versions and summarization pipelines.
* The annotation is based on rn6 and NetAffx release 36. Gene identifiers should be mapped to current identifiers before cross-species analysis.
* Differences between rat and mouse may result from both species biology and platform-specific measurement effects.

---

## Cross-dataset technical comparison

### Sample sizes

| Dataset       | Control | Aire-deficient | Important replicate detail                                   |
| ------------- | ------: | -------------: | ------------------------------------------------------------ |
| GSE85         |       3 |              3 | Each array pools five mice                                   |
| GSE2585       |       2 |              2 | Two independent experiments; selected CD80-high samples      |
| GSE14365      |       3 |              2 | Unbalanced; original model included three additional samples |
| GSE33878      |       3 |              3 | Individual mice                                              |
| rat_mtec_high |       3 |              3 | Internal rat dataset; individual samples                     |

### Platforms

| Dataset       | Platform         | Array                | Probe design              | Processing                   |
| ------------- | ---------------- | -------------------- | ------------------------- | ---------------------------- |
| GSE85         | GPL81            | Mouse U74Av2         | Older 3′ expression array | MAS5 / GEO normalized values |
| GSE2585       | GPL81            | Mouse U74Av2         | Older 3′ expression array | MAS5                         |
| GSE14365      | GPL1261          | Mouse Genome 430 2.0 | 3′ expression array       | RMA                          |
| GSE33878      | GPL6246          | Mouse Gene 1.0 ST    | Whole-transcript array    | RMA                          |
| rat_mtec_high | Clariom S Rat HT | Clariom S Rat HT     | Whole-transcript array    | Gene-level SST-RMA           |
