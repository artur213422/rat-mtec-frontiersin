GSE2585 already has mtec low replicas, but its 2 vs 2.


| Dataset                   | Control vs Aire-KO | Assay                             | mTEC-low definition        | Assessment                                       |
| ------------------------- | -----------------: | --------------------------------- | -------------------------- | ------------------------------------------------ |
| **GSE222285 / GSE205410** |         **8 vs 4** | Bulk RNA-seq, Illumina HiSeq 2500 | CD45⁻EpCAM⁺Ly51⁻MHC-II-low | Best dataset; duplicate GEO depositions          |
| **GSE65617**              |         **3 vs 3** | Bulk RNA-seq, Illumina HiSeq 2000 | Explicit `mTEClo` sorting  | Good independent dataset                         |
| **GSE221114**             |         **3 vs 3** | Bulk RNA-seq, Ion Torrent Proton  | Explicit `mTEClo` label    | Usable, but metadata are incomplete/inconsistent |
| **GSE2585**               |         **2 vs 2** | Affymetrix U74Av2 microarray      | CD80-low mTEC              | Valid, but too small for your preferred analysis |






| Feature            | GSE205410             | GSE222285               |
| ------------------ | --------------------- | ----------------------- |
| Sex                | **Female**            | **Male**                |
| Strain             | C57BL/6               | C57BL/6                 |
| Cell population    | MHCII-low mTEC        | MHCII-low mTEC          |
| WT replicates      | 8                     | 8                       |
| Aire-KO replicates | 4                     | 4                       |
| Internal labels    | `Aire_1`–`Aire_12`    | Also `Aire_1`–`Aire_12` |
| BioProject         | PRJNA844982           | PRJNA918978             |
| GEO samples        | GSM6211564–GSM6211575 | GSM6919454–GSM6919465   |
| Publication        | PMID 40244172         | PMID 36695731           |


GSE222285 and GSE205410 are different mice, but sama research group, thus methology is same.

GSE2585, which is only mTEC low dataset made with affy, was released in 2005 and it has about 10k genes. 
Newer datasets have 20k, both in affy and RNA seq.
GSE2585 genes overlap with other rna seq datasets. 
