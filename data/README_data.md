# Data Directory

This directory intentionally contains no raw data files.

## Automatic Download

The 16S rRNA sequencing data used in this analysis are sourced programmatically from the
Bioconductor experiment data package `HMP2Data` (Stansfield et al., 2024) and cached
locally by the ExperimentHub infrastructure upon first access.

To download the data, install the package and call the accessor:

```r
BiocManager::install("HMP2Data")
library(HMP2Data)
ps <- HMP2Data::IBD16S()  # downloads and caches automatically
```

## Source

- **Dataset:** Inflammatory Bowel Disease Multi-omics Database (IBDMDB)
- **Study:** integrative Human Microbiome Project (iHMP / HMP2)
- **Reference:** Lloyd-Price et al. (2019) Nature, 569, pp. 655-662.
  doi: 10.1038/s41586-019-1237-9
- **Package:** HMP2Data v1.25.0 (Artistic-2.0 licence)
  https://bioconductor.org/packages/HMP2Data

## Ethical Notes

Data are de-identified and publicly released. Original studies obtained
participant consent and institutional approval. No attempt is made to
re-identify any participant. See the HMP DACC terms at https://ibdmdb.org.
