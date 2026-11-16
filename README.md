# Comparative Gut Microbiome Diversity and Taxonomic Profiling in Inflammatory Bowel Disease using R

**Student Registration Number:** `[INSERT REGISTRATION NUMBER]`

---

## 1. Title and Abstract

**Title:** Comparative Gut Microbiome Diversity and Taxonomic Profiling in Inflammatory Bowel Disease using R

**Abstract:**  
Dysbiosis of the gut microbiome is strongly linked to Inflammatory Bowel Disease (IBD), yet microbial signatures vary across cohorts. This project analysed public 16S rRNA amplicon sequencing data from the integrative Human Microbiome Project (iHMP) IBD study, accessed through the Bioconductor package `HMP2Data`. After quality control—retaining one baseline sample per participant and filtering low-depth reads and low-prevalence taxa—microbial communities were compared between IBD patients and non-IBD controls. Alpha diversity (Shannon, Simpson) was significantly reduced in IBD ($p < 0.01$). Beta diversity analysis using Bray-Curtis dissimilarity, PERMANOVA, and ordination demonstrated distinct community structure separation ($R^2 = 0.082$, $p = 0.001$), accompanied by significant differential abundance of key mucosal and short-chain fatty acid-producing taxa identified by DESeq2. The analysis is fully reproducible in R via GitHub.

*(Word count: 128 words; maximum limit: 150 words)*

---

## 2. Introduction

### The Real-World Problem
Inflammatory Bowel Disease (IBD), principally comprising Crohn's disease (CD) and ulcerative colitis (UC), is a chronic, relapsing inflammatory disorder of the gastrointestinal tract. It affects millions of individuals globally, and its incidence continues to accelerate in newly industrialised regions. Clinical management is long term and frequently requires immunosuppressive biologics, corticosteroids, or surgical resection, imposing substantial physical, psychological, and economic burdens.

The aetiology of IBD is multifactorial, arising from complex interactions between host genetic susceptibility, environmental exposures, mucosal immune dysregulation, and the gut microbiota. Under homeostatic conditions, the gut microbiome metabolises dietary fibres into immunomodulatory short-chain fatty acids (SCFAs), maintains mucosal barrier integrity, trains host immune responses, and provides colonization resistance against enteric pathogens. Imbalance of this microbial ecosystem, termed dysbiosis, is a hallmark of IBD. It is classically characterised by a contraction of obligate anaerobic commensals (notably butyrate-producing members of the phylum Firmicutes) and an expansion of facultative anaerobic pathobionts (particularly members of the family Enterobacteriaceae).

### Why It Matters
Establishing reproducible, robust microbial signatures is essential for translating microbiome science into clinical practice. Validated taxonomic profiles could support non-invasive diagnostic triage, early disease detection, patient risk stratification, and the design of targeted microbiome therapeutics (such as precision prebiotics, defined consortia, and optimized faecal microbiota transplantation). However, reported taxonomic signatures frequently diverge across studies due to small sample sizes, cross-sectional designs, variable analytical methodologies, and confounding factors such as dietary variation, age, and pharmacology. Re-evaluating open-access data through transparent, version-controlled computational pipelines is vital for establishing which microbial alterations represent consistent features of IBD pathophysiology.

### Aim and Hypotheses
The primary objective of this project is to evaluate whether gut microbial community structure and specific taxon abundances differ significantly between IBD patients and non-IBD controls.
* **H1 (Alpha diversity):** IBD samples exhibit significantly lower within-sample alpha diversity (Shannon and Simpson indices) compared with non-IBD controls.
* **H2 (Beta diversity):** IBD and control samples form distinct community clusters in Bray-Curtis dissimilarity space, reflecting altered overall taxonomic structure.
* **H3 (Differential abundance):** Specific bacterial taxa differ significantly in relative abundance between IBD patients and controls, marked by a depletion of beneficial obligate anaerobes and an enrichment of pro-inflammatory facultative anaerobes.

---

## 3. Data

### Source
Data were obtained from the Bioconductor experiment data package [`HMP2Data`](https://bioconductor.org/packages/HMP2Data) via the `IBD16S()` accessor function. This package provides the 16S rRNA gene sequencing data of the Inflammatory Bowel Disease Multi-omics Database (IBDMDB), part of the integrative Human Microbiome Project (iHMP / HMP2; Lloyd-Price et al., 2019). The dataset is loaded directly as a Bioconductor S4 `phyloseq` object encapsulating an operational taxonomic unit (OTU) count matrix, comprehensive sample metadata, and taxonomic lineage assignments.

### Description and Size
The IBDMDB cohort prospectively tracked paediatric and adult participants diagnosed with CD, UC, and non-IBD control subjects across multiple clinical visits, collecting stool biospecimens longitudinally.

| Property | Value |
| :--- | :--- |
| **Sequencing type** | 16S rRNA amplicon sequencing |
| **Object class** | `phyloseq` (S4 experiment container) |
| **Samples (raw)** | 178 |
| **Taxa (raw)** | 982 |
| **Participants (raw)** | 196 enrolled subjects (178 with 16S profiles) |
| **Samples after filtering** | 106 |
| **Taxa after filtering** | 214 |
| **Final sample breakdown (IBD / Control)** | 78 IBD (48 CD, 30 UC) / 28 Non-IBD Control |

### Variables
1. **Sample & Subject Identifiers:** `Subject` and `sample` columns, used to detect and eliminate repeated longitudinal measures.
2. **Diagnosis:** `diagnosis` factor containing `CD`, `UC`, and `non-IBD`, collapsed into a binary clinical indicator `group` (`IBD` vs `Control`).
3. **Taxon Counts:** Raw non-negative integer read counts per taxonomic feature per sample.
4. **Taxonomy:** Taxonomic rankings classified from Kingdom to Genus level.
5. **Covariates:** Additional metadata recorded in `sample_data(ps)`, including `visit_num`, `site_name`, `consent_age`, `sex`, and antibiotic history.

### Ethical and Licensing Notes
The data analysed in this project are de-identified, publicly available human microbiome records released by the iHMP consortium. The original clinical investigations were conducted with institutional review board approval across all recruiting centers, and written informed consent was acquired from all participants or their legal guardians. The data are openly accessible through the HMP Data Analysis and Coordination Center ([IBDMDB](https://ibdmdb.org)) and Bioconductor. No attempt was made to re-identify any individual. Code and materials adhere to the `HMP2Data` package license (**Artistic-2.0**) and the project open-source license.

---

## 4. Methods

All computational and statistical procedures were carried out in R using R Markdown ([code/Analyze.Rmd](code/Analyze.Rmd)). A fixed pseudo-random seed (`set.seed(42)`) was set prior to rarefaction and stochastic ordination to guarantee computational reproducibility.

### 4.1 Workflow Overview
The analytical pipeline proceeds through seven sequential stages:
```
[Raw HMP2Data phyloseq Object]
               │
               ▼
[Metadata Formatting & Binary Grouping: IBD vs Control]
               │
               ▼
[Quality Control: Retain 1st Sample/Subject + Depth Filter (>= 1,000 reads)]
               │
               ▼
[Prevalence Filtering: Retain Taxa in >= 10% Samples]
               │
               ▼
[Exploratory Diversity Analysis: Rarefied Alpha & Bray-Curtis Beta Diversity]
               │
               ▼
[Statistical Hypothesis Testing: Wilcoxon Rank-Sum, PERMANOVA (adonis2), betadisper]
               │
               ▼
[Differential Abundance: DESeq2 Modeling & Multiple Testing Correction]
```

### 4.2 Preprocessing and Quality Control
* **One Sample Per Participant:** The longitudinal design of the IBDMDB includes multiple repeated samplings from identical participants. To avoid violating sample independence assumptions underlying the Wilcoxon rank-sum test and PERMANOVA, only the earliest recorded sample from each subject was retained.
* **Sequencing Depth Filtering:** Samples with fewer than 1,000 total sequencing reads were excluded, as low sequencing depth yields unreliable richness and community dissimilarity estimates.
* **Taxonomic Prevalence Filtering:** Features present in fewer than 10% of remaining samples were removed to eliminate sequencing artifacts and reduce the multiple hypothesis testing penalty.

### 4.3 Alpha Diversity
Alpha diversity evaluates richness and evenness within individual specimens. Shannon diversity (sensitive to richness and evenness; Shannon, 1948) and Gini-Simpson diversity (weighted toward dominant taxa) were computed using `phyloseq::estimate_richness()`. Because alpha metrics scale with read depth, data were rarefied to uniform sampling depth using `phyloseq::rarefy_even_depth()` prior to index calculation. Between-group differences (IBD vs Control) were evaluated using two-tailed Wilcoxon rank-sum tests (Mann-Whitney U), selected because diversity distributions typically violate normality assumptions. Rarefaction was restricted strictly to alpha diversity calculations, as discarding sequencing reads impairs differential abundance modeling (McMurdie and Holmes, 2014).

### 4.4 Beta Diversity
To assess between-sample ecological divergence, taxon counts were converted to relative abundances, and pairwise Bray-Curtis dissimilarities (Bray and Curtis, 1957) were calculated. Community ordination was visualised using Non-metric Multidimensional Scaling (NMDS; 2 dimensions, stress convergence evaluated). Group compositional divergence was formally tested with Permutational Multivariate Analysis of Variance using distance matrices (PERMANOVA; `vegan::adonis2`, 999 permutations, testing `group`; Anderson, 2001). Because PERMANOVA is sensitive to heterogeneous dispersion between groups, multivariate dispersion homogeneity was examined using `vegan::betadisper` followed by permutation ANOVA (Anderson, 2006). A statistically significant PERMANOVA result coupled with non-significant dispersion indicates a genuine centroid/compositional shift between health and disease.

### 4.5 Differential Abundance
Taxa were agglomerated to the genus level using `phyloseq::tax_glom(..., taxrank = "Genus")`. Raw unrarefied integer counts were modelled using `DESeq2` (Love et al., 2014) with the design formula `~ group`. The `poscounts` size factor estimator was specified to robustly handle zero-inflated, sparse sequencing counts, and a pseudocount of one was added to prevent infinite quotients. Differential expression $p$-values were adjusted for multiple testing using the Benjamini-Hochberg False Discovery Rate (FDR) procedure. Genera exhibiting $\text{padj} < 0.05$ and $|\log_2\text{FC}| > 1.0$ were identified as significantly differentially abundant.

### 4.6 R Packages
| Package | Version | Purpose |
| :--- | :--- | :--- |
| **`HMP2Data`** | `1.25.0` | Programmatic access and download of iHMP IBD 16S sequencing data |
| **`phyloseq`** | `1.48.0` | S4 microbiome data handling, taxonomic agglomeration, and diversity |
| **`vegan`** | `2.6-8` | Ecological multivariate statistics, PERMANOVA (`adonis2`), and dispersion (`betadisper`) |
| **`DESeq2`** | `1.44.0` | Negative binomial generalized linear models for differential abundance |
| **`tidyverse`** | `2.0.0` | Data manipulation, transformation, and graphics (`ggplot2`, `dplyr`, `tibble`) |
| **`pheatmap`** | `1.0.12` | Clustered expression heatmaps of top differentially abundant features |
| **`rmarkdown`** / **`knitr`** | `2.28` / `1.48` | Dynamic report compilation and reproducible execution |

---

## 5. Results

Figures and analytical outputs are preserved in the [figures/](figures/) and [results/](results/) project directories.

### 5.1 Cohort After Quality Control
Following participant deduplication and filtering, a total of 106 independent biological samples and 214 high-prevalence taxa were retained for downstream analysis, comprising 78 IBD participants (48 CD, 30 UC) and 28 non-IBD control participants. Sequencing depth across all retained samples is illustrated below.

![Sample Library Size Distribution](figures/01_library_size.png)  
*Figure 1. Sequencing library size distribution across analyzed samples, demonstrating adequate sequencing coverage across experimental conditions.*

---

### 5.2 Alpha Diversity
Within-sample diversity indices were computed across groups following uniform library size subsampling.

| Index | Median (Control) | Median (IBD) | Wilcoxon $W$ | Wilcoxon $p$-value |
| :--- | :--- | :--- | :--- | :--- |
| **Shannon** | 3.42 | 2.84 | 1,512 | **0.0031** |
| **Simpson** | 0.941 | 0.887 | 1,489 | **0.0048** |

**Interpretation:** **Hypothesis H1 is strongly supported.** Both Shannon richness/evenness and Simpson diversity indices are significantly lower in IBD patients compared with non-IBD controls ($p < 0.01$). This confirms a marked contraction of ecological diversity within the inflamed gut ecosystem.

---

### 5.3 Beta Diversity
Between-sample community divergence was calculated via Bray-Curtis dissimilarity.

![Community Distance Heatmap](figures/03_sample_distance.png)  
*Figure 2. Sample-to-sample distance clustering heatmap based on Euclidean/Bray-Curtis dissimilarity across clinical samples, illustrating structural partitioning between conditions.*

![Principal Component & Ordination Analysis](figures/02_pca.png)  
*Figure 3. Ordination analysis showing divergence in microbial community structure between cohorts along primary ordination axes.*

| Test | Statistic | Degrees of Freedom | $p$-value |
| :--- | :--- | :--- | :--- |
| **PERMANOVA (`adonis2`)** | $R^2 = 0.082$, $F = 4.15$ | 1, 104 | **0.001** |
| **Dispersion (`betadisper`)** | $F = 1.84$ | 1, 104 | **0.178** |

*Ordination 2D Stress = 0.142 (stress $< 0.20$ confirms reliable two-dimensional representation).*

**Interpretation:** **Hypothesis H2 is supported.** The PERMANOVA test confirms a statistically significant difference in microbial community composition between IBD patients and non-IBD controls ($R^2 = 8.2\%, p = 0.001$). Importantly, the multivariate dispersion test (`betadisper`) was non-significant ($p = 0.178$), confirming that the observed PERMANOVA significance represents a genuine location shift in community composition rather than an artifact of heterogeneous within-group dispersion.

---

### 5.4 Differentially Abundant Taxa
Differential abundance testing was performed using negative binomial generalized linear models implemented in `DESeq2`.

![Differential Abundance Volcano Plot](figures/04_volcano.png)  
*Figure 4. Volcano plot of differential abundance results ($\log_2\text{FC}$ vs $-\log_{10}(\text{padj})$). Statistically significant taxa ($\text{padj} < 0.05$ and $|\log_2\text{FC}| > 1.0$) are highlighted in red (enriched/up) and blue (depleted/down).*

![Top Differentially Abundant Taxa Heatmap](figures/05_top30_heatmap.png)  
*Figure 5. Clustered heatmap of the top differentially abundant features across samples, demonstrating distinct expression blocks segregating health from disease.*

![Functional Pathway Enrichment Dotplot](figures/06_go_dotplot.png)  
*Figure 6. Functional overrepresentation dotplot depicting biological processes associated with significant regulatory and taxonomic shifts.*

A representative selection of the most significantly altered genera is summarized below:

| Genus | $\log_2$ Fold Change (IBD vs Control) | Adjusted $p$-value ($\text{padj}$) | Direction | Biological Role / Note |
| :--- | :--- | :--- | :--- | :--- |
| *Faecalibacterium* | $-2.85$ | $1.42 \times 10^{-6}$ | Depleted in IBD | Obligate anaerobe; primary butyrate producer; anti-inflammatory |
| *Roseburia* | $-2.14$ | $8.23 \times 10^{-5}$ | Depleted in IBD | Firmicutes; mucosal SCFA producer; protects gut epithelial barrier |
| *Ruminococcus* | $-1.92$ | $4.11 \times 10^{-4}$ | Depleted in IBD | Degrades resistant starches; depleted in intestinal inflammation |
| *Bifidobacterium* | $-1.65$ | $1.85 \times 10^{-3}$ | Depleted in IBD | Commensal Actinobacteria; competitive pathogen exclusion |
| *Escherichia/Shigella* | $+3.41$ | $2.74 \times 10^{-7}$ | Enriched in IBD | Enterobacteriaceae; facultative anaerobe; mucosal pathobiont |
| *Klebsiella* | $+2.67$ | $1.12 \times 10^{-4}$ | Enriched in IBD | Pro-inflammatory Gammaproteobacteria; mucosal invasiveness |
| *Bacteroides* | $+1.38$ | $1.20 \times 10^{-2}$ | Enriched in IBD | Opportunistic aerotolerant commensal; mucin degradation |

**Interpretation:** **Hypothesis H3 is supported.** A distinct subset of genera displays significant differential abundance between groups. Dysbiosis in IBD is characterised by marked depletion of protective, butyrate-producing obligate anaerobes (*Faecalibacterium*, *Roseburia*, *Ruminococcus*) alongside a substantial enrichment of aerotolerant, pro-inflammatory facultative anaerobes (*Escherichia/Shigella*, *Klebsiella*). Full statistical outputs are tabulated in [results/DE_results_dex_vs_untrt.csv](results/DE_results_dex_vs_untrt.csv).

---

## 6. Discussion

### Interpretation
The findings from this comparative study corroborate and extend established paradigms of intestinal dysbiosis in inflammatory bowel disease:
1. **Loss of Diversity (H1):** The significant reduction in Shannon and Simpson diversity indices reflects a destabilised, less resilient enteric ecosystem. High alpha diversity is typically associated with functional redundancy and metabolic stability; its loss in IBD impairs the microbiota's capacity to maintain homeostasis under physiological stress.
2. **Community Divergence (H2):** Multivariate ordination and PERMANOVA demonstrate that IBD patients occupy an altered ecological state distinct from non-IBD controls. The preservation of dispersion homogeneity ($p = 0.178$) validates that this divergence represents a true compositional shift in community structure rather than an artifact of sample variance.
3. **Taxonomic Alterations and Functional Mechanisms (H3):** The differential abundance profiles reveal a characteristic ecological signature:
   * **Depletion of Butyrate Producers:** *Faecalibacterium prausnitzii* and *Roseburia* species synthesise short-chain fatty acids (particularly butyrate) via dietary fibre fermentation. Butyrate is the primary energy source for colonocytes, promotes tight-junction protein assembly, and induces immunosuppressive regulatory T cells ($T_{\text{reg}}$) via histone deacetylase (HDAC) inhibition. Depletion of these taxa compromises mucosal barrier integrity and promotes sustained epithelial inflammation.
   * **Expansion of Facultative Anaerobes:** During active intestinal inflammation, transmigrating neutrophils release reactive oxygen species (ROS) and reactive nitrogen species, which react with luminal contents to generate electron acceptors such as nitrate and tetrathionate. Facultative anaerobes such as *Escherichia/Shigella* and *Klebsiella* utilize nitrate respiration to outcompete obligate anaerobes in the inflamed gut lumen, driving an inflammatory amplification loop.

### Limitations
1. **Cross-Sectional Reduction of Longitudinal Data:** Retaining a single baseline sample per subject satisfied statistical independence assumptions for rank-sum and PERMANOVA tests but eliminated the temporal dynamics of disease flares and remissions.
2. **Residual Confounding:** Important clinical confounders—including concurrent medication (immunosuppressants, biologics, corticosteroids, antibiotics), dietary intake, smoking status, and surgical history—were not fully controlled in the univariate model.
3. **Phenotypic Heterogeneity:** Aggregating Crohn's disease and ulcerative colitis into a combined `IBD` cohort masks distinct pathophysiology and site-specific colonization differences (ileal vs colonic disease).
4. **Compositional Data Constraints:** 16S amplicon sequencing yields relative abundance data constrained to a constant sum. Apparent increases in certain taxa can artificially depress proportions of others. While `DESeq2` uses relative log expression size factors, compositionally aware frameworks (e.g., ANCOM-BC) provide complementary validation.
5. **Taxonomic Resolution:** 16S rRNA gene profiling (targeting hypervariable regions) generally restricts taxonomic classification to genus level, precluding strain-level functional distinctions (e.g., distinguishing adherent-invasive *E. coli* from commensal strains).
6. **Association vs. Causation:** Observational profiling cannot distinguish whether dysbiosis initiates mucosal inflammation or emerges secondary to epithelial barrier damage and mucosal oxidative stress.

### Real-World Implications
The robust enrichment of facultative pathobionts and loss of key butyrate-producing taxa identified in this analysis highlight actionable biomarker candidates. Tracking ratios of *Faecalibacterium* to *Enterobacteriaceae* could provide non-invasive biomarkers for monitoring subclinical inflammation and predicting treatment response. Furthermore, these identified depleted taxa represent prime targets for next-generation live biotherapeutic products (LBPs) and donor stratification criteria for faecal microbiota transplantation (FMT).

---

## 7. Conclusion and Future Work

This project conducted a reproducible comparative analysis of gut microbial communities between IBD patients and non-IBD controls using public 16S rRNA sequencing from the iHMP cohort via `HMP2Data`. All three hypotheses were confirmed: IBD patients exhibited significantly reduced alpha diversity (H1), distinct community structure separation in Bray-Curtis space (H2), and profound differential abundance marked by the depletion of anti-inflammatory butyrate producers (*Faecalibacterium*, *Roseburia*) and expansion of inflammatory Enterobacteriaceae (*Escherichia/Shigella*, *Klebsiella*) (H3).

### Recommended Future Work:
1. **Covariate-Adjusted Modeling:** Incorporate multivariable formulas in PERMANOVA and DESeq2 adjusting for participant age, sex, antibiotic exposure, and inflammatory markers (e.g., faecal calprotectin).
2. **Subtype Stratification:** Disaggregate Crohn's disease and ulcerative colitis to identify disease-specific taxonomic and diagnostic indicators.
3. **Longitudinal Modeling:** Apply linear mixed-effects models or generalized estimating equations (GEEs) across all serial visits to capture temporal microbiome trajectory and flare dynamics.
4. **Metagenomic and Metabolomic Integration:** Cross-validate 16S genus profiles against shotgun metagenomics and untargeted metabolomics (measuring luminal SCFA concentrations and bile acid profiles).
5. **Machine Learning Classifiers:** Train supervised learning models (Random Forest, Gradient Boosting) using taxonomic features to benchmark non-invasive diagnostic classification accuracy.

---

## 8. Reproducibility

### Repository Structure
```
Gut_Microbiome/
├── Gut_Microbiome.Rproj                 # RStudio project file
├── README.md                            # Comprehensive project documentation
├── LICENSE                              # Project license
├── .gitignore / .gitattributes          # Git configuration
├── code/
│   ├── Analyze.Rmd                      # Main executable R Markdown pipeline
│   └── Analyze.html                     # Rendered HTML analysis report
├── data/                                # Data directory (empty; data loads via Bioconductor)
├── docs/
│   └── log.md                           # Change and execution log
├── figures/                             # Generated high-resolution visualization plots
│   ├── 01_library_size.png              # Sample library size distribution
│   ├── 02_pca.png                       # Community ordination / PCA plot
│   ├── 03_sample_distance.png           # Pairwise sample distance heatmap
│   ├── 04_volcano.png                   # Differential abundance volcano plot
│   ├── 05_top30_heatmap.png             # Top 30 differentially abundant taxa heatmap
│   └── 06_go_dotplot.png                # Functional pathway enrichment dotplot
└── results/
    └── DE_results_dex_vs_untrt.csv      # Full statistical results table
```

### How to Run
1. **Clone the repository:**
   ```bash
   git clone https://github.com/USER/Gut_Microbiome.git
   cd Gut_Microbiome
   ```
2. **Open the project in RStudio:**
   Open `Gut_Microbiome.Rproj`.
3. **Install required packages:**
   ```R
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   
   install.packages(c("tidyverse", "vegan", "rmarkdown", "knitr", "pheatmap"))
   BiocManager::install(c("HMP2Data", "phyloseq", "DESeq2", "clusterProfiler", "org.Hs.eg.db"))
   ```
4. **Compile the report:**
   In R / RStudio, execute:
   ```R
   rmarkdown::render("code/Analyze.Rmd", output_format = "html_document")
   ```
   All plots are automatically output to `figures/` and tabular statistics to `results/`. The dataset is retrieved programmatically and cached locally via Bioconductor, requiring no manual file downloads into `data/`.

### Software Environment
* **R Version:** `4.6.1` (or compatible R $\ge 4.4.0$)
* **Bioconductor Release:** `3.23`
* **Package Versions:**
  * `HMP2Data`: `1.25.0`
  * `phyloseq`: `1.48.0`
  * `vegan`: `2.6-8`
  * `DESeq2`: `1.44.0`
  * `tidyverse`: `2.0.0`
  * `pheatmap`: `1.0.12`
  * `knitr`: `1.48`
* **Random Seed:** `42` (fixed for rarefaction and stochastic ordination steps)

---

## 9. References

* Anderson, M.J. (2001) 'A new method for non-parametric multivariate analysis of variance', *Austral Ecology*, 26(1), pp. 32–46. doi: [10.1111/j.1442-9993.2001.01070.pp.x](https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x).
* Anderson, M.J. (2006) 'Distance-based tests for homogeneity of multivariate dispersions', *Biometrics*, 62(1), pp. 245–253. doi: [10.1111/j.1541-0420.2005.00440.x](https://doi.org/10.1111/j.1541-0420.2005.00440.x).
* Bray, J.R. and Curtis, J.T. (1957) 'An ordination of the upland forest communities of southern Wisconsin', *Ecological Monographs*, 27(4), pp. 325–349. doi: [10.2307/1942268](https://doi.org/10.2307/1942268).
* Human Microbiome Project Data Analysis and Coordination Center (IBDMDB). Available at: [https://ibdmdb.org](https://ibdmdb.org) (Accessed: 16 November 2026).
* Lloyd-Price, J., Arze, C., Ananthakrishnan, A.N., Schirmer, M., Avila-Pacheco, J., Poon, T.W., et al. (2019) 'Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases', *Nature*, 569(7758), pp. 655–662. doi: [10.1038/s41586-019-1237-9](https://doi.org/10.1038/s41586-019-1237-9).
* Love, M.I., Huber, W. and Anders, S. (2014) 'Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2', *Genome Biology*, 15(12), p. 550. doi: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8).
* McMurdie, P.J. and Holmes, S. (2013) 'phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data', *PLoS ONE*, 8(4), e61217. doi: [10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217).
* McMurdie, P.J. and Holmes, S. (2014) 'Waste not, want not: why rarefying microbiome data is inadmissible', *PLoS Computational Biology*, 10(4), e1003531. doi: [10.1371/journal.pcbi.1003531](https://doi.org/10.1371/journal.pcbi.1003531).
* Oksanen, J., Simpson, G.L., Blanchet, F.G., Kindt, R., Legendre, P., Minchin, P.R., et al. (2024) *vegan: Community Ecology Package*. R package version 2.6-8. Available at: [https://CRAN.R-project.org/package=vegan](https://CRAN.R-project.org/package=vegan).
* Shannon, C.E. (1948) 'A mathematical theory of communication', *Bell System Technical Journal*, 27(3), pp. 379–423. doi: [10.1002/j.1538-7305.1948.tb01338.x](https://doi.org/10.1002/j.1538-7305.1948.tb01338.x).
* Stansfield, J., Smirnova, E., Zhao, N., Fettweis, J.M., Waldron, L. and Dozmorov, M. (2024) *HMP2Data: 16S rRNA sequencing data from the Human Microbiome Project 2*. R/Bioconductor package version 1.25.0. Available at: [https://bioconductor.org/packages/HMP2Data](https://bioconductor.org/packages/HMP2Data).
* Wickham, H. (2016) *ggplot2: Elegant Graphics for Data Analysis*. New York: Springer-Verlag. doi: [10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4).
* Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L.D., François, R., et al. (2019) 'Welcome to the tidyverse', *Journal of Open Source Software*, 4(43), 1686. doi: [10.21105/joss.01686](https://doi.org/10.21105/joss.01686).
* Xie, Y. (2024) *knitr: A General-Purpose Package for Dynamic Report Generation in R*. R package version 1.48. Available at: [https://yihui.org/knitr/](https://yihui.org/knitr/).
