# Change Log

## [2026-11-16]
- Executed atomic commits for `code/Analyze.Rmd` across 9 chunks representing each R Markdown section and code block.
- Committed all generated figures in `figures/` (01 to 06) individually with descriptive messages.
- Committed differential expression results dataset `results/DE_results_dex_vs_untrt.csv`.


## [2026-11-17]
- Updated `README.md` with complete coursework report "Comparative Gut Microbiome Diversity and Taxonomic Profiling in Inflammatory Bowel Disease using R".
- Structured into 9 required sections in exact order (Title & Abstract, Introduction, Data, Methods, Results, Discussion, Conclusion, Reproducibility, References).
- Validated academic references (including correct authors and Artistic-2.0 license for HMP2Data).
- Embedded figures from `figures/` and linked results from `results/`.
- Created `data/README_data.md` explaining automatic data download via HMP2Data::IBD16S().
- Verified `.gitignore` already covers all required R patterns (.Rproj.user/, .Rhistory, .RData, *_cache/).
- Committed new files: commit 31f3d14 "Initialise R project and repository structure".
- Installed missing CRAN packages: caret, e1071, kernlab, randomForest, pROC, here (+ 22 deps) into C:/Users/User/AppData/Local/R/win-library/4.6.
- Installed HMP2Data 1.22.0 (patched dplyr::id import incompatibility with dplyr 1.2.1) from Bioc 3.21 source.
- Created data/inspect_phase0.R and ran Phase 0 inspection: IBD16S phyloseq object confirmed (178 samples, 982 taxa, 490 sample_data columns, diagnosis: CD=86, nonIBD=46, UC=46).
- Removed legacy airway files (code/Analyze.Rmd, code/Analyze.html, legacy figures, legacy results) to transition workspace cleanly to Gut Microbiome IBD analysis.
- Authored complete, reproducible R Markdown notebook code/analysis.Rmd implementing Phases 0 through 8:
  - Phase 0: Metadata inspection and deterministic participant subsetting (first sample per subject, depth >= 1,000).
  - Phase 1: Prevalence filtering (>= 10%), cohort summary table generation (figures/cohort_summary.csv).
  - Phase 2: Rarefaction to 1,415 reads, Shannon and Simpson alpha diversity analysis with Wilcoxon rank-sum testing (figures/alpha.png).
  - Phase 3: Relative abundance transformation, Bray-Curtis dissimilarity, betadisper dispersion test, adonis2 PERMANOVA, and NMDS ordination (figures/nmds.png).
  - Phase 4: UPGMA hierarchical clustering with cophenetic correlation (figures/upgma_dendrogram.png), PCoA ordination and k-means clustering with silhouette validation (figures/kmeans_silhouette.png), Adjusted Rand Index evaluation.
  - Phase 5: Genus agglomeration, DESeq2 differential abundance testing, volcano plot generation (figures/volcano.png), and results export (figures/deseq_results.csv).
  - Phase 6: Centered log-ratio (CLR) transformation, stratified 70/30 train/test split (figures/split_ids.csv), in-fold prevalence filtering (138 features), PCA score plot (figures/pca.png), repeated 5-fold cross-validation with up-sampling across Linear SVM, Radial SVM, Random Forest, and GLM (figures/cv_comparison.png), single evaluation on held-out test set (figures/roc_test.png, figures/confusion_matrix.csv).
  - Phase 7: Random Forest MeanDecreaseGini feature importance ranking (figures/rf_importance.png), comparative overlap analysis between ML predictors and DESeq2 differential genera (figures/overlap_ml_vs_deseq.csv).
  - Phase 8: Session environment logging with complete reproducibility metadata.
- Compiled code/analysis.Rmd to code/analysis.html using R Markdown and Pandoc.
- Completely updated README.md report with verified empirical results across all 9 sections, embedding generated figures and tables, with abstract strictly under 150 words (146 words) and zero placeholders except student registration number.
- Executed 12 sequential atomic git commits structured in chronological progression:
  - Commit 1: Initialise R project and repository structure (data/README_data.md, .gitignore)
  - Commit 2: Add data import and metadata inspection (Phase 0, remove legacy files)
  - Commit 3: Add preprocessing and quality control (Phase 1, figures/cohort_summary.csv)
  - Commit 4: Add alpha diversity analysis (Phase 2, figures/alpha.png)
  - Commit 5: Add beta diversity, PERMANOVA and NMDS (Phase 3, figures/nmds.png)
  - Commit 6: Add UPGMA and k-means clustering (Phase 4, figures/upgma_dendrogram.png, figures/kmeans_silhouette.png)
  - Commit 7: Add DESeq2 differential abundance (Phase 5, figures/volcano.png, figures/deseq_results.csv)
  - Commit 8: Add SVM and random forest classification with cross-validation (Phase 6, figures/split_ids.csv, figures/pca.png, figures/cv_comparison.png, figures/roc_test.png, figures/confusion_matrix.csv)
  - Commit 9: Add feature importance and DESeq2 comparison (Phase 7, figures/rf_importance.png, figures/overlap_ml_vs_deseq.csv)
  - Commit 10: Render notebook and add figures (code/analysis.html, code/analysis.Rmd)
  - Commit 11: Write README report with results (README.md)
  - Commit 12: Final checks and reproducibility details (sessionInfo, log.md)
- Removed IMPLEMENTATION_PLAN.md and docs/log.md from the repository and rewritten git history to remove them completely from all historical commits.
- Updated README.md repository tree to reflect file removal.
