# Supplementary material for the PDX Hi-C paper

Work in progress.

- [00_preprocessing.Rmd](00_preprocessing.Rmd) - Separates loop data in any text format into PR/CR-specific and common loops and anchors. If `merge_adjacent <- TRUE`, nearby (flanking) regions are considered common. The code chunks in the "Load data" section create two GInteraction objects non-overlapping with excludable regions (centromeres, telomeres, etc.), which are processed using the same downstream code. Only one code chunk corresponding to a specific data type should be set to `eval = TRUE`. Data from Mustache, HiCcompare, NeoLoopFinder,  GENOVA, SpectralTAD, and hicFindTADs can be loaded. Data loading chunks for other data formats can be written. 
    - Input: Loop data. Set the corresponding code chunk to `eval=TRUE`
    - Output: Separated loops and anchors in the `preprocessing` subfolder. PR/CR `all`, `common`, and `unique` loops/anchors are outputted as BEDPE/BED files. Color-specifying columns are included to emphasize PR/CR/common loops/anchors (green/red/blue colors). Loops/anchors are not collapsed (reduced); therefore, `PR_common` and `CR_common` loops/anchors may differ. `log.txt` - counts and width of loops and anchors.
