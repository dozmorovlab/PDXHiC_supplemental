# Condition-specific and common loops and anchors

<!-- Desktop/PDXHiC/Maggie/30_Anchor_Loop_Barplots.Rmd -->
- `01_Anchor_Loop_Barplots.Rmd` - anchor and loop counts in each condition visualized with stacked barplots
  - Input: loop and anchor metrics file containing number of loops and anchors for each condition
  - Output: stacked barplot showing counts of loops and anchors for each condition separated by common and unique
  
<!-- HiC_files/results/Maggie/Mustache_results/preprocessing_any/Loops/Aggregated_Peak_Analysis/run_apa_mustache.qsub -->
- `02_APA_Mustache.qsub` - aggregated peak analysis done on Mustache loops over all resolutions
  - Input: `UCD52PR.hic`, `UCD52CR.hic` HiC matrices for each condition, loop files for each condition and loop type (Common/Unique) for all resolutions
  - Output: APA plots for each loop & HiC matrix pairwise combinations
  
<!-- Desktop/PDXHiC/Maggie/13_Loop_Width.sh -->
- `03_Loop_Width.sh` - loop size piecharts
  - Input: loop bedpe files for each condition and resolution
  - Output: piecharts for each condition and resolution showing loop size

<!--Desktop/PDXHiC/Maggie/18_CTCF_orientation.sh -->
- `04_CTCF_Orientation.sh` - CTCF orientation piecharts
  - Input: `fimo.bed` hg38 CTCF motif file, loop bedpe files for each condition and resolution
  - Output: CTCF orientation piecharts for each condition and resolution
  