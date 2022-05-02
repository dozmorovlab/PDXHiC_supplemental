# Replicability and compartmentalization changes in drug resistance

<!-- /Users/mdozmorov/Documents/mount/juicer/Mikhail/PDXproject_analysis/04.Arima_replicates_HiCrep/run_HiCrep.qsub -->
- `01_HiCrep.sh` - HiCrep analysis of replicates
    - Input: Hi-C replicates in .mcool format, converted using https://github.com/4dn-dcic/hic2cool
    - Output: pairwise chromosome-specific HiCrep measures, files like `sample1_sample2.txt`

<!-- /Users/mdozmorov/Documents/Work/GitHub/Katarzyna/PDXHiC/Mikhail/04_replicates_HiCrep.Rmd -->
- `02_replicates_HiCrep.Rmd` - HiCrep using individual replicates. 
    - Input: pairwise chromosome-specific HiCrep measures, files like `sample1_sample2.txt`. Correlation summarized as mean/median.
    - Output: MDS plots

<!-- /Users/mdozmorov/Documents/Work/GitHub/Katarzyna/PDXHiC/Mikhail/02_replicates_straw.Rmd -->
- `03_replicates_correlation.Rmd` - correlation using individual replicates. 
    - Input: straw-extracted sparse chromosome-specific matrices at 1Mb and 100kb resolutions. Correlation summarized as mean/median.
    - Output: correlation plots

<!-- /Users/mdozmorov/Documents/Work/GitHub/Katarzyna/PDXHiC/Mikhail/GENOVA_distance.Rmd -->
- `04_distance.Rmd` - distance-dependent and differential analysis
    - Input: Sample-specific .hic files
    - Output: Genome-wide and chromosome-specific difference decay plots

<!-- /Users/mdozmorov/Documents/Work/GitHub/Katarzyna/PDXHiC/Mikhail/GENOVA_ABsaddle.Rmd -->
- `05_ABsaddle.Rmd` - AB compartment and saddle analysis
    - Input: Sample-specific 500kb and 50kb mcool files
    - Output: pdf with AB compartment and saddle plots