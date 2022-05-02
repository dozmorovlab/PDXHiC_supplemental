# Whole Genome Sequencing coverage differences between the CR and PR conditions

<!-- /Users/mdozmorov/Documents/Work/GitHub/Katarzyna/PDXHiC/WGS/07_DNAcopy.Rmd -->
- `07_DNAcopy.Rmd` - log2ratio segmentation
    - Input: bigWig files `UCD52CR_vs_PR_hic_10kb_log2ratio.bw` from either wgs or hic coverage folders
    - Output: diagnostic plots and a BED file with the detected DUP and DEL segments

<!-- /Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/WGS/submit08_htseq.sh -->
- `submit08_htseq.sh` - Counting reads in genes
    - Input: BAM files
    - Output: `UCD52PR.txt`, `UCD52CR.txt` gene counts

<!-- /Users/mdozmorov/Documents/Work/GitHub/Katarzyna/PDXHiC/WGS/09_gene_coverage_diff_04.06.22.pdf -->
- `09_gene_coverage_diff_04.06.22.Rmd` - correlating gene expression and coverage
    - Input: `UCD52PR.txt`, `UCD52CR.txt` gene counts, `04.06.22_DEGs_edgeR_UCD52PR_CR_annotated.xlsx` from RNA-seq analysis
    - Output: `04.06.22_UCD52PRCR_coverage_differences_all.svg`, `04.06.22_UCD52PRCR_coverage_differences_ABC.svg` - coverage-expression correlation plots; `04.06.22_Tables_WGS_Gene_Coverage.xlsx`, `04.06.22_UCD52PRCR_coverage_differences.xlsx` - gene coverage differences
