# Differentially expressed transcripts and their functional significance.

<!-- /Users/mdozmorov/Documents/Work/VCU_work/ChuckHarrell/ChuckHarrell_BrainMetastasis_12-2016/2021.10.28_NewBulkRNASeq_Preprocessing/10.28.21_TPMCalc_DataPreproc.R -->
- `01_TPMCalc_DataPreproc.R` - Preprocessing of Salmon gene counts
    - Input: `*.quant.sf` sample-specific Salmon output files
    - Output: `Count_human_AllGenes_2021_BulkRNASeq_10.28.21.tsv`, `TPM_human_AllGenes_2021_BulkRNASeq_10.28.21.tsv`

<!-- /Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/RNA-seq/2021.11.19_UCD52_DEG/04.06.22_UCD52_DEG_Analysis_edgeR.R -->
- `02_DEG_Analysis_edgeR.R` - Differentially expressed genes detection
    - Input: `Count_human_AllGenes_2021_BulkRNASeq_10.28.21.tsv`, `TPM_human_AllGenes_2021_BulkRNASeq_10.28.21.tsv`
    - Output: `Data_Human_DESeq2_UCD52_DEGs_Untreated.vs.CarboplatinResistant_04.08.22.txt`

<!-- /Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/RNA-seq/2021.11.19_UCD52_DEG/04.06.22_GSEA_edgeR.Rmd -->
- `03_GSEA_edgeR.Rmd` - GSEA analysis of differentially expressed genes
    - Input: `Data_Human_DESeq2_UCD52_DEGs_Untreated.vs.CarboplatinResistant_04.06.22.txt`
    - Output: `GSEA_edgeR_UCD52PR_CR_0.1_1.xlsx` - GSEA results; `04.06.22_DEGs_edgeR_UCD52PR_CR_annotated.xlsx` - annotated DEGs

<!-- /Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/RNA-seq/2021.11.19_UCD52_DEG/Figure_Volcano_protein_coding.Rmd -->
- `04a_Figure_Volcano_protein_coding.Rmd`, `04b_Figure_Volcano_noncoding.Rmd` - volcano plot figures
    - Input: `04.06.22_DEGs_edgeR_UCD52PR_CR_annotated.xlsx`
    - Output: `Figure_Volcano_protein_coding.svg`, `Figure_Volcano_noncoding.svg`

<!-- /Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/RNA-seq/2021.11.19_UCD52_DEG/Pathview.Rmd -->
- `05_Pathview.Rmd` - KEGG pathway plots overlaying gene expression
    - Input: `GSEA_edgeR_UCD52PR_CR_0.1_1.xlsx`, `04.06.22_DEGs_edgeR_UCD52PR_CR_annotated.xlsx`
    - Output: `pathways_selected.pdf`
    
<!-- Desktop/PDXHiC/Maggie/34_LncSEA_Visualizations.Rmd -->
- `06_LncSEA_Visualizations.Rmd` - upregulated lncRNAs in CR condition.
  - Input: `Supplementary_Table_1.xlsx` sheets with output from running LncSEA on lncRNAs upregulated in CR condition.
  - Output: barplot showing top lncRNA sets for different LncSEA classes. 
  
<!-- Desktop/PDXHiC/Maggie/40_GSEA_figures_RNAseq_dcHiC.Rmd -->
- `07_GSEA_Enrichment_Figure.Rmd` - barplots of most significant pathways from running GSEA analysis.
  - Input: KEGG & MSigDb GSEA enrichment produced from `03_GSEA_edgeR.Rmd`. 
  - Output: barplots of most significant KEGG & MSigDb pathways enriched in differentially expressed genes. 
  
  

