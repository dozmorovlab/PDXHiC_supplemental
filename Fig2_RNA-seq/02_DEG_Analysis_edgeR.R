set.seed(1)
library(dplyr)
library(readr)
library(writexl)
library(edgeR)
library(annotables)

# Prepare annotations
gene_annotations <- grch38[!(grepl("_", grch38$chr) | grepl("GL", grch38$chr)), c("ensgene", "symbol", "description", "biotype")]
gene_annotations <- gene_annotations[!is.na(gene_annotations$symbol), ]
gene_annotations <- gene_annotations[!duplicated(gene_annotations), ]

# Settings
dir_data       <- "/Users/mdozmorov/Documents/Work/VCU_work/ChuckHarrell/ChuckHarrell_BrainMetastasis_12-2016/2021.10.28_NewBulkRNASeq_Preprocessing" # Working directory
fileNameIn1    <- file.path(dir_data, "processed_data/Count_human_AllGenes_2021_BulkRNASeq_10.28.21.tsv")
fileNameIn2    <- file.path(dir_data, "processed_data/TPM_human_AllGenes_2021_BulkRNASeq_10.28.21.tsv")
dir_results    <- "/Users/mdozmorov/Documents/Data/GoogleDrive/HiC_files/results/RNA-seq/2021.11.19_UCD52_DEG/data"
# fileNameOut1   <- file.path(dir_data, "Data_Human_DESeq2_UCD52_DEGs_Untreated.vs.CarboplatinResistant_04.06.22.txt") # File to save all results of differential expression analysis
fileNameOut1   <- file.path(dir_data, "Data_Human_DESeq2_UCD52_DEGs_Untreated.vs.CarboplatinResistant_04.08.22.txt") # File to save all results of differential expression analysis

# Read the full count table
counts <- read_tsv(fileNameIn1)
# Select UCD52-specific samples
counts <- counts %>% dplyr::select("gene_id", "symbol", starts_with("UCD52"))
# Make annotation data frame
sample_annotation <- data.frame(Sample = colnames(counts)[grepl("UCD52", colnames(counts))],
                                Group  = ifelse(grepl("CR", colnames(counts)[grepl("UCD52", colnames(counts))]), "CR", "PR"))
# Read the TPM counts
TPM <- read_tsv(fileNameIn2)
# Select UCD52-specific samples
TPM <- TPM %>% dplyr::select("gene_id", starts_with("UCD52"))
# Add average per condition
TPM <- TPM %>% dplyr::mutate(log2Average.PR = rowMeans(TPM %>% dplyr::select(starts_with("UCD52."))), .before = starts_with("UCD52")) %>%
               dplyr::mutate(log2Average.CR = rowMeans(TPM %>% dplyr::select(starts_with("UCD52CR"))), .before = starts_with("UCD52"))
# Convert to log2
TPM <- cbind(TPM %>% dplyr::select("gene_id"), log2((TPM %>% dplyr::select(-c("gene_id"))) + 1))

## Analysis 1
# "CR" vs. "PR" comparison, so positive logFC correspond to genes UP in "CR"
# Define reference and treatment groups
group_ref <- "PR"
group_trt <- "CR"
# Define analysis ID
analysis_ID <- paste0(group_trt, "_vs_", group_ref)
# Subset the data, adjust manually
index <- sample_annotation$Group == group_trt | sample_annotation$Group == group_ref
# Or, use all the data
# index <- sample_annotation$Sample == sample_annotation$Sample

sample_annotation_subset <- sample_annotation[index, ]
counts_subset <- counts[, 3:ncol(counts)]
counts_subset <- data.frame(Geneid = counts$gene_id, counts_subset[, index])
all.equal(sample_annotation_subset$Sample, colnames(counts_subset[, 2:ncol(counts_subset)]))

# Adjust manually
Group <- factor(sample_annotation_subset$Group)
Group <- relevel(Group, ref = group_ref)
Group
design <- model.matrix(~Group, data = Group)

# Create edgeR object
edgeR.dgelist = DGEList(counts = counts_subset[, 2:ncol(counts_subset)], genes = counts_subset$Geneid, group = Group)
# Filtering
keep <- rowSums(cpm(edgeR.dgelist)>1) >= 2
edgeR.dgelist <- edgeR.dgelist[keep, , keep.lib.sizes=FALSE]
# Normalization
edgeR.dgelist = calcNormFactors((edgeR.dgelist), method = "TMM")
edgeR.dgelist = estimateDisp(edgeR.dgelist, design)

# One-way ANOVA analysis
fit <- glmFit(edgeR.dgelist, design)
# Individual comparisons
lrt <- glmLRT(fit, coef = 2)

res.full <- topTags(object = lrt, n = Inf)
res.full_to_save <- left_join(res.full$table, gene_annotations, by = c("genes" = "ensgene"))
# Append TPMs
res.full_to_save <- left_join(res.full_to_save, TPM, by = c("genes" = "gene_id"))
# Save the data
write_tsv(res.full_to_save, fileNameOut1)
