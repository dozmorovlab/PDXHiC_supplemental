## Amy Olex
## 10.28.21
## Salmon and tximport pipeline for new BulkRNASeq data with 35 samples.  Code was directly copied from the original bulk experiment with Data Freeze on 8/18/17
## Much of this code was copied from http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
##
## This script imports the human and mouse data for all samples from the Salmon files. 
## It extracts the read counts and recalculates the TPM values for human and mouse seperatly.
## Gene names and the prefered sample names are appended to the TPM and raw count data frames.
## Two data files are then written out, one with mouse expression and one with human expression for each sample.
## 
## This script does nothing else except format the raw expression values.
## Old update log:
## UPDATE 11/15/17: To include correction for percent human and saving of the extra TPM file.
## UPDATE 11/28/17: To use percentage values that only include human and mouse mapped reads.  Unmapped reads were not included in the total.
## UPDATE 04/11/18: To also output the unfiltered gene list that is normalized per human percent.  All of these lines are commented out.
## UPDATE 08/22/18: To save raw, unprocessed TPM values for all genes for GEO upload.
## UPDATE 10/28/21: Modified data import to pull in file with both species and seperate in R.

# Location of .quant.sf files with sample-specific gene counts, from Salmon
setwd("/Volumes/GoogleDrive/My Drive/Active Collaborations/CHarrell/2021_BulkRNASeq/")

library(tximport)
library(readr)
library(AnnotationHub)
library(ensembldb)
library(RNASeqBits) # https://github.com/AmyOlex/RNASeqBits
library(NMF)
library(limma)

## Get all the names and the path of the files.
files <- file.path("quant_files", list.files("quant_files"))
## have to name files if I want my data columns to be named.
names(files) <- sub("_001.trimmed.quant.sf", "", sub("quant_files/", "", files))


## check to ensure they exist
stopifnot(all(file.exists(files)))

## Build a tx2gene object that will summarize all transcript expressions to the gene level
ah <- AnnotationHub()

annot_h <- query(ah, patter=c("Homo","EnsDb", "87"))
EnsDb_h <- annot_h[[1]]
df_human <- transcripts(EnsDb_h, return.type = "DataFrame")
tx2gene_h <- df_human[,c("tx_id","gene_id")]

annot_m <- query(ah, patter=c("Mus","EnsDb", "88"))
EnsDb_m <- annot_m[[1]]
df_mouse <- transcripts(EnsDb_m, return.type = "DataFrame")
tx2gene_m <- df_mouse[,c("tx_id","gene_id")]

## Import Salmon files and summarize to the gene level
## Note it will only recognize the transcript IDs that are in the txi list.  Thus, importing the same files that contain both
## human and mouse we can easily seperate based on the txi lists of transcripts.
txi_h <- tximport(files, type = "salmon", tx2gene = tx2gene_h, ignoreTxVersion = TRUE, importer = function(x) read_tsv(x, col_types="cnnnn"))
txi_m <- tximport(files, type = "salmon", tx2gene = tx2gene_m, ignoreTxVersion = TRUE, importer = function(x) read_tsv(x, col_types="cnnnn"))

#### Re-calculate the TPM values from the given counts and lengths
source("https://raw.githubusercontent.com/AmyOlex/WCCTR_RNASeq_Pipeline/d5e637e21b962a2e7354873f1f04edf39c38be04/WorkstationScripts/calc.tpm.fromSalmon.R")
# source("/Users/alolex/Desktop/CCTR_Git_Repos/WCCTR_RNASeq_Pipeline/WorkstationScripts/calc.tpm.fromSalmon.R")
TPM_h <- calc.tpm.fromSalmon(txi_h, colnames(txi_h$counts))
TPM_m <- calc.tpm.fromSalmon(txi_m, colnames(txi_m$counts))

## Get raw count data
counts_h <- round(as.data.frame(txi_h$counts))
#names(counts_h) <- names(TPM_h)

counts_m <- round(as.data.frame(txi_m$counts))
#names(counts_m) <- names(TPM_m)

################### Save raw TPM values for GEO ####################

###### 
##For this data I don't need the prefered samples name because it is already labeled as such.
## Load in sample meta data
#samples <- read.delim(file="~/Desktop/CCTR/Data/ChuckHarrell_BrainMetastasis_12-2016/8-18-17_Salmon_RawCountDataFreeze/Sample_Metadata_Freeze_8-18-17_updated_11-27-17.txt", header=TRUE, row.names=2)

## Ensure all sample metadata names are in the raw data and vice versa
#length(row.names(samples)) == length(names(counts_h))
#stopifnot(all(row.names(samples) %in% names(counts_h)))
#stopifnot(all(names(counts_h) %in% row.names(samples)))
## Reorder the sample data to be in the same order as the count data matrix
#samples <- samples[names(counts_h),]

## Sanity check that sample are in the same order
## This should return TRUE for all of these
#stopifnot(all(names(TPM_h) == row.names(samples)))
#stopifnot(all(names(TPM_m) == row.names(samples)))

#TPM_h2 <- TPM_h
#TPM_m2 <- TPM_m
## Now set the col names of each data matrix to the rownames of the samples matrix
#names(TPM_h2) <- samples$preferred.sample.name
#names(TPM_m2) <- samples$preferred.sample.name

######
## Now add in the Gene Symbols as the first column
df_human_gene <- genes(EnsDb_h, return.type = "DataFrame")[,c("gene_id","symbol")]
df_mouse_gene <- genes(EnsDb_m, return.type = "DataFrame")[,c("gene_id","symbol")]

TPM_h2_annot <- merge(df_human_gene, TPM_h, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)
TPM_m2_annot <- merge(df_mouse_gene, TPM_m, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)

counts_h_annot <- merge(df_human_gene, counts_h, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)
counts_m_annot <- merge(df_mouse_gene, counts_m, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)


## Now save the files!
write.table(TPM_h2_annot, file="processed_data/TPM_human_AllGenes_2021_BulkRNASeq_10.28.21.tsv", quote=FALSE, sep="\t", row.names=FALSE)
write.table(TPM_m2_annot, file="processed_data/TPM_mouse_AllGenes_2021_BulkRNASeq_10.28.21.tsv", quote=FALSE, sep="\t", row.names=FALSE)
write.table(counts_h_annot, file="processed_data/Count_human_AllGenes_2021_BulkRNASeq_10.28.21.tsv", quote=FALSE, sep="\t", row.names=FALSE)
write.table(counts_m_annot, file="processed_data/Count_mouse_AllGenes_2021_BulkRNASeq_10.28.21.tsv", quote=FALSE, sep="\t", row.names=FALSE)



##################### Process TPM values for Analyses ###################

## Remove all genes with zero expression across all samples
counts_h_noZero <- filter.removeZeroSumRows(counts_h)
TPM_h_noZero <- TPM_h[row.names(counts_h_noZero),]

counts_m_noZero <- filter.removeZeroSumRows(counts_m)
TPM_m_noZero <- TPM_m[row.names(counts_m_noZero),]

## convert to preferred sample names


## Load in sample meta data
#samples <- read.delim(file="~/Desktop/CCTR/Data/ChuckHarrell_BrainMetastasis_12-2016/8-18-17_Salmon_RawCountDataFreeze/Sample_Metadata_Freeze_8-18-17_updated_11-27-17.txt", header=TRUE, row.names=2)

## Ensure all sample metadata names are in the raw data and vice versa
#length(row.names(samples)) == length(names(counts_h_noZero))
#stopifnot(all(row.names(samples) %in% names(counts_h_noZero)))
#stopifnot(all(names(counts_h_noZero) %in% row.names(samples)))
## Reorder the sample data to be in the same order as the count data matrix
#samples <- samples[names(counts_h_noZero),]

### Ensure all sample metadata names are in the raw data without filtering and vice versa
#length(row.names(samples)) == length(names(counts_h))
#stopifnot(all(row.names(samples) %in% names(counts_h)))
#stopifnot(all(names(counts_h) %in% row.names(samples)))
### Reorder the sample data to be in the same order as the count data matrix
#samples <- samples[names(counts_h),]
### Sanity check
#stopifnot(all(names(counts_h) == row.names(samples)))
#stopifnot(all(names(TPM_h) == row.names(samples)))
### Now set the col names of each data matrix to the rownames of the samples matrix
#names(counts_h) <- samples$preferred.sample.name
#names(TPM_h) <- samples$preferred.sample.name


## Sanity check that sample are in the same order
## This should return TRUE for all of these
#stopifnot(all(names(counts_h_noZero) == row.names(samples)))
#stopifnot(all(names(counts_m_noZero) == row.names(samples)))
#stopifnot(all(names(TPM_h_noZero) == row.names(samples)))
#stopifnot(all(names(TPM_m_noZero) == row.names(samples)))

## Now set the col names of each data matrix to the rownames of the samples matrix
#names(counts_h_noZero) <- samples$preferred.sample.name
#names(counts_m_noZero) <- samples$preferred.sample.name
#names(TPM_h_noZero) <- samples$preferred.sample.name
#names(TPM_m_noZero) <- samples$preferred.sample.name

## Now add in the Gene Symbols as the first column
#df_human_gene <- genes(EnsDb_h, return.type = "DataFrame")[,c("gene_id","symbol")]
#df_mouse_gene <- genes(EnsDb_m, return.type = "DataFrame")[,c("gene_id","symbol")]

counts_h_noZero_annot <- merge(df_human_gene, counts_h_noZero, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)
counts_m_noZero_annot <- merge(df_mouse_gene, counts_m_noZero, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)
TPM_h_noZero_annot <- merge(df_human_gene, TPM_h_noZero, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)
TPM_m_noZero_annot <- merge(df_mouse_gene, TPM_m_noZero, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)

#TPM_h_annot <- merge(df_human_gene, TPM_h, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)



TPMlog2_h_noZero <- log2(TPM_h_noZero+1)
TPMlog2_m_noZero <- log2(TPM_m_noZero+1)
TPMlog2_h_noZero_annot <- merge(df_human_gene, TPMlog2_h_noZero, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)
TPMlog2_m_noZero_annot <- merge(df_mouse_gene, TPMlog2_m_noZero, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)

#TPMlog2_h <- log2(TPM_h+1)
#TPMlog2_h_annot <- merge(df_human_gene, TPMlog2_h, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)

processed_data/Count_mouse_AllGenes_2021_BulkRNASeq_10.28.21.tsv
## Now save the files!
write.table(counts_h_noZero_annot, file="processed_data/Count_human_ZeroExpRemoved_2021_BulkRNASeq_10.28.21.tsv", quote=FALSE, sep="\t", row.names=FALSE)
write.table(counts_m_noZero_annot, file="processed_data/Count_mouse_ZeroExpRemoved_2021_BulkRNASeq_10.28.21.tsv", quote=FALSE, sep="\t", row.names=FALSE)
write.table(TPM_h_noZero_annot, file="processed_data/TPM_human_ZeroExpRemoved_2021_BulkRNASeq_10.28.21.tsv", quote=FALSE, sep="\t", row.names=FALSE)
write.table(TPM_m_noZero_annot, file="processed_data/TPM_mouse_ZeroExpRemoved_2021_BulkRNASeq_10.28.21.tsv", quote=FALSE, sep="\t", row.names=FALSE)


## Also save the Log2 Transformed TPM values
write.table(TPMlog2_h_noZero_annot, file="processed_data/log2TPM_human_ZeroExpRemoved_2021_BulkRNASeq_10.28.21.tsv", quote=FALSE, sep="\t", row.names=FALSE)
write.table(TPMlog2_m_noZero_annot, file="processed_data/log2TPM_mouse_ZeroExpRemoved_2021_BulkRNASeq_10.28.21.tsv", quote=FALSE, sep="\t", row.names=FALSE)


## Upper Quantile Normalize the Log2 values using only expressed genes.
## this only uses data that is expressed.  Otherwise we end up dividing by zero for some samples.
TPMlog2_h_noZero.quantileAll <- apply(TPMlog2_h_noZero, 2, function(x){quantile(x[x>0], 0.75)})
TPMlog2_h_noZero.norm <- as.data.frame(t(t(TPMlog2_h_noZero) / TPMlog2_h_noZero.quantileAll))

TPMlog2_m_noZero.quantileAll <- apply(TPMlog2_m_noZero, 2, function(x){quantile(x[x>0], 0.75)})
TPMlog2_m_noZero.norm <- as.data.frame(t(t(TPMlog2_m_noZero) / TPMlog2_m_noZero.quantileAll))

TPMlog2_h_noZero.norm_annot <- merge(df_human_gene, TPMlog2_h_noZero.norm, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)
TPMlog2_m_noZero.norm_annot <- merge(df_mouse_gene, TPMlog2_m_noZero.norm, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)

## Also save the Log2 Transformed TPM values
write.table(TPMlog2_h_noZero.norm_annot, file="processed_data/log2TPM_human_ZeroExpRemoved_QuantileNorm_2021_BulkRNASeq_10.28.21.tsv", quote=FALSE, sep="\t", row.names=FALSE)
write.table(TPMlog2_m_noZero.norm_annot, file="processed_data/log2TPM_mouse_ZeroExpRemoved_QuantileNorm_2021_BulkRNASeq_10.28.21.tsv", quote=FALSE, sep="\t", row.names=FALSE)




## Normalize the human and mouse genes by human and mouse percent 
##### UPDATE 6/13/18: We don't want to normalize by human percent anymore.
#TPMLog2Norm_humanPercentCorrected <- removeBatchEffect(TPMlog2_h_noZero, covariates=samples$human.percent, design=model.matrix(~samples$tissue+samples$PDX.line))
#TPMLog2Norm_mousePercentCorrected <- removeBatchEffect(TPMlog2_m_noZero, covariates=samples$mouse.percent, design=model.matrix(~samples$tissue+samples$PDX.line))

#TPMLog2Norm_humanPercentCorrected_notFiltered <- removeBatchEffect(TPMlog2_h, covariates=samples$human.percent, design=model.matrix(~samples$tissue+samples$PDX.line))
#TPMlog2_h_corrected_annot <- merge(df_human_gene, TPMLog2Norm_humanPercentCorrected_notFiltered, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)


#TPMlog2_h_noZero_corrected_annot <- merge(df_human_gene, TPMLog2Norm_humanPercentCorrected, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)
#TPMlog2_m_noZero_corrected_annot <- merge(df_mouse_gene, TPMLog2Norm_mousePercentCorrected, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)

## Save the Log2 Transformed and corrected TPM values
#write.table(TPMlog2_h_noZero_corrected_annot, file="log2TPM_human_ZeroExpRemoved_HumanPercentCorrected_DataFreeze_8.18.17_updated_11.27.17.txt", quote=FALSE, sep="\t", row.names=FALSE)
#write.table(TPMlog2_m_noZero_corrected_annot, file="log2TPM_mouse_ZeroExpRemoved_MousePercentCorrected_DataFreeze_8.18.17_updated_11.27.17.txt", quote=FALSE, sep="\t", row.names=FALSE)

#write.table(TPMlog2_h_corrected_annot, file="log2TPM_human_NotFiltered_HumanPercentCorrected_DataFreeze_8.18.17_updated_04.11.18.txt", quote=FALSE, sep="\t", row.names=FALSE)



################## This section is no longer needed.  The Unmapped reads are no longer in the percent calculations.#################
###### I need to add in the unmapped percent to the metadata file and also use that in the covariates.  It needs to be a matrix with 2 columns.
#human_cov = samples[,c("human.percent","unmapped.percent"), drop=FALSE]
#mouse_cov = samples[,c("mouse.percent","unmapped.percent"), drop=FALSE]

## Normalize the human and mouse genes by human and mouse percent with the unmapped read percent as well.
#TPMLog2Norm_humanUnmappedPercentCorrected <- removeBatchEffect(TPMlog2_h_noZero, covariates=human_cov, design=model.matrix(~samples$tissue+samples$PDX.line))
#TPMLog2Norm_mouseUnmappedPercentCorrected <- removeBatchEffect(TPMlog2_h_noZero, covariates=mouse_cov, design=model.matrix(~samples$tissue+samples$PDX.line))

#TPMlog2_h_noZero_corrected2_annot <- merge(df_human_gene, TPMLog2Norm_humanUnmappedPercentCorrected, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)
#TPMlog2_m_noZero_corrected2_annot <- merge(df_mouse_gene, TPMLog2Norm_mouseUnmappedPercentCorrected, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)

## Save the Log2 Transformed and corrected TPM values
#write.table(TPMlog2_h_noZero_corrected2_annot, file="log2TPM_human_ZeroExpRemoved_HumanUnmappedPercentCorrected_DataFreeze_8-18-17.txt", quote=FALSE, sep="\t", row.names=FALSE)
#write.table(TPMlog2_m_noZero_corrected2_annot, file="log2TPM_mouse_ZeroExpRemoved_MouseUnmappedPercentCorrected_DataFreeze_8-18-17.txt", quote=FALSE, sep="\t", row.names=FALSE)


