library(tximport) ## tximport require for read kallisto abundance file
metadata<- read.csv("treatment.csv", header = TRUE) ## read meatadata file for sample
files <-file.path("SS",metadata$sample, "abundance.h5")
## load multiple file from path SS/sample_name folder/abundance file
names(files) <- paste0(metadata$sample, 1:21) ## Name file according to sample
all(file.exists(files))
## tximport imports transcript-level estimates from various external software and 
## optionally summarizes abundances, counts, and transcript lengths to the gene-level
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
head(txi.kallisto$counts)
names(txi.kallisto)
## Make sample table for conditions 
sampleTable <- data.frame(treatment = metadata$treatment)
## Chnage sample table rownames with sample SRR number 
rownames(sampleTable) <- colnames(txi.kallisto$counts)
library(DESeq2)
## DESeqDataSet is a subclass of RangedSummarizedExperiment, used to store the input values, 
## intermediate calculations and results of an analysis of differential expression.
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~treatment)
## This function performs a default analysis 1- estimate size factor 2- estimate dispersion
## 3- Negative Binomial GLM fitting
dds <- DESeq(dds)
resultsNames(dds)

library(apeglm) ##apeglm provides Bayesian shrinkage estimators
## results extracts a result table from a DESeq analysis giving base means across samples
res <- results(dds, name="treatment_Stv_vs_Dex")
## Adds shrunken log2 fold changes (LFC) and SE to a results table from DESeq run without LFC shrinkage.
res <- lfcShrink(dds, coef="treatment_Stv_vs_Dex", type="apeglm")
res
## Save results into CSV file
write.csv(res, 'DESeq2_kallisto.csv')


library(EnhancedVolcano) #EnhancedVolcano will attempt to fit as many variable names in the plot window as possible
EnhancedVolcano(res,
                lab = rownames(res), ## Remove label by assign space, add label by rownames(res)
                x = 'log2FoldChange',
                y = 'pvalue')
#############################################

#####################################################
library( org.Hs.eg.db) 
library(AnnotationDbi)

