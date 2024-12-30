library(sleuth) ##load the sleuth package 
library(dplyr) #load dplyr package 
#read metadata file, that contain sample information and conditions
metadata<- read.csv("metadata.csv", header = TRUE) 
##add path variable that contain path to kallisto abundance.h5 file
metadata <- mutate(metadata, path = file.path("SS",Run, "abundance.h5"))
## Rename the Run variable (Contain SRR number of sample) in metadata file to sample 
## because sleuth_prep require sample column.
metadata <- dplyr::rename(metadata, sample = Run)
head(metadata)
## This method takes a list of samples with kallisto results and returns a sleuth object
## with the defined normalization of the data across samples
so <- sleuth_prep(metadata, ~Genotype)
## This behavior will fit the full model initially specified and store it in the
## sleuth object under 'full'.
so <- sleuth_fit(so, ~Genotype, 'full')
## This function computes the Wald test on one specific 'beta' coefficient on every transcript.
so_wt<- sleuth_wt(so, 'GenotypeWT')
## This function extracts Wald test results from a sleuth object
results <- sleuth_results(so_wt, 'GenotypeWT', test_type = "wt")
head(results, 20)
## A volcano plot is a p
## Ideally, it looks like a volcano; more significance typically results in higher beta
sleuth::plot_volcano(so_wt, 'GenotypeWT')
## Make an 'MA plot' for a given test. MA plots display, for each transcript, 
## the mean of abundances across samples on the x-axis and fold change on the y-axis.
plot_ma(so_wt, 'GenotypeWT', 'wt')

## Save the sleuth results into CSV file. 
write.csv(results, 'sleuth_results.csv')

##replace code here. principal component graphs showing clustering and divergence of genes
plot_pca(so_wt, text_labels = TRUE, color_by = "Assay.Type")
#heat map also showing clusters,based on genotype and treatment
plot_sample_heatmap(so_wt,annotation_cols = setdiff(colnames(so_wt$sample_to_covariates), 'Genotype'))

## Download the s.cerevisiae_gene_ensembl information using biomaRt
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "scerevisiae_gene_ensembl",
                         host = "www.ensembl.org")
# host = "ensembl.org")
## Get require attribute and filter refseq DNA according to sleuth target id
extracted_results <- biomaRt::getBM(
  attributes = c( "external_gene_name", "description",
                  "transcript_biotype", "refseq_dna"),"refseq_dna",results$target_id,  
  mart = mart)
## Save CSV file of s.cerevisiae_gene_ensembl information that is filter according sleuth IDs
write.csv(extracted_results, "Sleuth genes.csv")
