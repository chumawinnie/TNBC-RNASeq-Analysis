#install.packages('ballgown')
#install.packages('RColorBrewer')
#install.packages('genefilter')
#install.packages('dplyr')
#install.packages('devtools')

##Required for Linux users
#sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev


BiocManager::install("ballgown")

library(ballgown) #for differential expression analysis
library(RColorBrewer) #or RSkittleBrewer for plot coloring
library(genefilter) #calc of means and variances
library(dplyr) #sort and arrange results
#library(devtools) #to ensure reproducibility and install relevant packages
library(gplots)
library(plyr)


pheno_data <- read.csv("treatment.csv", sep=",")

bg_data <- ballgown(dataDir = "E:/ballgown/cvrg", samplePattern = "SRR", pData=pheno_data)

#rowVars calculated vatiance of each row
bg_filt = subset(bg_data,"rowVars(texpr(bg_data)) >1",genomesubset=TRUE) #subset subset only specific samples og genomic locations
bg_table = texpr(bg_filt, 'all') #texpr (extract transcrip-level expressions from bg objects )
bg_gene_names = unique(bg_table[, 9:10])  #takes unique values(all rows) from column 9 and 10 which gene id and gene name

gene_expression = as.data.frame(gexpr(bg_filt))
colnames(gene_expression) <- gsub(pattern = "FPKM.", replacement = "", x = names(gene_expression))


#results_transcripts <- stattest(bg_filt, feature="transcript", covariate="phenotype", getFC=TRUE, meas="FPKM", adjustvars = NULL)
results_genes <- stattest(bg_filt, feature="gene", covariate="treatment", getFC=TRUE, meas="FPKM",  adjustvars = NULL)

##optional --> isoform switching 
#results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg_chrX_filt), geneIDs=ballgown::geneIDs(bg_chrX_filt), results_transcripts)
#results_transcripts <- arrange(results_transcripts,pval)

write.csv(bg_gene_names, "bg_gene_names.csv")
write.csv(bg_table, "bg_table.csv")
write.csv(results_genes, "results_genes.csv")
#write.csv(results_transcripts, "transcripts.csv")

#results_transcriptsmer = data.frame(geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt), transcriptNames=ballgown::transcriptNames(bg_filt), results_transcripts)
results_genes <- arrange(results_genes,qval)

#results_transcripts$logFC <- log2(results_transcripts$fc)
results_genes$logFC <- log2(results_genes$fc) #taking log of Fold Chanage Values

#write.csv(results_transcripts, "transcript_results.csv",row.names=FALSE)
write.csv(results_genes, "gene_results_log_fc_q_value_arranged.csv", row.names=FALSE) 

#subset(results_transcripts,results_transcripts$qval<0.05)
#subset(results_genes,results_genes$pval<0.05)

# subset of result transcripts where pvalue is less than 0.05
#tra <- subset(results_transcripts,results_transcripts$pval<0.05)
gen <- subset(results_genes,results_genes$qval<0.05) #confirm it first

#write.csv(tra, "filtered transcripts.csv")
write.csv(gen, "filtered_genes_qval_0.05.csv") #confirm it first


#####################Box plot#######################
tropical= c('darkorange', 'dodgerblue',
            'hotpink', 'limegreen', 'yellow')
palette(tropical)
#extract transcript-level expression measurements from ballgown objects
fpkm = texpr(bg_filt,meas="FPKM")
fpkm = log2(fpkm+1)
boxplot (fpkm,col=as.numeric(pheno_data$phenotype),las=2,ylab='log2(FPKM+1)')
 

#plotTranscripts(ballgown::geneIDs(bg_data)[ballgown::geneNames(bg_data) == "PCL9"], bg_data, main=c('Gene TPRG1L in sample SRR066423'), sample=c('SRR14585769'))


my_palette <- colorRampPalette(c("green", "yellow", "red"))(n = 299)
mat_data <- data.matrix(fpkm)
heatmap.2(mat_data,
          main = "Heatmap Between Samples Dex vs. Stv",
          notecol="black",
          # density.info="none",
          trace="none",
          margins =c(12,9),
          col=my_palette,
          # breaks=col_breaks,
          dendrogram="both")

write.csv(gene_expression, "raw_gene_expression_results.csv", row.names=FALSE)
write.csv(merge(results_genes, gene_expression, by=0, all=T), "merged_all_results.csv", row.names=TRUE)
#write.csv(results_transcripts, "transcript_results.csv", row.names=FALSE)
write.csv(results_genes, "final_logFCgene_results.csv", row.names=FALSE) 


transcript_gene_table = indexes(bg_filt)$t2g
head(transcript_gene_table)
#Each row of data represents a transcript. Many of these transcripts represent the same gene. Determine the numbers of transcripts and unique genes
length(row.names(transcript_gene_table)) 
length(unique(transcript_gene_table[,"g_id"]))

counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)
#Plot #2 - the distribution of transcript sizes as a histogram
full_table <- texpr(bg_filt , 'all')
hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")


#########################################
sig=which(results_genes$qval<0.05)
hist(results_genes[sig,"logFC"], breaks=50, col="seagreen", xlab="log2(Fold change Stv vs. Dex)", main="Distribution of DEGs Values in Stv vs. Dex")
abline(v=-1, col="black", lwd=2, lty=2)
abline(v=1, col="black", lwd=2, lty=2)
legend("topright", "Fold-change > 1", lwd=2, lty=2)

################################# enhanced volcano###########################
library(EnhancedVolcano)

indices <- match(results_genes$id, texpr(bg_data, 'all')$gene_id)
gene_names_for_result <- texpr(bg_data, 'all')$gene_name[indices]
final_data_Frame <- data.frame(geneNames=gene_names_for_result, results_genes)
write.csv(final_data_Frame,"deg_result_genes_with_gene_names.csv")

final_without_dots <- subset(final_data_Frame, geneNames!=".")

EnhancedVolcano(final_without_dots,
                lab = final_without_dots$geneNames,
                x = "logFC",
                y = "pval",
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Stv vs. Dex")

#trans <- subset(results_transcriptsmer, pval < 0.05)
#trans <- subset(results_transcriptsmer, fc < -0.5 | fc > 0.5)
#final_set <- subset(results_transcriptsmer, pval < 0.05)
final_set <- subset(final_without_dots, logFC < -1 | logFC > 1 & qval < 0.05)
write.csv(final_set,"final_set_with_cutoff_0.5.csv")
#write.csv(final_without_dots,"final_without_dots_0.5.csv")
