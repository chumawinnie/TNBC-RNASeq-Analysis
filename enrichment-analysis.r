data <- read.csv("final_set_with_cutoff_0.5.csv", header = TRUE)#reading signifacnt genes files

#applying parameter for upregulated and download egenes that pass threshold 
up.idx <- which(data$qval < 0.05 & data$logFC > 1)
dn.idx <- which(data$qval < 0.05 & data$logFC < -1)

#getting gene names from the data
all.genes <- data$geneNames
up.genes <- data[up.idx,]$geneNames
dn.genes <- data[dn.idx,]$geneNames


library(enrichR) #loading enrichR library
db <- listEnrichrDbs() #listing all enrichr databases
setEnrichrSite("YeastEnrichr") #setting the enrichR site to yeast
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", 
         "GO_Biological_Process_2018", "KEGG_2018") #list of databases that we are employing



BP <- "GO_Biological_Process_2018" #naming the database to be used 
enriched <- enrichr(up.genes, BP) #passing up/dn genes to enrichr function with the database name, bp/cc/mf/kegg
plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value", title = "Upregulated Genes - GO Biological Process")
MF <- "GO_Molecular_Function_2018"
enriched <- enrichr(up.genes, MF)
plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value", title = "Upregulated Genes - GO Molecular Function")
CC <- "GO_Cellular_Component_2018"
enriched <- enrichr(up.genes, CC)
plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value", title = "Upregulated Genes - GO Cellular Component")
KEGG <- "KEGG_2018"
enriched <- enrichr(up.genes, KEGG)
plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value", title = "Upregulated Genes - KEGG Pathways")



BP <- "GO_Biological_Process_2018"
enriched <- enrichr(dn.genes, BP)
plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value", title = "Downregulated Genes - GO Biological Process")
MF <- "GO_Molecular_Function_2018"
enriched <- enrichr(dn.genes, MF)
plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value", title = "Downregulated Genes - GO Molecular Function")
CC <- "GO_Cellular_Component_2018"
enriched <- enrichr(dn.genes, CC)
plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value", title = "Downregulated Genes - GO Cellular Component")
KEGG <- "KEGG_2018"
enriched <- enrichr(dn.genes, KEGG)
plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value", title = "Downregulated Genes - KEGG Pathways")

