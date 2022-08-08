library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library("berryFunctions")
library(collections)
library(tibble)
library(janitor)
library(Metrics)
library(reshape)
library(scales)
library('tidyr')

exTregSigs <- as.data.frame(read.csv("/Users/jmakings/Desktop/GSEA_Antoine/Gene_Signatures.gmx", sep="\t",header = TRUE, row.names = 1))
exTregSigs <- rownames_to_column(exTregSigs, 'exTregs')

# significant overlap between top 100 
intersect(exTregSigs$exTregs[1:50],exTregSigs$Tregs[1:50])

geneExpression <- as.data.frame(read.csv("/Users/jmakings/Desktop/GSEA_Antoine/GSEA_cluster_matrix.gct", sep="\t",header = TRUE, row.names = 1))

geneExp <- geneExpression[-c(1),]
geneExp <- row_to_names(geneExp,1)
geneExp <- subset(geneExp, select = -c(1))







