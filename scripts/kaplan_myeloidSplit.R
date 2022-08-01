library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library("berryFunctions")

kaplan <- as.sparse(read.csv("/Users/jmakings/Desktop/CAVA_Data/Kaplan_Data/Kaplan.RSEMv1.3.1_gene_expected_count.txt", sep="\t", header = TRUE, row.names = 1))
kaplan <- CollapseSpeciesExpressionMatrix(kaplan)
kaplan <- data.frame(kaplan)
kaplanT <- t(as.matrix(kaplan))

kapSrt <- CreateSeuratObject(counts = kaplanT, row.names = 0)

annotation <- data.frame(read.csv("/Users/jmakings/Desktop/CAVA_Data/Kaplan_Data/Annotation_File.csv", sep=",", header = TRUE, row.names = 1))

# function converts counts dataframe to counts per million
counts_To_cpm <- function(df) {
  for (i in 1:ncol(df)) {
    if (i %% 25000 == 0) {
      print(paste0("Samples converted: ",i))
    }
    df[,i] <- (df[,i]/sum(as.vector(df[,i])) * 1000000)
  }
  return(df)
}

cpm <- counts_To_cpm(kaplan)

# get transcript id and remove its isoform identifier
transcripts <- row.names(kaplan)
trans2 <- gsub("\\..*", "", transcripts)
cpm$transcript_id <- trans2

# reorder
cpm <- select(cpm, transcript_id, A0,A1,A2,A3,A4,A5,A6,A7,A8,A9)

# reduce annotation file to what we need
anno <- annotation[c("gene_name", "gene_id", "gene_biotype")]

# merge annotations with cpm file
merged <- merge(anno, cpm, by.x = "gene_id", by.y = "transcript_id")

# drop duplicated gene names and non-protein coding genes (I believe this dataset is all protein coding genes)
merged <- merged[!duplicated(merged$gene_name),]
merged <- subset(merged, gene_biotype = "protein_coding")

#sigMtxRef <- data.frame(read.csv("/Users/jmakings/Desktop/CAVA_Data/InputFiles/CAVA_myeloidSplit_ss.txt",sep="\t", header = TRUE, row.names = 1))

genes <- ssFinalT[,1]
genes <- genes[-1]

merg2 <- subset(merged, select = -c(gene_id, gene_biotype))
colnames(merg2)[1] <- "Gene"
merg2final <- merg2[-c(1599),] # this gene undefined

write.table(merg2final, "/Users/jmakings/Desktop/CAVA_Data/InputFiles/kaplanBulkBig.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

merg3 <- subset(merg2, merg2$Gene %in% genes)

write.table(merg3, "/Users/jmakings/Desktop/CAVA_Data/InputFiles/kaplanBulkSmall.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

# some issue with larger file, maybe check later 
