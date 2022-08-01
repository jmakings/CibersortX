library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library("berryFunctions")

# Load the PBMC dataset
pbmc.rna <- as.sparse(read.csv("/Users/jmakings/Desktop/CAVA_Data/raw.csv", sep=",", header = TRUE, row.names = 1))

metaDF <- data.frame(read.csv("/Users/jmakings/Desktop/CAVA_Data/final_metadata.csv", sep=",", header = TRUE, row.names = 1))

df <- data.frame(pbmc.rna)

# Transpose of dataframe (columns = samples, rows = genes)
df <- t(df)
df <- as.data.frame(df)

# make copy (for safety, compute time is so long)
df2 <- df

# Convert raw counts to CPM for each sample
for (i in 1:ncol(df)) {
  if (i %% 10000 == 0) {
    print(i)
  }
  df[,i] <- (df[,i]/sum(as.vector(df[,i])) * 1000000)
}

# make copy (for safety, compute time is so long)
df2 <- df

# getting the cell types from the metadata
cellTypes <- metaDF$Cell_Type

#cava_genes %in% rownames(bulk)

# setting the cell type as the column names, and deleting the redundant first row of cell types
colnames(df) <- cellTypes

# write final matrix to csv file for cibersort
write.csv(df, "/Users/jmakings/Desktop/CAVA_Data/cpm2.csv", row.names = TRUE)

write.table(df,"/Users/jmakings/Desktop/CAVA_Data/cpm2.txt", sep = "\t",row.names = TRUE)




