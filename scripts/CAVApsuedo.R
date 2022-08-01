library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library("berryFunctions")

# Jeff way of doing it (ofc did not work)
# Load the PBMC dataset
# pbmc.rna <- as.sparse(read.csv("/Users/jmakings/Desktop/CAVA_Data/raw.csv", sep=",", header = TRUE, row.names = 1))
# pbmc.rna <- CollapseSpeciesExpressionMatrix(pbmc.rna)
# pbmc.rna <- t(as.matrix(pbmc.rna))
# 
# metaDF <- data.frame(read.csv("/Users/jmakings/Desktop/CAVA_Data/final_metadata.csv", sep=",", header = TRUE, row.names = 1))


#Seurat Object Setup 
# pbmc <- CreateSeuratObject(counts=pbmc.rna, row.names=0)
# 
# rownames(pbmc@meta.data)
# order(metaDF$Cell_Index, rownames(pbmc@meta.data))
# rownames(pbmc@meta.data)
# pbmc@meta.data<- cbind(metaDF, pbmc@meta.data)
# pbmc@meta.data
# x <- as.data.frame(pbmc@meta.data)
# 
# select_meta <- x %>% select(Cluster)
# 
# select_meta %>% group_by(Cluster) %>% summarize(n())
# 
# new_meta <- t(select_meta)
# 
# df_merged <- merge(df, select_meta, by=0)
# df_merged %>% dplyr::group_by(Cluster) %>% summarize_all(funs(sum)) -> df_summed

# restart using Sujit's method

# Load in data and create Seurat object 
seurat <- as.sparse(read.csv("/Users/jmakings/Desktop/CAVA_Data/raw.csv", sep=",", header = TRUE, row.names = 1))
seurat <- CollapseSpeciesExpressionMatrix(seurat)
seurat <- t(as.matrix(seurat))

seuratObj <- CreateSeuratObject(counts = seurat, row.names = 0)

# get counts
counts <- seuratObj@assays$RNA@counts 

# dataframe from counts
df <- as.data.frame(seuratObj@assays$RNA@counts) %>% t()

# load metadata and add to Seurat object
metaDF <- data.frame(read.csv("/Users/jmakings/Desktop/CAVA_Data/final_metadata.csv", sep=",", header = TRUE, row.names = 1))
seuratObj@meta.data <- metaDF

# separate variable for metadata
x <- seuratObj@meta.data

# get cluster data
select_meta <- x %>% select(Cluster)

# transfrom "select_meta" to be grouped by cluster
select_meta %>% group_by(Cluster) %>% summarise(n())

# merge cluster data into original dataframe
df_merged <- merge(df, select_meta, by=0)

# summarise transcript counts by cluster
df_merged %>% dplyr::group_by(Cluster) %>% summarise_if(is.numeric, funs(sum)) -> df_summed

df_summedT <- t(df_summed)
df_summedT <- as.data.frame(df_summedT)

df_summedT <- tibble::rownames_to_column(df_summedT, "Gene")
df_summedT["Gene"][df_summedT["Gene"] == "Cluster"] <- "Gene"

names(df_summedT) <- NULL

write.table(df_summedT,"/Users/jmakings/Desktop/CAVA_Data/CAVA_pseudoBulk.txt", sep = "\t",row.names = FALSE)


