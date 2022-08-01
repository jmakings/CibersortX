library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library("berryFunctions")

# Using Sujit's method

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


