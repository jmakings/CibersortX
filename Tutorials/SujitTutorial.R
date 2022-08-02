#Loading packages 
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(slingshot)
library("scater")

#Reading the data
setwd("/Users/sujitsilasarmstrongsuthahar/Desktop/Ley Lab/Ryo_Trajectory_Analysis_CAVA")
pbmc.rna <- as.sparse(read.csv("/Users/sujitsilasarmstrongsuthahar/Desktop/Ley Lab/Ryo_Trajectory_Analysis_CAVA/raw.csv", sep=",", header = TRUE, row.names = 1))
pbmc.rna <- CollapseSpeciesExpressionMatrix(pbmc.rna)
pbmc.rna <- t(as.matrix(pbmc.rna))

pbmc.adt <- as.sparse(read.csv("/Users/sujitsilasarmstrongsuthahar/Desktop/Ley Lab/Ryo_Trajectory_Analysis_CAVA/raw_ab.csv", sep = ",", header = TRUE, row.names = 1))
pbmc.adt <- CollapseSpeciesExpressionMatrix(pbmc.adt)
pbmc.adt <- t(as.matrix(pbmc.adt))
all.equal(colnames(pbmc.adt), colnames(pbmc.adt)) 

#Seurat Object Setup 
pbmc <- CreateSeuratObject(counts=pbmc.rna, row.names=0)

#Adding metadata to the seurat object
pbmc@meta.data
data <- read.csv("/Users/sujitsilasarmstrongsuthahar/Desktop/Ley Lab/Ryo_Trajectory_Analysis_CAVA/final_metadata.csv")
rownames(pbmc@meta.data)
order(data$Cell_Index, rownames(pbmc@meta.data))
rownames(pbmc@meta.data)
pbmc@meta.data<- cbind(data, pbmc@meta.data)
x <- as.data.frame(pbmc@meta.data)


#Create an assay to store the Antibody Derivative Tag (ATG) information
adt_assay <- CreateAssayObject(counts=pbmc.adt)
pbmc[["ADT"]] <- adt_assay
Assays(pbmc)
rownames(pbmc[["ADT"]])

#Clustering cells based on scRNA-seq
DefaultAssay(pbmc) <- 'RNA'
DefaultAssay(pbmc)

#Visualization and Clustering Steps
?NormalizeData()
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", margin = 1, scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, verbose=FALSE)
pbmc <- FindNeighbors(pbmc, dims=1:30)
pbmc <- FindClusters(pbmc, resolution = 0.6, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims=1:30)
DimPlot(pbmc, label = T, group.by = "Group")

?NormalizeData()

#Visualizing multiple modalities side by side
DefaultAssay(pbmc) <- "ADT"
DefaultAssay(pbmc)
pbmc <- NormalizeData(pbmc, normalization.method = "CLR", margin = 2)

pbmc@meta.data

#Slingshot 
sce <- SingleCellExperiment(assays = list(counts=pbmc@assays$ADT))
geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]                     


#
pbmc
#Cluster Heatmaps
max(pbmc@assays[[what_assay]]@scale.data)
min(pbmc@assays[[what_assay]]@scale.data)

max(pbmc@assays[[what_assay]]@scale.data["ACKR3",])
min(pbmc@assays[[what_assay]]@scale.data["ACKR3",])

row.names(pbmc@assays[[what_assay]]@scale.data)
library(scales)
what_assay <- "RNA"
data <-data.frame(t(as.matrix(GetAssayData(pbmc, assay = what_assay, slot = "scale.data"))))
data[1:5, 1:5]
scale_limits <-  c(-2, 2)
# rescale every markers separatly
scaled_data <- data.frame(lapply(data, function(x) rescale(x, to = scale_limits)))
# rescale your data
#scaled_data <- rescale(as.matrix(data), to = scale_limits)
row.names(scaled_data) <- row.names(data)
scaled_data_trans <- t(scaled_data)
scaled_data_trans[1:5, 1:5]
# add the scaled data in Seurat object
pbmc@assays[[what_assay]]@scale.data <- scaled_data_trans
pbmc@assays[[what_assay]]@scale.data[1:5, 1:5]

#Compl
library(ComplexHeatmap)
Heatmap(pbmc@assays$RNA@scale.data)
DoHeatmap(adt.cluster.averages, features = marker_list, size = 3, draw.lines = FALSE,
          group.colors = my_colors,
          assay = “ADT”, slot = what_slot) + scale_fill_gradientn(colours=rev(brewer.pal(11,“RdYlBu”))
pbmc@assays$RNA@scale.data
