library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library("berryFunctions")
library(collections)

# Part 1: Create a matrix for use on whole blood (includes all major cell types) 

# Load in data and create Seurat object 
wholeBlood <- as.sparse(read.csv("/Users/jmakings/Desktop/CAVA_Data/GSE149938_umi_matrix.csv", sep=",", header = TRUE, row.names = 1))
wholeBlood <- as.data.frame(wholeBlood)
wholeBloodT <- as.data.frame(t(wholeBlood))

kaplanTypes <- as.data.frame(read.csv("/Users/jmakings/Desktop/CAVA_Data/InputFiles/kaplanBulkBig.txt", sep="\t", header = TRUE))
kaplanGenes <- kaplanTypes$Gene

pbmc <- CreateSeuratObject(counts = wholeBloodT, project = "Xie_WholeBlood")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)

plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

pbmc$CellType <- Idents(pbmc)
table(Idents(pbmc))

FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

wholeBloodNames <- rownames(wholeBlood)

# intersection of genes between this dataset and Kaplan's
length(intersect(wholeBloodNames, kaplanGenes))

cellDict <- dict(list("HSC"="Hematopoietic Stem Cell", "MPP"="Multi-Potent Progenitor", 
                             "CMP"= "Common Myeloid Progenitor", "MEP"="Megakaryocyte/erythrocyte Progenitor", 
                      "LMPP"="Lympho-Myeloid Primed Progenitor", "MLP"="Lympho-myeloid Progenitor", 
                      "BNK"="BNK", "GMP"="Granulocyte/macrophage Progenitor", "proB"="proB", "preB"="preB", 
                      "immB"="immatureB", "regB"="regulatory B", "naiB"="Naive B", "memB"="memB", "plasma"="plasma",
                      "CLP"="Common Lymphoid Progenitor", "NKP"="NK Progenitor", "toxiNK"= "cytotoxic NK", 
                      "kineNK"="cytokine NK", "CD4T"="CD4T", "CD8T"="CD8T", "hMDP"= "h monocyte-dendritic cell progenitors", 
                      "cMOP"="Common Monocyte Progenitors", "preM"="pre-monocyte", "claM"="classical Monocyte", 
                      "interM"="intermediate Monocyte","nonM"="nonclassical Monocyte", "proN"="pro-myelocyte", 
                      "myeN"="myelocyte", "metaN"="Meta-myelocyte", "matureN"= "Mature Neutrophil", 
                      "ery"="erythrocyte"))

majCellTypes <- function(pbmc) {
  largerTypes <- list()
  for (i in 1:length(pbmc$orig.ident)) {
    t <- pbmc$orig.ident[i]
    if (t == "HSC" || t == "MPP" || t == "CMP" || t == "MEP" || t == "LMPP" || t == "MLP" || t == "BNK" || t == "GMP") {
      largerTypes <- append(largerTypes, "HSPC")
    } else if (t == "immB" || t == "memB" || t == "naiB" || t == "preB" || t == "proB" || t == "regB" ) {
      largerTypes <- append(largerTypes, "B") 
    } else if (t == "plasma") {
      largerTypes <- append(largerTypes, "plasma")
    } else if (t == "ery") {
      largerTypes <- append(largerTypes, "ery")
    } else if (t == "CLP" || t == "toxiNK" || t == "kineNK" || t == "NKP") {
      largerTypes <- append(largerTypes, "NK")
    } else if (t == "CD4T") {
      largerTypes <- append(largerTypes, 'CD4T') 
    } else if (t == "CD8T") {
      largerTypes <- append(largerTypes, 'CD8T')
    } else if (t == "proN" || t == "myeN" || t == "matureN" || t == "metaN") {
      largerTypes <- append(largerTypes, 'Neutrophil') 
    } else if (t == "hMDP" || t == "cMOP" || t == "preM" || t == "claM" || t == "interM" || t == "nonM") {
      largerTypes <- append(largerTypes, "Monocyte")
    } else {
      print(t)
    }
  }
  return(largerTypes)
}

majorTypes <- majCellTypes(pbmc) 
names(majorTypes) <- colnames(pbmc)
majorTypes <- as.data.frame(majorTypes)
majorTypes <- t(majorTypes)
pbmc <- AddMetaData(pbmc, majorTypes, "MajorTypes")

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

cpm <- counts_To_cpm(wholeBloodT)
cpm <- rbind(c(majorTypes), cpm)
genes <- rownames(cpm)
cpm$Gene <- genes
final <- cpm %>% select(Gene, everything())
final["Gene"][final["Gene"] == 1] <- "Gene"

write.table(final,"/Users/jmakings/Desktop/CAVA_Data/InputFiles/XieSS.txt", sep = "\t",row.names = FALSE, col.names = FALSE)

# This portion to create a matrix without erythrocytes, plasma, HSPCs, and Neutrophils for testing on the CAVA dataset

cpmReduce <- cpm
indices <- list()
for (i in 1:ncol(cpmReduce)) {
  x <- cpmReduce[1,i]
  if (x == "HSPC" || x == "ery" || x == "plasma" || x == "Neutrophil") {
    indices <- append(indices, i)
  }
}
indices <- unlist(indices)
cpmReduce %>% select(-indices) -> cpmReduce

cpmReduce$Gene <- rownames(cpmReduce)
cpmReduce <- cpmReduce %>% select(Gene, everything())
cpmReduce[1,1] <- "Gene"

write.table(cpmReduce,"/Users/jmakings/Desktop/CAVA_Data/InputFiles/XieSSreduce.txt", sep = "\t",row.names = FALSE, col.names = FALSE)

