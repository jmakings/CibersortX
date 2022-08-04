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
library("DESeq2")



# read in Seurat Objects for each cell type
B <- readRDS("WIHS/B.RDS")
CD4 <- readRDS("WIHS/CD4T.RDS")
CD8 <- readRDS("WIHS/CD8T_CellsRm.RDS")
CM <- readRDS("WIHS/CM.RDS")
INT <- readRDS("WIHS/INT.RDS")
NCM <- readRDS("WIHS/NCM.RDS")
NK <- readRDS("WIHS/NK.RDS")

# merge Seurat objects into one large seurat object
wihs <- merge(B, y = c(CD4,CD8,CM,INT,NCM,NK), project = "WIHS")

# to view the number of different called cell types
table(wihs$Type)

