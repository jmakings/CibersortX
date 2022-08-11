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
library(HSMMSingleCell)
library(SingleCellExperiment)
library('scater')
library(destiny)
library(BisqueRNA)
library('CATALYST')
library(plotly)
library("scatterplot3d")
library(monocle)

# load in gene signatures matrix for exTreg, Treg, Th1, Tn
exTregSigs <- as.data.frame(read.csv("/Users/jmakings/Desktop/GSEA_Antoine/Gene_Signatures.gmx", sep="\t",header = TRUE, row.names = 1))
exTregSigs <- rownames_to_column(exTregSigs, 'exTregs')

# significant overlap between top 100 in the gene sets for exTregs and Tregs
intersect(exTregSigs$exTregs[1:50],exTregSigs$Tregs[1:50])

# load in cluster information
geneExpression <- as.data.frame(read.csv("/Users/jmakings/Desktop/GSEA_Antoine/GSEA_cluster_matrix.gct", sep="\t",header = TRUE, row.names = 1))

geneExp <- geneExpression[-c(1),]
geneExp <- row_to_names(geneExp,1)
geneExp <- subset(geneExp, select = -c(1))

sort7 <- geneExp %>% arrange(desc(CD4T_7))
sort7 <- rownames(sort7)

sort17 <- geneExp %>% arrange(desc(CD4T_17))
sort17 <- rownames(sort17)

# intersection between top expressed Cluster 7 and Cluster 17 genes
length(intersect(sort17[1:50],sort7[1:50]))

# 31 shared genes from the top 50 of both Cluster 7 and 17

length(intersect(sort17[1:25],sort7[1:25]))
# 15 between the top 25 

length(intersect(sort17[1:10],sort7[1:10]))
# 6 between the top 10


# Trajectory analysis for CAVA

seuratcava@assays$RNA@counts
DimPlot(object =  seuratcava)

seuratcava <- FindVariableFeatures(seuratcava)
seuratcava <- ScaleData(seuratcava)
seuratcava <- RunPCA(seuratcava)
seuratcava <- RunUMAP(seuratcava, dims = 1:10)
DimPlot(seuratcava, reduction='umap')

rna <- read.csv("/Users/jmakings/Desktop/CibersortX/raw/raw.csv", sep="\t",header = TRUE, row.names = 1)
adt <- read.csv("/Users/jmakings/Desktop/CibersortX/raw/raw_ab.csv", sep="\t",header = TRUE, row.names = 1)
rna <- CollapseSpeciesExpressionMatrix(rna)

seuratcava <- CreateSeuratObject(counts = rna, project = "CAVA")

adt_assay <- CreateAssayObject(counts = adt)
seuratcava[['ADT']] <- adt_assay

# perform visualization and clustering steps
seuratcava <- NormalizeData(seuratcava)
seuratcava <- FindVariableFeatures(seuratcava)
seuratcava <- ScaleData(seuratcava)
seuratcava <- RunPCA(seuratcava, verbose = FALSE)
seuratcava <- FindNeighbors(seuratcava, dims = 1:30)
seuratcava <- FindClusters(seuratcava, resolution = 0.8, verbose = FALSE)
seuratcava <- RunUMAP(seuratcava, dims = 1:30)
DimPlot(seuratcava, label = TRUE)

# figure this out later 


# Destiny trajectory analysis

# importing Sujit's Seurat object

CD4_cava <- readRDS("/Users/jmakings/Desktop/CibersortX/AntoineGSEA/CD4T_CAVA.rds")

# converting to expression set 
CD4_es <- SeuratToExpressionSet(CD4_cava)
CD4_es

exprs(CD4_es)[1:9,1:9]

dm <- DiffusionMap(CD4_es, n_pcs = 50)

plot(dm, 2:3,legend_main = 'Cluster')

destiny1 <- qplot(DC1,DC3,data = dm2)


set.seed(2)
# converting to single cell experiment 
CD4_ss <- as.SingleCellExperiment(CD4_cava)

matrix <- as.data.frame(assay(CD4_ss, i='logcounts'))

dm2 <- DiffusionMap(t(matrix), n_pcs = 50)
reducedDim(CD4_ss, type = 'DC') <- dm2@eigenvectors

destiny2<- plotReducedDim(CD4_ss,dimred = 'DC',colour_by = 'Cluster')

destiny3 <- plotReducedDim(CD4_ss, dimred = 'DC', ncomponents = 2:3,colour_by = 'Cluster')

# maybe filter for the clusters that we want? 

reduced <- filterSCE(CD4_ss, CD4_ss$Cluster=='CD4T_7' | CD4_ss$Cluster=='CD4T_17' |
                       CD4_ss$Cluster=='CD4T_1' | CD4_ss$Cluster=='CD4T_2' | CD4_ss$Cluster=='CD4T_14' )

plotReducedDim(reduced, dimred = 'DC', ncomponents = 1:2,colour_by = 'Cluster')

seven17 <- filterSCE(CD4_ss, CD4_ss$Cluster=='CD4T_7' | CD4_ss$Cluster=='CD4T_17' )

plotReducedDim(seven17, dimred = 'DC', ncomponents = 1:2,colour_by = 'Cluster')

plotReducedDim(seven17, dimred = 'DC', ncomponents = 3,colour_by = 'Cluster')

dpt <- DPT(dm2)

branchPlot <- plot(dpt, root = 2, paths_to = c(1,3), col_by = 'branch')

clusterPlot <- plot(dpt, root = 2, paths_to = c(1,3), col_by = 'cluster')

dptTimePlot <- ggplot(as.data.frame(colData(CD4_ss)),
       aes(x = dpt,
           y = Cluster, color = Cluster)) + geom_point() + theme_classic() +
  xlab("DPT") + ylab("Cluster") +
  ggtitle("Cells ordered by DPT")

ggplot(as.data.frame(colData(CD4_ss)),
       aes(x = dpt,
           y = Cluster, color = branch)) + geom_point() + theme_classic() +
  xlab("DPT") + ylab("Cluster") +
  ggtitle("Cells ordered by DPT")

dptTimePlot

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

palette16 <- c("#000000","#999999","#004949","#009292","#ff6db6","#ffb6db",
               "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
               "#920000","#924900","#db6d00","#24ff24","#ffff6d")

CD4_ss$dc1 <- dm2$DC1
CD4_ss$dc2 <- dm2$DC2
CD4_ss$dc3 <- dm2$DC3

reduced <- filterSCE(CD4_ss, CD4_ss$Cluster=='CD4T_7' | CD4_ss$Cluster=='CD4T_17' |
                       CD4_ss$Cluster=='CD4T_1' | CD4_ss$Cluster=='CD4T_2' | CD4_ss$Cluster=='CD4T_14' )

fiveClusters <- ggplot(as.data.frame(colData(reduced)),
       aes(x = dc1,
           y = dc2, color = Cluster)) + geom_point() + theme_classic() +
  xlab("DC1") + ylab("DC2") + scale_colour_manual(values=cbPalette) + 
  ggtitle("Destiny DC Plot by Cluster")
  
allClusters <- ggplot(as.data.frame(colData(CD4_ss)),
         aes(x = dc1,
             y = dc2, color = Cluster)) + geom_point() + theme_classic() +
    xlab("DC1") + ylab("DC2") + scale_colour_manual(values=palette16) + 
  ggtitle("Destiny DC Plot by Cluster")

plot(dpt, col_by='branch', divide=3, dcs=c(3,1,-2), pch=20 )

plotReducedDim(reduced, dimred = 'DC', ncomponents = 3,colour_by = 'Cluster')

axx=list(title='DC3')
axy=list(title='DC2')
axz=list(title='DC1')

plotly1 <- plot_ly(x=reduced$dc3, y=reduced$dc2, z=reduced$dc1, type="scatter3d", mode="markers", color=reduced$Cluster,
        marker= list(size =3)) %>% layout(showlegend=TRUE, scene = list(xaxis = list(title = 'DC3'),
                                                   yaxis = list(title = 'DC2'),
                                                   zaxis = list(title = 'DC1')))

colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73")
colors <- colors[as.numeric(reduced$Cluster)]

dc1<- CD4_ss$dc1
dc2 <- CD4_ss$dc2
dc3 <- CD4_ss$dc3
cluster <- CD4_ss$Cluster
df <- data.frame(dc1,dc2,dc3,cluster)

DC1 <- reduced$dc1
DC2 <- reduced$dc2
DC3 <- reduced$dc3
cluster <- gsub('CD4T_','',reduced$Cluster)
df2 <- data.frame(DC1,DC2,DC3,cluster)

cbPalette <- c("#E5F5F9","#1D91C0","#67001F",
               "#F7FCFD", "#D4B9DA","#ffff6d")
colors <- cbPalette[as.numeric(df2$cluster)]

scatter1 <- scatterplot3d(x=df2$DC3, y=df2$DC1, z=-df2$DC2, pch=16,color=colors, cex.symbols = 0.34)

plotly1
