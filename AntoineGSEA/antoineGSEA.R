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
#library('CATALYST')
library(plotly)
library("scatterplot3d")
library(monocle)
library(ComplexHeatmap)
library(stringr)
library(reshape)
library(data.table)

##### Load datasets and preprocessing #####

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

##### Exploratory Analysis ##### 
# intersection between top expressed Cluster 7 and Cluster 17 genes
length(intersect(sort17[1:50],sort7[1:50]))

# 31 shared genes from the top 50 of both Cluster 7 and 17

length(intersect(sort17[1:25],sort7[1:25]))
# 15 between the top 25 

length(intersect(sort17[1:10],sort7[1:10]))
# 6 between the top 10


##### DESTINY TRAJECTORY ANALYSIS###### 

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

# different color palettes for colorblind viewing
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

palette16 <- c("#000000","#999999","#004949","#009292","#ff6db6","#ffb6db",
               "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
               "#920000","#924900","#db6d00","#24ff24","#ffff6d")

# loading top DCs into single cell experiment
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

# 3D scatter plot of first 3 DCs 

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

# 8/12/21 Work on changing clusters to gray

grayPalette <- c("gray","gray","gray","gray","gray","gray","gray","gray","blue","gray","gray","gray","gray","gray","red", "gray")

alpha_vec <- vector()
for (x in CD4_ss$Cluster) {
  if (x == "CD4T_7" || x == "CD4T_17") {
    alpha_vec <- append(alpha_vec, 1)
  }
  else 
    alpha_vec <- append(alpha_vec, 0.4)
}

grayPlot <- ggplot(as.data.frame(colData(CD4_ss)),
       aes(x = dc1,
           y = dc2, color = Cluster)) + geom_point(alpha = alpha_vec) + theme_classic() +
  xlab("DC1") + ylab("DC2") + scale_colour_manual(values=grayPalette) + 
  scale_fill_manual(values = alpha(alpha_vec)) + 
  ggtitle("Destiny DC Plot by Cluster")
grayPlot

alpha_vec2 <- vector()
for (x in reduced$Cluster) {
  if (x == "CD4T_7" || x == "CD4T_17") {
    alpha_vec2 <- append(alpha_vec2, 1)
  }
  else 
    alpha_vec2 <- append(alpha_vec2, 0.4)
}

gray5palette <- c("gray","gray", "blue","gray","red")
gray5Plot <- ggplot(as.data.frame(colData(reduced)),
                   aes(x = dc1,
                       y = dc2, color = Cluster)) + geom_point(alpha = alpha_vec2) + theme_classic() +
  xlab("DC1") + ylab("DC2") + scale_colour_manual(values=gray5palette) + geom_abline()
  ggtitle("Destiny DC Plot by Cluster")

gray5Plot

# function to get bin lines, include the plot to put over and the number of bins
plotBinLines <- function(plot, bins) {
  binCD4_ss <- cut(CD4_ss$dc1, breaks=bins)
  binChange <- str_replace_all(binCD4_ss, pattern = ".*,([^-]*)]*", replacement = "\\1")
  binChange2 <- gsub(']','',binChange)
  binChange2 <- sort(as.numeric(unique(binChange2)))
  dif <- binChange2[2] -binChange2[1]
  
  bin20plot <- plot + geom_vline(xintercept = c(unique(binChange2),binChange2[1]-dif))
  return(bin20plot)
}

# plot of bin lines
binPlot <- plotBinLines(grayPlot,10)

# get bins then create to Seurat object
binCD4_ss <- cut(CD4_ss$dc1, breaks=10, labels=c(1:10))
CD4_ss$Bins <- binCD4_ss

newSrt <- as.Seurat(CD4_ss)

# setting ident to cluster, but will later set it to bin and cluster together 
newSrt <- SetIdent(newSrt, value=newSrt@meta.data$Cluster)
newSrt <- SetIdent(newSrt, value = newSrt@meta.data$Bins)

# scaling data (for heatmap later)
all.genes <- rownames(newSrt)
newSrt <- ScaleData(newSrt, features = all.genes)
genes <- rownames(newSrt@assays$RNA@counts)
newSrt <- ScaleData(newSrt, assay = 'RNA', features = genes)

exTregBinsDE <- function(Srt) {
  l <- list()
  for (i in 1:7) {
    print(i)
    l[[i]] <- FindMarkers(subset(x = Srt, Cluster == 'CD4T_7'), ident.1 =i,indent.2=c(8,9,10), min.pct=0, logfc.threshold = 0)
  }
  return(l)
}

#Cluster 7 DE of each bin vs exTreg cluster on the right (Bins 8,9,10)
binsList <- exTregBinsDE(newSrt)

#Cluster 7 DE of each bin vs all other bins
allBins <- FindAllMarkers(subset(x = newSrt, Cluster == 'CD4T_7'), min.pct=0, logfc.threshold = 0)

newSrt <- SetIdent(newSrt, value = newSrt@meta.data$Bins)
FindAllMarkers(newSrt)

#Cluster 7 DE Bins 1-7 vs Bins 8-10
leftvsRight <- FindMarkers(subset(x = newSrt, Cluster == 'CD4T_7'), ident.1 =1:7,indent.2=c(8,9,10), min.pct=0, logfc.threshold = 0)

# add bin and cluster together
newSrt$binCluster <- str_c(newSrt$Cluster, "_Bin_",newSrt$Bins)
counts <- data.frame(t(newSrt@assays$RNA@counts), newSrt$Cell_Index, newSrt$binCluster)

# set it as ident 
newSrt <- SetIdent(newSrt, value = newSrt@meta.data$binCluster)

# finding DEs for each bin/cluster pairing
binClustervsRest <- FindAllMarkers(newSrt, min.pct=0, logfc.threshold = 0,min.cells.feature = 0,min.cells.group = 0)

frame <- as.data.frame(binClustervsRest)

# subset for clusters 7 and 17
DE.7.17 <- subset(binClustervsRest, cluster=='CD4T_7_Bin_1' | cluster== 'CD4T_7_Bin_2' |
                    cluster=='CD4T_7_Bin_3' | cluster=='CD4T_7_Bin_4' | cluster== 'CD4T_7_Bin_5' |
                    cluster=='CD4T_7_Bin_6' | cluster=='CD4T_7_Bin_7' | cluster== 'CD4T_7_Bin_8' |
                    cluster=='CD4T_7_Bin_9' | cluster=='CD4T_7_Bin_10' | cluster== 'CD4T_17_Bin_1' |
                    cluster=='CD4T_17_Bin_2' | cluster== 'CD4T_17_Bin_3' |
                    cluster=='CD4T_17_Bin_4' | cluster== 'CD4T_17_Bin_7' |
                    cluster=='CD4T_17_Bin_8' )

df_bin <- subset(binClustervsRest, gene=="CD56"| gene=="CD16" |gene=="CD25")

#### Less efficient Pseudobulk method, use one below provided by Sujit ####

counts %>% group_by(newSrt$binCluster) %>% summarise_if(is.numeric, funs(sum)) -> countsBulk
countsBulk <- rename(countsBulk, binCluster = 'newSrt$binCluster')
countsBulk <- as.data.frame(countsBulk)


countsBulk <- as.data.frame(t(countsBulk))
countsBulk <- countsBulk %>% row_to_names(row_number = 1)

bulkTreg <- subset(countsBulk[c(69:75,124:132)] )

bulkTreg <- as.data.frame(t(bulkTreg))
# convert this to numeric

bulkTreg <- matrix(as.numeric(bulkTreg), ncol = ncol(bulkTreg))

# converts dataframe to numeric
convert_Numeric <- function(df) {
  newDf <- df
  for (i in colnames(df)) {
    #print(df[i])
    int <- df[i]
    newDf[i] <- as.numeric(unlist(int))
  }
  return(newDf)
}

numBulkTreg <- convert_Numeric(bulkTreg)

longform <- numBulkTreg %>% rownames_to_column() %>% gather(colname, value, -rowname)

ggplot(longform, aes(x = rowname, y = colname, fill = value)) +
  geom_tile() + theme(axis.text=element_text(size=4.5),
                      axis.title=element_text(size=1))

# from counts dataframe, create a bulk from bin and cluster info (use 'fun; to specify function to create bulk)
bulkTransform <- function(counts, newSrt, fun) {
  counts %>% group_by(newSrt$binCluster) %>% summarise_if(is.numeric, funs(fun)) -> countsBulk
  countsBulk <- rename(countsBulk, binCluster = 'newSrt$binCluster')
  countsBulk <- as.data.frame(countsBulk)
  
  countsBulk <- as.data.frame(t(countsBulk))
  countsBulk <- countsBulk %>% row_to_names(row_number = 1)
  
  bulkTreg <- subset(countsBulk[c(69:75,124:132)] )
  
  bulkTreg <- as.data.frame(t(bulkTreg))

  numBulkTreg <- convert_Numeric(bulkTreg)
  return(numBulkTreg)
}

meanTreg <- bulkTransform(counts, newSrt, mean)

longform2 <- meanTreg %>% rownames_to_column() %>% gather(colname, value, -rowname)

ggplot(longform2, aes(x = rowname, y = colname, fill = value)) +
  geom_tile() + theme(axis.text=element_text(size=4.5),
                      axis.title=element_text(size=4))

TregBinsDE <- function(Srt) {
  l <- list()
  for (i in c(2,3,4,7,8,9)) {
    print(i)
    l[[i]] <- FindMarkers(subset(x = Srt, Cluster == 'CD4T_17'), ident.1 =i,indent.2=c(1), min.pct=0, logfc.threshold = 0)
  }
  return(l)
}

TregBins <- TregBinsDE(newSrt)

#### Creating Single Cell heatmaps of DE markers from each bin/cluster####

binClustervsRest %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(newSrt, features = top10$gene)

# subset for just clusters 7 and 17 
Cluster.7.17 <- subset(x = newSrt, idents = c('CD4T_7_Bin_1','CD4T_7_Bin_2', 
                          'CD4T_7_Bin_3','CD4T_7_Bin_4'
                      ,'CD4T_7_Bin_5','CD4T_7_Bin_6', 
                         'CD4T_7_Bin_7', 'CD4T_7_Bin_8', 
                         'CD4T_7_Bin_9', 'CD4T_17_Bin_1', 'CD4T_17_Bin_2',
                         'CD4T_17_Bin_3',  'CD4T_17_Bin_4', 
                         'CD4T_17_Bin_7', 'CD4T_17_Bin_8', 'CD4T_17_Bin_9', 
                      'CD4T_7_Bin_10'))

# order the clusters accordingly
clusterLevels <- c('CD4T_7_Bin_1','CD4T_17_Bin_1','CD4T_7_Bin_2',
                   'CD4T_17_Bin_2','CD4T_7_Bin_3','CD4T_17_Bin_3',
                   'CD4T_7_Bin_4', 'CD4T_17_Bin_4', 'CD4T_7_Bin_5',
                   'CD4T_7_Bin_6', 'CD4T_7_Bin_7','CD4T_17_Bin_7',
                   'CD4T_7_Bin_8', 'CD4T_17_Bin_8','CD4T_7_Bin_9', 
                   'CD4T_17_Bin_9', 'CD4T_7_Bin_10')

# change active idents and assay
Cluster.7.17@active.ident <- factor(x=Cluster.7.17@active.ident, levels = clusterLevels)
Cluster.7.17@active.assay <- 'ADT'

DoHeatmap(Cluster.7.17, features = top10$gene, angle = 90, size = 1.5)

smallerCluster <- subset(x = Cluster.7.17, idents = c('CD4T_7_Bin_1', 'CD4T_7_Bin_2', 
                                                'CD4T_7_Bin_3','CD4T_7_Bin_4'
                                                ,'CD4T_7_Bin_5','CD4T_7_Bin_6', 
                                                'CD4T_7_Bin_7', 'CD4T_7_Bin_8', 
                                                'CD4T_7_Bin_9', 'CD4T_7_Bin_10','CD4T_17_Bin_1', 'CD4T_17_Bin_2',
                                                'CD4T_17_Bin_3'))


DoHeatmap(smallerCluster, features = top10$gene, angle = 90, size = 2)


tregCluster <- subset(x = Cluster.7.17, idents = c('CD4T_17_Bin_1','CD4T_17_Bin_2',
                                               'CD4T_17_Bin_3',  'CD4T_17_Bin_4', 
                                               'CD4T_17_Bin_7', 'CD4T_17_Bin_8', 'CD4T_17_Bin_9'))

DoHeatmap(tregCluster, features = top10$gene, angle = 90, size = 2)

tregClusterminus <- subset(x = Cluster.7.17, idents = c('CD4T_17_Bin_2',
                                             'CD4T_17_Bin_3',  'CD4T_17_Bin_4', 
                                             'CD4T_17_Bin_7', 'CD4T_17_Bin_8', 'CD4T_17_Bin_9'))

DoHeatmap(tregClusterminus, features = top10$gene,size = 2)

exTregCluster <- subset(x = Cluster.7.17, idents = c('CD4T_7_Bin_1','CD4T_7_Bin_2', 
                                               'CD4T_7_Bin_3','CD4T_7_Bin_4'
                                               ,'CD4T_7_Bin_5','CD4T_7_Bin_6', 
                                               'CD4T_7_Bin_7', 'CD4T_7_Bin_8', 
                                               'CD4T_7_Bin_9', 'CD4T_7_Bin_10'))
DoHeatmap(exTregCluster, features = top10$gene, size = 2)


##### Create Pseudobulk/ heatmap from Sujit's way (standard Ley workflow)#####

norm_df <- as.data.frame(newSrt@assays$ADT@data) %>% t()
norm_meta <- newSrt@meta.data %>% select(binCluster)
norm_merge <- merge(norm_df, norm_meta, by=0)
norm_avg <- norm_merge %>% group_by(binCluster) %>% summarise_all(funs(mean)) %>% data.table()
norm_avg$Row.names <- NULL
norm_avg$binCluster <- gsub('CD4T','C',norm_avg$binCluster)
norm_avg <- norm_avg %>% t()
norm_avg <- norm_avg %>% row_to_names(row_number = 1)

write.csv(norm_avg,'AntoineGSEA/AntoineBinClusterMatrix.csv')

normScaled <- CreateSeuratObject(norm_avg)
data <- as.data.frame(t(as.matrix(GetAssayData(normScaled))))
scale_limits <- c(-2,2)
scaled_data <- data.frame(lapply(data,function(x)rescale(x,to=scale_limits)))
rownames(scaled_data) <-  rownames(data)
scaled_data <- t(scaled_data)
normScaled@assays$RNA@scale.data <- scaled_data
normScaled = as.matrix(normScaled@assays$RNA@scale.data)

pal_cols <- colorRampPalette(c("blue","yellow","red"))

normScaled <- as.data.frame(normScaled)

clusters <- c('C_7_Bin_1','C_7_Bin_2', 
              'C_7_Bin_3','C_7_Bin_4'
              ,'C_7_Bin_5','C_7_Bin_6', 
              'C_7_Bin_7', 'C_7_Bin_8', 
              'C_7_Bin_9', 
              'C_7_Bin_10', 'C_17_Bin_1', 'C_17_Bin_2',
              'C_17_Bin_3',  'C_17_Bin_4', 
              'C_17_Bin_7', 'C_17_Bin_8', 'C_17_Bin_9')
scaled.7.17 <- normScaled[co_1]

newnames_row <- lapply(
  rownames(scaled.7.17) ,
  function(x) bquote(bold(.(x))))

newnames_cols <- lapply(
  colnames(scaled.7.17),
  function(x) bquote(bold(.(x))))

pheatmap::pheatmap(scaled.7.17, color = pal_cols(100), cluster_rows = FALSE,
                   main = "Scaled Heatmap for Cluster 7 and 17 Bulk Marker Expression, by Bin", fontsize = 10, fontsize_col = 10, fontsize_row = 10, fontsize_number=10, 
                   labels_row = as.expression(newnames_row),
                   labels_col = as.expression(newnames_cols))

scaled.7.17 <- as.matrix(scaled.7.17)

y = Heatmap(scaled.7.17, column_order = co_1, cluster_rows=T)
g <- column_order(y)
h <- row_order(y)

y = Heatmap(scaled.7.17, column_order = clusters, cluster_rows=T)

##### Filtering for interested clusters, creating table of Cluster 7 pct1 values #####

binClusterGenes <- FindAllMarkers(newSrt, assay = "RNA",min.pct=0, logfc.threshold = 0,min.cells.feature = 0,min.cells.group = 0)

# filter for only the DEGs for the clusters we are interested in

binClusterDF <- as.data.frame(binClusterGenes)

#just cluster7 genes as a function of bin, for Klaus (8/30/22) 
cluster7genes <- binClusterDF[(binClusterDF$cluster=='CD4T_7_Bin_1' |
                                 binClusterDF$cluster=='CD4T_7_Bin_2' |
                                 binClusterDF$cluster=='CD4T_7_Bin_3' |
                                 binClusterDF$cluster=='CD4T_7_Bin_4' |
                                 binClusterDF$cluster=='CD4T_7_Bin_5' |
                                 binClusterDF$cluster=='CD4T_7_Bin_6' |
                                 binClusterDF$cluster=='CD4T_7_Bin_7' |
                                 binClusterDF$cluster=='CD4T_7_Bin_8' |
                                 binClusterDF$cluster=='CD4T_7_Bin_9' |
                                 binClusterDF$cluster=='CD4T_7_Bin_10'),]

# adding log2(pct1/pct2)*log2fc to dataframe
cluster7genes$'log2(pct1/pct2)*log2fc' = log2(cluster7genes$pct.1/cluster7genes$pct.2)*cluster7genes$avg_log2FC

# genes with log2(pct1/pct2)*log2fc == inf, which is invalid
infs <- cluster7genes[is.infinite(cluster7genes$`log2(pct1/pct2)*log2fc`),]

# removing infinite values from df
cluster7genes <- cluster7genes[!is.infinite(cluster7genes$`log2(pct1/pct2)*log2fc`),]

# sort by cluster, pct.1, and log2(pct1/pct2)*log2fc 
cluster7genes <- cluster7genes[order(cluster7genes$cluster, -cluster7genes$pct.1,-cluster7genes$'log2(pct1/pct2)*log2fc'),]

# specify cluster order
cluster7genes <- cluster7genes[order(factor(cluster7genes$cluster,levels=c("CD4T_7_Bin_1", "CD4T_7_Bin_2", "CD4T_7_Bin_3", "CD4T_7_Bin_4", "CD4T_7_Bin_5", 
                                                                     "CD4T_7_Bin_6", "CD4T_7_Bin_7", "CD4T_7_Bin_8", "CD4T_7_Bin_9", "CD4T_7_Bin_10"))),]

notExpressed <- subset(cluster7genes, pct.1 > 0.1 | pct.2 > 0.1)

write.csv(notExpressed,'AntoineGSEA/cluster7_gene_pct.csv', row.names = F)

##### Creating line plot for specified genes over each bin #####
genesForPlot <- c('NKG7','CCL5','CST7','GZMA','HOPX','GZMK','DUSP2','KLRB1','GZMB','FGFBP2','SP1R5')

cluster7genes$bin <- sub('CD4T_7_Bin_', '', cluster7genes$cluster)
cluster7genes <- transform(cluster7genes, bin = as.numeric(bin))

lineplotDF <- subset(cluster7genes, gene %in% genesForPlot)

gene_lineplot <- function(lineplotDF) {
  lineplot <- lineplotDF %>% group_by(bin) %>% ggplot(., aes(x=bin, y=pct.1, color=gene, shape=gene, group=gene)) + geom_point(size=5)  + geom_line(size=1)+
    theme(axis.title = element_text(face="bold", size = 14)) +
    theme(legend.title=element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90),
          panel.grid.minor = element_blank())+ scale_x_continuous(breaks = round(seq(min(lineplotDF$bin), max(lineplotDF$bin), by = 1),1)) + 
    scale_shape_manual(values=seq(0,20))+ 
    scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73"
                                 , "#0072B2", "#D55E00", "#CC79A7", "#DDCC77", "#117733", "#332288", "#AA4499"))+
    labs(title = "PCT1 of selected genes in Cluster 7 across Bins", x ='Bin', y = "pct1") +
    theme(title = element_text(size=15, face = "bold"), legend.title = element_text(size=0), legend.text=element_text(size=15), 
          axis.line =element_line(size = 1), axis.ticks = element_line(size=1), axis.text = element_text(size=12))
  
  return(lineplot) 
}

lineplot <- gene_lineplot(lineplotDF)


##### Differential Expression/single cell heatmap for each gene #####

filtered <- binClusterDF[(binClusterDF$cluster=='CD4T_7_Bin_1' |
                            binClusterDF$cluster=='CD4T_7_Bin_2' |
                            binClusterDF$cluster=='CD4T_7_Bin_3' |
                            binClusterDF$cluster=='CD4T_7_Bin_4' |
                            binClusterDF$cluster=='CD4T_7_Bin_5' |
                            binClusterDF$cluster=='CD4T_7_Bin_6' |
                            binClusterDF$cluster=='CD4T_7_Bin_7' |
                            binClusterDF$cluster=='CD4T_7_Bin_8' |
                            binClusterDF$cluster=='CD4T_7_Bin_9' |
                            binClusterDF$cluster=='CD4T_7_Bin_10' |
                            binClusterDF$cluster=='CD4T_17_Bin_1' |
                            binClusterDF$cluster=='CD4T_17_Bin_2' |
                            binClusterDF$cluster=='CD4T_17_Bin_3' |
                            binClusterDF$cluster=='CD4T_17_Bin_4' |
                            binClusterDF$cluster=='CD4T_17_Bin_7' |
                            binClusterDF$cluster=='CD4T_17_Bin_8' |
                            binClusterDF$cluster=='CD4T_17_Bin_9'),]


filtered %>% 
  group_by(cluster) %>% 
  top_n(n = 15, wt = avg_log2FC) -> topgenes

DefaultAssay(Cluster.7.17) <- 'RNA'

# all of cluster 7 and 17
DoHeatmap(Cluster.7.17, features = topgenes$gene, angle = 90, size = 2) +
  theme(text=element_text(size = 7), legend.text = element_text(size = 10))

DefaultAssay(exTregCluster) <- 'RNA'
# all of cluster 7
DoHeatmap(exTregCluster, features = topgenes$gene, angle = 90, size = 2) +
  theme(text=element_text(size = 7), legend.text = element_text(size = 10))

DefaultAssay(tregCluster) <- 'RNA'
# cluster 17 
DoHeatmap(tregCluster, features = topgenes$gene, angle = 90, size = 2) + 
  theme(text=element_text(size = 7), legend.text = element_text(size = 10))

DefaultAssay(tregClusterminus) <- 'RNA'
# cluster 17 minus bin 1
DoHeatmap(tregClusterminus, features = topgenes$gene,size = 2) +
  theme(text=element_text(size = 7), legend.text = element_text(size = 10))

##### Make heatmap for gene expression data##### 

gene_df <- as.data.frame(newSrt@assays$RNA@data) %>% t()
gene_meta <- newSrt@meta.data %>% select(binCluster)
gene_merge <- merge(gene_df, gene_meta, by=0)
gene_avg <- gene_merge %>% group_by(binCluster) %>% summarise_all(funs(mean)) %>% data.table()
gene_avg$Row.names <- NULL
gene_avg$binCluster <- gsub('CD4T','C',gene_avg$binCluster)
gene_avg <- gene_avg %>% t()
gene_avg <- gene_avg %>% row_to_names(row_number = 1)

write.csv(norm_avg,'AntoineGSEA/AntoineBinClusterMatrixGenes.csv')

geneScaled <- CreateSeuratObject(gene_avg)
genedata <- as.data.frame(t(as.matrix(GetAssayData(geneScaled))))
scale_limits <- c(-2,2)
scaled_gene <- data.frame(lapply(genedata,function(x)rescale(x,to=scale_limits)))
rownames(scaled_gene) <-  rownames(genedata)
scaled_gene <- t(scaled_gene)
geneScaled@assays$RNA@scale.data <- scaled_gene
geneScaled = as.matrix(geneScaled@assays$RNA@scale.data)

pal_cols <- colorRampPalette(c("blue","yellow","red"))

geneScaled <- as.data.frame(geneScaled)

clusters <- c('C_7_Bin_1','C_7_Bin_2', 
              'C_7_Bin_3','C_7_Bin_4'
              ,'C_7_Bin_5','C_7_Bin_6', 
              'C_7_Bin_7', 'C_7_Bin_8', 
              'C_7_Bin_9', 'C_17_Bin_1', 'C_17_Bin_2',
              'C_17_Bin_3',  'C_17_Bin_4', 
              'C_17_Bin_7', 'C_17_Bin_8', 'C_17_Bin_9', 
              'C_7_Bin_10')
scaled.7.17 <- geneScaled[clusters]

newnames_row1 <- lapply(
  rownames(scaled.7.17) ,
  function(x) bquote(bold(.(x))))

newnames_cols1 <- lapply(
  colnames(scaled.7.17),
  function(x) bquote(bold(.(x))))

pheatmap::pheatmap(scaled.7.17, color = pal_cols(100), cluster_rows = FALSE, 
                   main = "Scaled Heatmap for Cluster 7 and 17 Bulk Transcriptome, by Bin", 
                   fontsize = 10, fontsize_col = 10, fontsize_row = 2, fontsize_number=2, 
                   show_rownames = T, 
                   labels_row = as.expression(newnames_row1),
                   labels_col = as.expression(newnames_cols1))

Heatmap(scaled.7.17, column_order = co_1, cluster_rows=T, show_row_names = F,
        row_names_gp = gpar(fontsize=6))

scaledDown <- scaled.7.17[row.names(scaled.7.17) %in% topgenes$gene,]

newnames_row <- lapply(
  rownames(scaledDown) ,
  function(x) bquote(bold(.(x))))

newnames_cols <- lapply(
  colnames(scaledDown),
  function(x) bquote(bold(.(x))))

pheatmap::pheatmap(scaledDown, color = pal_cols(100), cluster_rows = FALSE, 
                   main = "Scaled Heatmap for Cluster 7 and 17 Bulk Transcriptome, by Bin", 
                   fontsize = 10, fontsize_col = 10, fontsize_row = 5.5, fontsize_number=2, 
                   labels_row = as.expression(newnames_row),
                   labels_col = as.expression(newnames_cols))

scaledDown <- as.matrix(scaledDown)

Heatmap(scaledDown, cluster_columns=T, cluster_rows=T, 
            row_names_gp = gpar(fontsize=6))
g <- column_order(y)
h <- row_order(y)

co_1 <- c('C_7_Bin_1','C_17_Bin_1','C_7_Bin_2',
          'C_17_Bin_2','C_7_Bin_3','C_17_Bin_3',
          'C_7_Bin_4', 'C_17_Bin_4', 'C_7_Bin_5',
          'C_7_Bin_6', 'C_7_Bin_7','C_17_Bin_7',
          'C_7_Bin_8', 'C_17_Bin_8','C_7_Bin_9', 
          'C_17_Bin_9', 'C_7_Bin_10')

Heatmap(scaledDown, column_order = co_1, cluster_rows=T, 
        row_names_gp = gpar(fontsize=6))

cbPalette3 = c('gray','gray','blue','lightgreen','red')
fiveClusters <- ggplot(as.data.frame(colData(reduced)),
                       aes(x = dc1,
                           y = dc2, color = Cluster)) + geom_point() + theme_classic() +
  xlab("DC1") + ylab("DC2") + scale_colour_manual(values=cbPalette3) + 
  ggtitle("Destiny DC Plot by Cluster")


plotBinLines(fiveClusters, 10)

###### Add in Cluster 2 ######

# Markers:

table(newSrt$Bins)

Cluster.2.7.17 <- subset(x = newSrt, idents = c('CD4T_7_Bin_1','CD4T_7_Bin_2', 
                                              'CD4T_7_Bin_3','CD4T_7_Bin_4'
                                              ,'CD4T_7_Bin_5','CD4T_7_Bin_6', 
                                              'CD4T_7_Bin_7', 'CD4T_7_Bin_8', 
                                              'CD4T_7_Bin_9', 'CD4T_17_Bin_1', 'CD4T_17_Bin_2',
                                              'CD4T_17_Bin_3',  'CD4T_17_Bin_4', 
                                              'CD4T_17_Bin_7', 'CD4T_17_Bin_8', 'CD4T_17_Bin_9', 
                                              'CD4T_7_Bin_10', 'CD4T_2_Bin_1',  
                                              'CD4T_2_Bin_10', 'CD4T_2_Bin_2', 'CD4T_2_Bin_3',
                                              'CD4T_2_Bin_4', 'CD4T_2_Bin_5',  
                                              'CD4T_2_Bin_6', 'CD4T_2_Bin_7', 'CD4T_2_Bin_8',
                                              'CD4T_2_Bin_9'))

# order the clusters accordingly
clusterLevels2 <- c('CD4T_7_Bin_1','CD4T_17_Bin_1','CD4T_2_Bin_1',
                   'CD4T_7_Bin_2','CD4T_17_Bin_2', 'CD4T_2_Bin_2',
                   'CD4T_7_Bin_3','CD4T_17_Bin_3','CD4T_2_Bin_3',
                   'CD4T_7_Bin_4', 'CD4T_17_Bin_4', 'CD4T_2_Bin_4',
                   'CD4T_7_Bin_5', 'CD4T_2_Bin_5',
                   'CD4T_7_Bin_6', 'CD4T_2_Bin_6',
                   'CD4T_7_Bin_7','CD4T_17_Bin_7', 'CD4T_2_Bin_7',
                   'CD4T_7_Bin_8', 'CD4T_17_Bin_8','CD4T_2_Bin_8',
                   'CD4T_7_Bin_9', 'CD4T_17_Bin_9', 'CD4T_2_Bin_9',
                   'CD4T_7_Bin_10', 'CD4T_2_Bin_10')

# change active idents and assay
Cluster.2.7.17@active.ident <- factor(x=Cluster.2.7.17@active.ident, levels = clusterLevels2)
Cluster.2.7.17@active.assay <- 'ADT'

# Cluster 2,7,17 marker single cell heatmap
DoHeatmap(Cluster.2.7.17, features = top10$gene, angle = 90, size = 1.5)

# Cluster 2,7,17 transcript single cell heatmap
Cluster.2.7.17@active.assay <- 'RNA'
DoHeatmap(Cluster.2.7.17, features = topgenes$gene, angle = 90, size = 2) +
  theme(text=element_text(size = 7), legend.text = element_text(size = 10))


Cluster.2.7 <- subset(x = Cluster.2.7.17, idents = c('CD4T_7_Bin_1','CD4T_7_Bin_2', 
                                                'CD4T_7_Bin_3','CD4T_7_Bin_4'
                                                ,'CD4T_7_Bin_5','CD4T_7_Bin_6', 
                                                'CD4T_7_Bin_7', 'CD4T_7_Bin_8', 
                                                'CD4T_7_Bin_9',
                                                'CD4T_7_Bin_10', 'CD4T_2_Bin_1',  
                                                'CD4T_2_Bin_10', 'CD4T_2_Bin_2', 'CD4T_2_Bin_3',
                                                'CD4T_2_Bin_4', 'CD4T_2_Bin_5',  
                                                'CD4T_2_Bin_6', 'CD4T_2_Bin_7', 'CD4T_2_Bin_8',
                                                'CD4T_2_Bin_9'))

# Cluster 2,7 transcript single cell heatmap
DoHeatmap(Cluster.2.7, features = topgenes$gene, angle = 90, size = 2) +
  theme(text=element_text(size = 7), legend.text = element_text(size = 10))

# Cluster 2,7 marker single cell heatmap
Cluster.2.7@active.assay <- 'ADT'
DoHeatmap(Cluster.2.7, features = top10$gene, angle = 90, size = 2)
# do heatmap if need be

# Cluster 2,17 heatmaps

Cluster.2.17 <- subset(x = Cluster.2.7.17, idents = c('CD4T_17_Bin_1', 'CD4T_17_Bin_2',
                                                'CD4T_17_Bin_3',  'CD4T_17_Bin_4', 
                                                'CD4T_17_Bin_7', 'CD4T_17_Bin_8', 'CD4T_17_Bin_9', 
                                                'CD4T_2_Bin_1',  
                                                'CD4T_2_Bin_10', 'CD4T_2_Bin_2', 'CD4T_2_Bin_3',
                                                'CD4T_2_Bin_4', 'CD4T_2_Bin_5',  
                                                'CD4T_2_Bin_6', 'CD4T_2_Bin_7', 'CD4T_2_Bin_8',
                                                'CD4T_2_Bin_9'))

# Cluster 2,17 marker single cell heatmap
Cluster.2.17@active.assay <- 'ADT'
DoHeatmap(Cluster.2.17, features = top10$gene, angle = 90, size = 2)

Cluster.2.17reduced <- subset(x = Cluster.2.17, idents = c( 'CD4T_17_Bin_2',
                                                      'CD4T_17_Bin_3',  'CD4T_17_Bin_4', 
                                                      'CD4T_17_Bin_7', 'CD4T_17_Bin_8', 'CD4T_17_Bin_9', 
                                                      'CD4T_2_Bin_10', 'CD4T_2_Bin_2', 'CD4T_2_Bin_3',
                                                      'CD4T_2_Bin_4', 'CD4T_2_Bin_5',  
                                                      'CD4T_2_Bin_6', 'CD4T_2_Bin_7', 'CD4T_2_Bin_8',
                                                      'CD4T_2_Bin_9'))
DoHeatmap(Cluster.2.17reduced, features = top10$gene, angle = 90, size = 2)

Cluster.2.17@active.assay <- 'RNA'
DoHeatmap(Cluster.2.17, features = topgenes$gene, angle = 90, size = 2) +
  theme(text=element_text(size = 7), legend.text = element_text(size = 10))

Cluster.2.17reduced@active.assay <- 'RNA'
DoHeatmap(Cluster.2.17reduced, features = topgenes$gene, angle = 90, size = 2) +
  theme(text=element_text(size = 7), legend.text = element_text(size = 10))

# bulk marker expression heatmap

clusters2 <- c('C_7_Bin_1','C_17_Bin_1','C_2_Bin_1',
               'C_7_Bin_2','C_17_Bin_2', 'C_2_Bin_2',
               'C_7_Bin_3','C_17_Bin_3','C_2_Bin_3',
               'C_7_Bin_4', 'C_17_Bin_4', 'C_2_Bin_4',
               'C_7_Bin_5', 'C_2_Bin_5',
               'C_7_Bin_6', 'C_2_Bin_6',
               'C_7_Bin_7','C_17_Bin_7', 'C_2_Bin_7',
               'C_7_Bin_8', 'C_17_Bin_8','C_2_Bin_8',
               'C_7_Bin_9', 'C_17_Bin_9', 'C_2_Bin_9',
               'C_7_Bin_10', 'C_2_Bin_10')
scaled.2.7.17 <- normScaled[clusters2]
scaled.2.7.17 <- as.matrix(scaled.2.7.17)
Heatmap(scaled.2.7.17, column_order = clusters2, cluster_rows=T)

# bulk transcript expression heatmap 

gene.2.7.17 <- geneScaled[clusters2]
gene.2.7.17 <- as.matrix(gene.2.7.17)
scaledGenes2.7.17 <- gene.2.7.17[row.names(gene.2.7.17) %in% topgenes$gene,]
Heatmap(scaledGenes2.7.17, column_order = clusters2, cluster_rows=T,
        row_names_gp = gpar(fontsize=6))




