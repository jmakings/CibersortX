library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library("berryFunctions")
library(collections)
library(tibble)

data <- as.sparse(read.csv("/Users/jmakings/Desktop/CAVA_Data/deletedgenes.csv", sep=",", header = TRUE, row.names = 1))
data <- as.data.frame(data)

genes <- rownames(data)

kaplanBulkBig <- as.data.frame(read.csv("/Users/jmakings/Desktop/CAVA_Data/InputFiles/kaplanBulkBig.txt", sep="\t", header = TRUE, row.names = 1))
kaplanBulkBig[1:10,1:10]

kaplanGenes <- rownames(kaplanBulkBig)

kaplanBulkBig$Gene <- kaplanGenes
kaplanBulkBig <- kaplanBulkBig %>% select(Gene, everything())

new <- subset(kaplanBulkBig, kaplanBulkBig$Gene %in% genes)

write.table(new, "/Users/jmakings/Desktop/CAVA_Data/InputFiles/kaplanPBMCs.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

cava <- as.data.frame(read.csv("Desktop/CAVA_Data/CAVA_myeloidSplit_ss.txt", sep="\t", header = TRUE, row.names = 1))

# first use gsub so the row names arent unique anymore

# restrict this to just the 74 genes
cava$Gene <- rownames(cava)
cava <- cava %>% select(Gene, everything())

newCava <- subset(cava, cava$Gene %in% genes)

remove <- list()
for (i in 2:ncol(newCava)) {
  if (i %% 25000 == 0) {
    print(paste0("Samples converted: ",i))
  }
  if (sum(newCava[,i]) == 0) {
    remove <- append(remove, i)
  }
}
remove <- as.numeric(unlist(remove))

newCava2 = select(newCava, -remove)

cells <- colnames(newCava2)
cells <- gsub("\\..*", "", cells)
newCava2 <- rbind(cells, newCava2)

write.table(newCava2, "/Users/jmakings/Desktop/CAVA_Data/InputFiles/cava74SS.txt", sep ="\t", col.names = FALSE, row.names = FALSE)

newCava2T <- t(newCava2)

newCava2T <- as.data.frame(newCava2T)

# to get counts of various cell types
table(newCava2T$`1`)

# next: train/test split for validation

wholeBlood <- as.sparse(read.csv("/Users/jmakings/Desktop/CAVA_Data/GSE149938_umi_matrix.csv", sep=",", header = TRUE, row.names = 1))
wholeBlood <- as.data.frame(wholeBlood)
wholeBloodT <- as.data.frame(t(wholeBlood))

pbmc <- CreateSeuratObject(counts = wholeBloodT, project = "Xie_WholeBlood")
pbmc$CellType <- Idents(pbmc)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)

plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

pbmc$CellType <- Idents(pbmc)
table(Idents(pbmc))

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

myeloidTypes <- function(pbmc) {
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
    } else if (t == "hMDP" || t == "cMOP" || t == "preM" ) {
      largerTypes <- append(largerTypes, "MonocytePrecursor")
    } else if (t == "claM") {
      largerTypes <- append(largerTypes, "CM")
    } else if (t == "interM") {
    largerTypes <- append(largerTypes, "INT")
    } else if (t == "nonM") {
      largerTypes <- append(largerTypes, "NCM")
    } else {
      print(t)
    }
  }
  return(largerTypes)
}

counts_To_cpm <- function(df) {
  for (i in 1:ncol(df)) {
    if (i %% 25000 == 0) {
      print(paste0("Samples converted: ",i))
    }
    df[,i] <- (df[,i]/sum(as.vector(df[,i])) * 1000000)
  }
  return(df)
}

types <- myeloidTypes(pbmc)
names(types) <- colnames(pbmc)
types <- as.data.frame(types)
types <- t(types)
pbmc <- AddMetaData(pbmc, types, "MajorTypes")
table(pbmc$MajorTypes)

x <- pbmc@meta.data
select_type <- x %>% select(MajorTypes)

wbloodcpm <- counts_To_cpm(wholeBloodT)
bt <- t(types)
wbloodcpm <- rbind(bt, wbloodcpm)

wbloodcpmT <- t(wbloodcpm)
wbloodcpmT <- as.data.frame(wbloodcpmT)

names(wbloodcpmT)[names(wbloodcpmT) == '1'] <- "Type"

filt <- subset(wbloodcpmT, Type != "HSPC" & Type != "MonocytePrecursor" &
                       Type != "Neutrophil" & Type != "ery" & Type != "plasma" )

filt <- t(filt)
filt <- as.data.frame(filt)

filt %>% dplyr::group_by(Type) %>% summarise_if(is.numeric, funs(mean)) -> filtFin




