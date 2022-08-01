library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library("berryFunctions")

# Load in data and create Seurat object 
data <- as.sparse(read.csv("/Users/jmakings/Desktop/CAVA_Data/raw.csv", sep=",", header = TRUE, row.names = 1))
data <- CollapseSpeciesExpressionMatrix(data)
data <- t(as.matrix(data))

seurat <- CreateSeuratObject(counts = data, row.names = 0)

# get counts
counts <- seurat@assays$RNA@counts

# dataframe from counts
df <- as.data.frame(counts) %>% t()

# load metadata (with myeloid splits) and add to Seurat object
metaDF <- data.frame(read.csv("/Users/jmakings/Desktop/CAVA_Data/final_metadata_myeloid.csv", sep=",", header = TRUE, row.names = 1))
seurat@meta.data <- metaDF

## Section 1: Creating the single- cell reference file for creating signature matrix

dfT <- t(df)

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
cpm <- counts_To_cpm(dfT)

# check_cpm <- function(cpm) {
#   for (i in 1:ncol(cpm)) {
#     val <- sum(cpm[,i])
#     if (val != 1e+06) {
#       print(i)
#       print(val)
#       return(FALSE) 
#     }
#   }
#   return(TRUE)
# }

# dataframe cpm has now been converted to counts per million. Can check by 
# confirming the counts in each sample column equal 1 million
rand <- sample.int(ncol(cpm),1)
sum(cpm[,rand]) == 1e+06

x <- seurat@meta.data
select_meta <- x %>% select(Cell_Type)

# merging cell type data into the cpm array
cpmT <- as.data.frame(t(cpm))
ssFinal <- merge(select_meta, cpmT, by=0)

# remove Rem and RemT
ssFinal2 <- ssFinal[ssFinal$Cell_Type != "Rem" & ssFinal$Cell_Type != "RemT",]
ssFinalT <- t(ssFinal2)

# setting the cell type to the column name, removing unneccessary rows
ssFinalT <- ssFinalT[-c(1),]
ssFinalT <- as.data.frame(ssFinalT)
rownames(ssFinalT)[1] <- "GeneSymbol"

# Writing to file
write.table(ssFinalT,"/Users/jmakings/Desktop/CAVA_Data/CAVA_myeloidSplit_ss2.txt", sep = "\t",row.names = TRUE, col.names = FALSE)


###

###

## Section 2: Creating Pseudo-bulk mixture file


# summarise transcript counts by cell Type

# This one is count sums
ssFinal2 %>% dplyr::group_by(Cell_Type) %>% summarise_if(is.numeric, funs(sum)) -> df_summed

df_summedT <- t(df_summed)
df_summedT <- as.data.frame(df_summedT)

df_summedT <- tibble::rownames_to_column(df_summedT, "Gene")
df_summedT["Gene"][df_summedT["Gene"] == "Cell_Type"] <- "Gene"

write.table(df_summedT,"/Users/jmakings/Desktop/CAVA_Data/CAVA_pseudo_myeloidSum.txt", sep = "\t",row.names = FALSE, col.names = FALSE)

# This one is count means
ssFinal2 %>% dplyr::group_by(Cell_Type) %>% summarise_if(is.numeric, funs(mean)) -> df_summed2

df_summed2T <- t(df_summed2)
df_summed2T <- as.data.frame(df_summed2T)

df_summed2T <- tibble::rownames_to_column(df_summed2T, "Gene")
df_summed2T["Gene"][df_summed2T["Gene"] == "Cell_Type"] <- "Gene"


write.table(df_summed2T,"/Users/jmakings/Desktop/CAVA_Data/CAVA_pseudo_myeloidMean.txt", sep = "\t",row.names = FALSE, col.names = FALSE)







