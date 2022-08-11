library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library("berryFunctions")

LM22 <- as.data.frame(read.csv("/Users/jmakings/Desktop/CibersortX/InputFiles/LM22.txt", 
                               sep="\t", header = TRUE, row.names = 1))

colnames(LM22)
LM22_edit <- subset(LM22, select = -c(Plasma.cells, Macrophages.M0, 
                                      Macrophages.M1, Macrophages.M2, 
                                      Dendritic.cells.resting, 
                                      Dendritic.cells.activated, 
                                      Mast.cells.resting, 
                                      Mast.cells.activated, 
                                      Eosinophils,
                                      Neutrophils))

LM22_edit <- rownames_to_column(LM22_edit)
colnames(LM22_edit)[which(names(LM22_edit) == 'rowname')] <- "Gene symbol"

write.table(LM22_edit,"/Users/jmakings/Desktop/CibersortX/InputFiles/LM22_edit.txt", sep = "\t",row.names = FALSE, col.names = TRUE)

B_cells <- LM22_edit$B.cells.naive + LM22_edit$B.cells.memory
CD4_T_cells <- LM22_edit$T.cells.CD4.naive + LM22_edit$T.cells.CD4.memory.resting +
  LM22_edit$T.cells.CD4.memory.activated + LM22_edit$T.cells.follicular.helper + 
  LM22_edit$T.cells.regulatory..Tregs. + LM22_edit$T.cells.gamma.delta

NK_cells <- LM22_edit$NK.cells.resting + LM22_edit$NK.cells.activated

LM22_clustered <- data.frame(LM22_edit$`Gene symbol`, B_cells, CD4_T_cells, LM22_edit$T.cells.CD8, NK_cells, LM22_edit$Monocytes)

# change column names then upload! 
colnames(LM22_clustered) <- c('Gene Symbol', 'B cells', 'CD4 T', 'CD8 T', 'NK cells', 'Monocytes')
LM22_clustered[1:5,]

write.table(LM22_clustered,"/Users/jmakings/Desktop/CibersortX/InputFiles/LM22_clustered.txt", sep = "\t",row.names = FALSE, col.names = TRUE)



