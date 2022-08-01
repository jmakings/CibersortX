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

cava <- as.sparse(read.csv("/Users/jmakings/Desktop/CibersortX/raw/raw.csv", sep=",", header = TRUE, row.names = 1))
cava <- CollapseSpeciesExpressionMatrix(cava)
cava <- as.data.frame(t(as.matrix(cava)))

# convert counts into cpm first
counts_To_cpm <- function(df) {
  for (i in 1:ncol(df)) {
    if (i %% 5000 == 0) {
      print(paste0("Samples converted: ",i))
    }
    df[,i] <- (df[,i]/sum(as.vector(df[,i])) * 1000000)
  }
  return(df)
}

# save to new dataframe
cava2 <- counts_To_cpm(cava)
# all columns should now add up to 1 million 

# create seurat object from this
seuratcava <- CreateSeuratObject(counts = cava2)

# counts
counts <- seuratcava@assays$RNA@counts

# load in metadata (with myeloid splits) 
metaDF <- data.frame(read.csv("/Users/jmakings/Desktop/CibersortX/raw/final_metadata_myeloid.csv", sep=",", header = TRUE, row.names = 1))
seuratcava@meta.data <- metaDF

# restart from here, except MERGE like in CAVA_pseduo_myeloid split instead
# of getting the counts from the metadata, this messes things up

# patient info into dataframe 
selectPat <- seuratcava@meta.data %>% select(Pat_ID)
cava3 <- merge(selectPat, t(cava2), by=0)

# cell type info into dataframe
selectCell <- seuratcava@meta.data %>% select(Cell_Type)
row.names(cava3) <- cava3$Row.names

cava4 <- merge(selectCell, cava3, by = 0)
row.names(cava4) <- cava4$Row.names
cava4 <- select(cava4, -c(Row.names, Row.names.y))

cellCounts <- seuratcava@meta.data %>% group_by(seuratcava@meta.data$Pat_ID) %>% summarise(n())

mean(cellCounts$'n()')
# 2284 cells per patients on average

#filter out Rem and RemT types

cava5 <- cava4[cava4$Cell_Type != "Rem" & cava4$Cell_Type != "RemT",]

# Find the mean counts after excluding Rem and RemT 
cellCounts2 <- cava5 %>% group_by(Pat_ID) %>% summarise(n())
mean(cellCounts2$'n()')
mean(cellCounts$`n()`) - mean(cellCounts2$`n()`)
# average 1920 cells per patient, so 363 cells lost on average

ssCava <- subset(cava5, select = -c(Pat_ID))

ssCava <- as.data.frame(t(ssCava))

row.names(ssCava)[1] <- "GeneSymbol"

# remove columns with all zeros
ssCava <- ssCava[, colSums(ssCava != 0) > 0 ]

# remove rows with all zeros
ssCava <- ssCava[rowSums(ssCava != 0) > 0,]

# this is to create single cell reference matrix 
write.table(ssCava, "/Users/jmakings/Desktop/CibersortX/InputFiles/allgenesCAVA.txt", 
            sep = "\t", col.names = FALSE, row.names = TRUE)

# this section to create pseudobulk for each patient

cava6 <- subset(cava5, select = -c(Cell_Type))

# this will create pseudobulk RNA cpm for each of the 61 patients
cava6 %>% group_by(Pat_ID) %>% summarise_if(is.numeric, funs(sum)) -> cava6g
cava6g$Pat_ID <- sub("^", "P", cava6g$Pat_ID)
 
cava6gt <- as.data.frame(t(cava6g))
cava6gt %>% row_to_names(row_number = 1) -> cava7
row.names(cava6gt)[1] <- "Gene"

# this file converts the counts to cpm before creating pseudobulk
write.table(cava6gt, "/Users/jmakings/Desktop/CibersortX/InputFiles/cavaPatientsLarge.txt", 
            sep = "\t", col.names = FALSE, row.names = TRUE)

# now convert the pseudo

df_numeric <- function(df) {
  for (i in colnames(df)) {
    print(df$i)
    df$i <- as.numeric(df$i)
  }
  return(df)
}

cava7 <- df_numeric(cava7)

# To see the percentages of each cell type
trueVals <- tabyl(cava5, Pat_ID, Cell_Type) %>%
  adorn_percentages("row") %>% 
  adorn_pct_formatting(digits = 1)

# divide by 100 to get decimal instead of percentage
tabyl_to_decimal <- function(df) {
  trueVals %>% mutate_each(funs(as.numeric(gsub("%","", ., fixed= TRUE)))) -> remove
  dec <- remove/100
  return(dec)
}

# these are the true CAVA fractions for each sample and cell type
trueVals <- tabyl_to_decimal(trueVals)

# calculate RMSE between real and predicted

# this one with batch correction (B-mode)
preds <- as.data.frame(read.csv("/Users/jmakings/Desktop/CibersortX/fractionsFiles/CIBERSORTx_Job151_Results.csv", header = TRUE, row.names = 1))

# this one without batch correction
preds2 <- as.data.frame(read.csv("/Users/jmakings/Desktop/CibersortX/fractionsFiles/CIBERSORTx_Job150_Results.csv", header = TRUE, row.names = 1))
 
# RMSE with actual and predicted dataframes
mse_summary <- function(actual, predict) {
  types <- colnames(actual)
  types <- types[types != "Pat_ID"]
  results <- as.data.frame(matrix(nrow=2, ncol=length(types)))
  colnames(results) <- types
  rownames(results) <- c("MSE", "MAE")
  numPreds <- length(row.names(actual))
  numActual <- length(row.names(predict))
  if (numActual != numPreds) {
    print("Error: Number of predictions does not match number of true values")
    return(NULL)
  }
  count <- 1
  for (i in types) {
    act <- actual[i]
    pred <- predict[i]
    mse <- mse(act[1:numActual,], pred[1:numActual,])
    results[1,count] <- mse
    mae <- mae(act[1:numActual,], pred[1:numActual,])
    results[2,count] <- mae
    count <- count + 1
  }
  return(results)
}

# Function for per-patient MSE and MAE. MSE dataframe first, then MAE
mse_perPatient <- function(actual, predict) {
  types <- colnames(actual)
  types <- types[types != "Pat_ID"]
  sampNames <- row.names(predict)
  numPreds <- length(row.names(actual))
  numActual <- length(row.names(predict))
  if (numActual != numPreds) {
    print("Error: Number of predictions does not match number of true values")
    return(NULL)
  }
  mseDF <- as.data.frame(matrix(nrow=numPreds, ncol=length(types)))
  maeDF <- as.data.frame(matrix(nrow=numPreds, ncol=length(types)))
  colnames(mseDF) <- types
  colnames(maeDF) <- types
  rownames(mseDF) <- sampNames
  rownames(maeDF) <- sampNames
  
  count <- 1
  for (i in types) {
    act <- actual[i]
    pred <- predict[i]
    for (j in 1:numPreds) {
      mse <- mse(act[j,], pred[j,])
      mae <- mae(act[j,], pred[j,])
      mseDF[j,count] <- mse
      maeDF[j,count] <- mae
    }
    count <- count + 1
  }
  ret <- c("MSE",mseDF, "MAE", maeDF)
  return(ret)
  #return(maeDF)
}

# this one only returns MAE, better for evaluating each patient individually
mae_perPatient <- function(actual, predict) {
  types <- colnames(actual)
  types <- types[types != "Pat_ID"]
  sampNames <- row.names(predict)
  numPreds <- length(row.names(actual))
  numActual <- length(row.names(predict))
  if (numActual != numPreds) {
    print("Error: Number of predictions does not match number of true values")
    return(NULL)
  }
  maeDF <- as.data.frame(matrix(nrow=numPreds, ncol=length(types)))
  colnames(maeDF) <- types
  rownames(maeDF) <- sampNames
  
  count <- 1
  for (i in types) {
    act <- actual[i]
    pred <- predict[i]
    for (j in 1:numPreds) {
      mae <- mae(act[j,], pred[j,])
      maeDF[j,count] <- mae
    }
    count <- count + 1
  }
  return(maeDF)
}

# Results for CAVA on CAVA, all genes used, with and without batch correction: 

# This is WITH batch correction
mse_summary(trueVals, preds)

mse_perPatient(trueVals,preds)
mae <- mae_perPatient(trueVals,preds)

# This is WITHOUT batch correction
mse_summary(trueVals, preds2)

mse_perPatient(trueVals,preds2)
mae2 <- mae_perPatient(trueVals,preds2)

# creates heatmap of MAE matrix, params: MAE dataframe, heatmap low color, heatmap high color, and title for heatmap
maeHeat <- function(mae, colorlow, colorhigh, title) {
  if ("Pat_ID" %in% colnames(mae)) {
    break
  } else {
    mae <- rownames_to_column(mae, "Pat_ID")
  }
  mae <- melt(mae, id.vars = "Pat_ID")
  colnames(mae) <- c("Pat_ID", "Cell_Type", "MAE")
  plot <- ggplot(mae, aes(x=Cell_Type, y = Pat_ID, fill=MAE)) + geom_tile() + 
    #scale_fill_gradient2(low = colorlow, high = colorhigh, limit = c(0,0.35)) + 
    ggtitle(title) + scale_fill_gradientn(colors = c("white", "yellow", "orange", "red"),values = c(0 ,0.15,0.3,1), limit = c(0,0.4))
    #scale_color_stepsn(colors = terrain.colors(10))
  return(plot)
}

heat1 <- maeHeat(mae, "White", "Red", "CIBERSORTx Mean Absolute Error\n on CAVA (all genes used\n in signature matrix)")

# then, test on 74 genes
# single cell file

genes74 <- as.sparse(read.csv("/Users/jmakings/Desktop/CibersortX/InputFiles/deletedgenes.csv", sep=",", header = TRUE, row.names = 1))
genes74 <- as.data.frame(genes74)

GeneList <- rownames(genes74)

cava74 <- subset(cava5, select = c("Cell_Type", "Pat_ID", GeneList))

cava74ss <- as.data.frame(t(subset(cava74, select = -c(Pat_ID))))
rownames(cava74ss)[1] <- "GeneSymbol"

# remove columns with all zeros
cava74ss <- cava74ss[, colSums(cava74ss != 0) > 0 ]

# remove rows with all zeros
cava74ss <- cava74ss[rowSums(cava74ss != 0) > 0,]

write.table(cava74ss, "/Users/jmakings/Desktop/CibersortX/InputFiles/cava74ssNew.txt", sep = "\t", col.names = FALSE, row.names = TRUE)

# pseudobulk file

cava74_2 <- subset(cava74, select = -c(Cell_Type))

# this will create pseudobulk RNA cpm for each of the 61 patients
cava74_2 %>% group_by(Pat_ID) %>% summarise_if(is.numeric, funs(sum)) -> cava74g
cava6g$Pat_ID <- sub("^", "P", cava6g$Pat_ID)

cava74gt <- as.data.frame(t(cava74g))
row.names(cava74gt)[1] <- "Gene"

write.table(cava74gt, "/Users/jmakings/Desktop/CibersortX/InputFiles/cava74pseudoNoNormal.txt", sep = "\t", col.names = FALSE, row.names = TRUE)


cava74cpm <- counts_To_cpm(row_to_names(cava74gt, row_number = 1))
cava74cpm <- rbind(colnames(cava74cpm), cava74cpm)
rownames(cava74cpm)[1] <- "Gene"

write.table(cava74cpm, "/Users/jmakings/Desktop/CibersortX/InputFiles/cava74pseudoNormal.txt", sep = "\t", col.names = FALSE, row.names = TRUE)

preds3 <- as.data.frame(read.csv("/Users/jmakings/Desktop/CibersortX/fractionsFiles/CIBERSORTx_Job153_Results.csv", header = TRUE, row.names = 1))

preds4 <- as.data.frame(read.csv("/Users/jmakings/Desktop/CibersortX/fractionsFiles/CIBERSORTx_Job154_Results.csv", header = TRUE, row.names = 1))
mse_myeloid(trueVals,preds3)

mse_summary(trueVals,preds4)
# normalizing the bulk data to counts per million has no effect on the model

mae74 <- mae_perPatient(trueVals,preds4)

heat2 <- maeHeat(mae74,"white","red", "CIBERSORTx Mean Absolute\n Error on CAVA (74 genes\n used in signature matrix)")


# up next: test with just myeloid? 

cavaMyeloid <- cava5 
cavaMyeloid$Cell_Type[cavaMyeloid$Cell_Type == "NCM" | cavaMyeloid$Cell_Type == "INT" | cavaMyeloid$Cell_Type == "CM"] <- "Myeloid"

sscavaM <- subset(cavaMyeloid, select = -c(Pat_ID))

sscavaM <- as.data.frame(t(sscavaM))

# remove columns with all zeros
sscavaM <- sscavaM[, colSums(sscavaM != 0) > 0 ]

# remove rows with all zeros
sscavaM <- sscavaM[rowSums(sscavaM != 0) > 0,]

row.names(sscavaM)[1] <- "GeneSymbol"

# this is to create single cell reference matrix 
write.table(sscavaM, "/Users/jmakings/Desktop/CibersortX/InputFiles/allgenesCAVAmyeloid.txt", 
            sep = "\t", col.names = FALSE, row.names = TRUE)

# create myeloid Pseudobulk

cavamB <- subset(cavaMyeloid, select = -c(Cell_Type))

# this will create pseudobulk RNA cpm for each of the 61 patients
cavamB %>% group_by(Pat_ID) %>% summarise_if(is.numeric, funs(sum)) -> cavamGroup

cavamGroup <- as.data.frame(t(cavamGroup))
row.names(cavamGroup)[1] <- "Gene"

write.table(cavamGroup, "/Users/jmakings/Desktop/CibersortX/InputFiles/allgenesCAVAmyeloidpseudo.txt", 
            sep = "\t", col.names = FALSE, row.names = TRUE)


 # evaluation of CAVA matrix trained on full dataset but only given 74 genes during testing

preds5 <- as.data.frame(read.csv("/Users/jmakings/Desktop/CibersortX/fractionsFiles/CIBERSORTx_Job157_Results.csv", header = TRUE, row.names = 1))

mse_summary(trueVals, preds5)

mse_summary(trueVals, preds2)

heat3 <- maeHeat(mae_perPatient(trueVals,preds5), "white","red", "CIBERSORTx model trained with \nall CAVA genes, tested on\n CAVA with just 74 genes")
heat1 + heat2 + heat3

mae74 <- mae_perPatient(trueVals,preds5)

summary74 <- mse_summary(trueVals, preds)

# next: Try this out on JUST myeloid cells










