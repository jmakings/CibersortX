# CibersortX Whole Blood Deconvolution

## The goal of this project is to deconvolve the relative proportions of different Peripheral Blood Mononuclear Cells (PBMC) from whole blood bulk RNA-seq samples using the CIBERSORTx analytical framework for the Ley Lab at the La Jolla Institute for Immunology.

In this repo contains the various R scripts for pre-processing of single cell and bulk RNA-seq data for use with the CIBERSORTx framework. 

In the case of this project, CIBERSORTx takes single cell RNA-seq expression data and uses it to create a signature matrix of "barcode genes" that discriminate each cell subset of interest (T-cells, B-cells, Monocytes, etc.). 

Then, this signature matrix can be applied to bulk RNA-seq data to deconvolve the cell type proportions in bulk cell mixtures.

Links:
- [CIBERSORTx website](https://cibersortx.stanford.edu/)
- CIBERSORT foundational paper, [Newman et al. 2015](https://www.nature.com/articles/nmeth.3337#Sec10)
- CIBERSORTx paper, which includes support for single cell RNA-seq, [Newman et al. 2019](https://www.nature.com/articles/s41587-019-0114-2#MOESM1)