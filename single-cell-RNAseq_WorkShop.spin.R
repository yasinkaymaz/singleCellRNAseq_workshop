
library(Seurat)
library(Matrix)
library(dplyr)
library(SIMLR)

#Read the 10X scRNAseq data into a R
pbmc.data <- Read10X(data.dir = "~/Downloads/scRNAseqWS/10Xdata/filtered_gene_bc_matrices/GRCh38/")

head(pbmc.data)

#Create a seurat object
pbmc <- CreateSeuratObject(raw.data = pbmc.data, 
                           min.cells = 3, min.genes = 200, project = "10X_PBMC")

#Seurat S4 object
class(pbmc)
head(pbmc@raw.data)
head(as.matrix(pbmc@raw.data)[1:10,1:10])
head(pbmc@meta.data)

#Quality check:
# 1. mitochondrial transcript expression content:
mito.genes <- grep(pattern = "^MT-", 
                   x= rownames(pbmc@raw.data), value=TRUE)

mito.genes
percent.mito <- Matrix:colSums()

