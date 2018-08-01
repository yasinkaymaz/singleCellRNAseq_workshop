###### - Input data and analysis tool requirements (setup)
#Load the R packages
library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)

# - Datasets and pre-processing steps
# Here, we are assuming that users already have a 10X dataset and run through Cell Ranger(TM) pipeline to generate barcode count matrices:

# - Loading data to R and initial quality checks (Seurat)
#Downloaded "4k PBMCs from a Healthy Donor" data from https://support.10xgenomics.com/single-cell-gene-expression/datasets
#http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz
# 4k PBMCs from a Healthy Donor
# Single Cell Gene Expression Dataset by Cell Ranger 2.1.0
# Peripheral blood mononuclear cells (PBMCs) from a healthy donor (same donor as pbmc8k). PBMCs are primary cells with relatively small amounts of RNA (~1pg RNA/cell).
# 
# 4,340 cells detected
# Sequenced on Illumina Hiseq4000 with approximately 87,000 reads per cell
# 26bp read1 (16bp Chromium barcode and 10bp UMI), 98bp read2 (transcript), and 8bp I7 sample barcode
# Analysis run with --expect-cells=5000
# Published on November 8, 2017
# 
# This dataset is licensed under the Creative Commons Attribution license.

###### - Loading data and initial quality checks
pbmc.data <- Read10X(data.dir = "~/Downloads/filtered_gene_bc_matrices/GRCh38/")

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "10X_PBMC")
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")

###### - Filtration, Normalization, and Scaling the data
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

###### - Finding variable genes and cell subtypes (with tSNE)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(pbmc@var.genes)

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = FALSE, pcs.print = 1:5)

PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
PCElbowPlot(object = pbmc)


###### - Detecting marker genes of cell clusters
# Cluster cells using Shared Nearest Neighbor (SNN) method. First find k-nearest neighbors for every cell, then, construct a SNN graph.
# For more information about the algorithms, read Waltman and van Eck (2013) The European Physical Journal B.
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = TRUE, save.SNN = TRUE)

head(pbmc@meta.data)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = pbmc)

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
markers <- pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

FeaturePlot(object = pbmc, features.plot = markers$gene, cols.use = c("grey", "blue"), reduction.use = "tsne")



#########################################################################
######     Second Part: Working with multiple scRNAseq datasets #########
#########################################################################
######   Alignment workflow for the four mouse brain datasets   #########
######----------------------------------------------------------#########

library(Seurat)
library(Matrix)

###### - Loading datasets and QC (with PCA)
setwd("/Users/yasinkaymaz/Documents/Harvard_Informatics/Data_Explore/mouse/")

#Load the metadata of each study.
#Zeisel: Single-cell RNA-seq of mouse cerebral cortex
Zeisel.meta <- read.delim("workshopData/Zeisel.RunTable.txt",header = TRUE)
rownames(Zeisel.meta) <- Zeisel.meta$Run_s
#Tasic: Adult mouse cortical cell taxonomy by single cell transcriptomics
Tasic.meta <- read.delim("workshopData/Tasic.RunTable.txt",header = TRUE)
rownames(Tasic.meta) <- Tasic.meta$Run_s
#Romanov: Single-cell RNA-seq of mouse hypothalamus
Romanov.meta <- read.delim("workshopData/Romanov.RunTable.txt",header = TRUE)
rownames(Romanov.meta) <- Romanov.meta$Run_s
#Marques: RNA-seq analysis of single cells of the oligodendrocyte lineage from nine distinct regions of the anterior-posterior and dorsal-ventral axis of the mouse juvenile central nervous system
Marques.meta <- read.delim("workshopData/Marques.RunTable.txt",header = TRUE)
rownames(Marques.meta) <- Marques.meta$Run_s

#Load the expression data of each study.
Zeisel.data <- read.delim("workshopData/Zeisel.expression.txt", header = TRUE, row.names = 1)
dim(Zeisel.data)
Tasic.data <- read.delim("workshopData/Tasic.expression.txt", header = TRUE, row.names = 1)
dim(Tasic.data)
Romanov.data <- read.delim("workshopData/Romanov.expression.txt", header = TRUE, row.names = 1)
dim(Romanov.data)
Marques.data <- read.delim("workshopData/Marques.expression.txt", header = TRUE, row.names = 1)
dim(Marques.data)

# Convert to sparse matrices for efficiency
zeisel.data <- as(as.matrix(Zeisel.data), "dgCMatrix")
romanov.data <- as(as.matrix(Romanov.data), "dgCMatrix")
tasic.data <- as(as.matrix(Tasic.data), "dgCMatrix")
marques.data <- as(as.matrix(Marques.data), "dgCMatrix")

# Create and setup Seurat objects for each dataset
zeisel <- CreateSeuratObject(raw.data = zeisel.data)
zeisel <- FilterCells(zeisel, subset.names = "nGene", low.thresholds = 2500)
zeisel <- NormalizeData(zeisel)
zeisel <- FindVariableGenes(zeisel, do.plot = F, display.progress = F)
zeisel <- ScaleData(zeisel)
zeisel@meta.data$tech <- "zeisel"

tasic <- CreateSeuratObject(raw.data = tasic.data)
tasic <- FilterCells(tasic, subset.names = "nGene", low.thresholds = 2500)
tasic <- NormalizeData(tasic)
tasic <- FindVariableGenes(tasic, do.plot = F, display.progress = F)
tasic <- ScaleData(tasic)
tasic@meta.data$tech <- "tasic"

romanov <- CreateSeuratObject(raw.data = romanov.data)
romanov <- FilterCells(romanov, subset.names = "nGene", low.thresholds = 2500)
romanov <- NormalizeData(romanov)
romanov <- FindVariableGenes(romanov, do.plot = F, display.progress = F)
romanov <- ScaleData(romanov)
romanov@meta.data$tech <- "romanov"

marques <- CreateSeuratObject(raw.data = marques.data)
marques <- FilterCells(marques, subset.names = "nGene", low.thresholds = 2500)
marques <- NormalizeData(marques)
marques <- FindVariableGenes(marques, do.plot = F, display.progress = F)
marques <- ScaleData(marques)
marques@meta.data$tech <- "marques"

#Check the PCA and tSNE prior to alignment:
batches <- rbind(Zeisel.meta[,c("Owner","SRA_Study_s")],
                 Tasic.meta[,c("Owner","SRA_Study_s")],
                 Romanov.meta[,c("Owner","SRA_Study_s")],
                 Marques.meta[,c("Owner","SRA_Study_s")])
combined.data <- cbind(Zeisel.data, Tasic.data, Romanov.data, Marques.data)
dim(combined.data)
fourDataset <- CreateSeuratObject(raw.data = combined.data, project = "4dataset.Pre")

fourDataset <- AddMetaData(object = fourDataset, metadata = batches, col.name = c("Owner","SRA_Study_s"))
fourDataset <- FilterCells(fourDataset, subset.names = "nGene", low.thresholds = 2500)
fourDataset <- NormalizeData(fourDataset)
fourDataset <- FindVariableGenes(fourDataset, do.plot = F, display.progress = F)
fourDataset <- ScaleData(fourDataset)

fourDataset <- RunPCA(object = fourDataset, pc.genes = fourDataset@var.genes, pcs.compute = 10, do.print = FALSE)
PCAPlot(fourDataset, pt.size=1, group.by ="Owner", dim.1 = 1, dim.2 = 2)

fourDataset <- RunTSNE(fourDataset, reduction.use = "pca", dims.use = 1:10)
TSNEPlot(fourDataset, do.label = T, group.by ="Owner")

###### - Finding common variable genes between multiple datasets
# Determine genes to use for CCA, must be highly variable in at least 2 datasets
ob.list <- list(zeisel, romanov, tasic, marques, E18Neurons)
genes.use <- c()
for (i in 1:length(ob.list)) {
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}


###### - Multi-dataset alignment with CCA and Clustering cells (tSNE)
# Run multi-set CCA
mouseBrain.integrated <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 15)

# CC Selection
MetageneBicorPlot(mouseBrain.integrated, grouping.var = "tech", dims.eval = 1:15)

# Run rare non-overlapping filtering
mouseBrain.integrated <- CalcVarExpRatio(object = mouseBrain.integrated, reduction.type = "pca",
                                         grouping.var = "tech", dims.use = 1:10)
mouseBrain.integrated <- SubsetData(mouseBrain.integrated, subset.name = "var.ratio.pca",
                                    accept.low = 0.5)

# Alignment
mouseBrain.integrated <- AlignSubspace(mouseBrain.integrated,
                                       reduction.type = "cca",
                                       grouping.var = "tech",
                                       dims.align = 1:10)

###### - Before and after comparison (pre/post-Alignment cell clustering)
# t-SNE and Clustering
mouseBrain.integrated <- FindClusters(mouseBrain.integrated, reduction.type = "cca.aligned",
                                      dims.use = 1:10, save.SNN = T, resolution = 0.4)
mouseBrain.integrated <- RunTSNE(mouseBrain.integrated,
                                 reduction.use = "cca.aligned",
                                 dims.use = 1:10)
# Visualization
TSNEPlot(mouseBrain.integrated, do.label = T)
TSNEPlot(mouseBrain.integrated, do.label = T,group.by ="tech")

###### - Pros/Cons of this method
