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



#############---------------------###################

###### 2)  Second Part, Working with multiple scRNAseq datasets:

###### - Loading datasets and QC (with PCA)

###### - Finding common variable genes between multiple datasets

###### - Multi-dataset alignment with CCA and Clustering cells (tSNE)

###### - Before and after comparison (pre/post-Alignment cell clustering)

###### - Pros/Cons of this method
