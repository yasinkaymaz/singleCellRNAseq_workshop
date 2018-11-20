library(Seurat)
library(Matrix)
library(dplyr)
library(SIMLR)

pbmc.data <- Read10X(data.dir = "~/Downloads/scRNAseqWS/10Xdata/filtered_gene_bc_matrices/GRCh38/")

pbmc <- CreateSeuratObject(raw.data = pbmc.data, 
                           min.cells = 3,
                           min.genes = 200,
                           project = "10X_PBMC")
head(pbmc@meta.data$nUMI)

mito.genes <- grep(pattern="^MT-", x =rownames(x=pbmc@raw.data), value = TRUE)
mito.genes

percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes,]) / Matrix::colSums(pbmc@raw.data)
head(percent.mito)

pbmc <- AddMetaData(pbmc,
                    metadata = percent.mito,
                    col.name = "percent.mito")

VlnPlot(pbmc, features.plot = c("nGene", "nUMI", "percent.mito"))

ggplot(pbmc@meta.data, aes(nUMI, nGene))+geom_point()

pbmc <- FilterCells(pbmc,
                    subset.names = c("nGene", "percent.mito"),
                    low.thresholds = c(200, -Inf),
                    high.thresholds = c(2500, 0.05))

VlnPlot(pbmc, features.plot = c("nGene", "percent.mito"))

pbmc <- NormalizeData(pbmc,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)

head(as.matrix(pbmc@data)[1:10,1:10])

pbmc <- FindVariableGenes(pbmc,
                              mean.function = ExpMean,
                              dispersion.function = LogVMR,
                              do.plot = TRUE)
head(pbmc@hvg.info)

hv.genes <- head(rownames(pbmc@hvg.info),1000)
hv.genes

pbmc <- ScaleData(pbmc,
                  genes.use = hv.genes,
                  vars.to.regress = c("nUMI", "percent.mito")
                  )


pbmc <- RunPCA(pbmc,
               pc.genes = hv.genes)

PCAPlot(pbmc, dim.1 = 2, dim.2 =3)

PCElbowPlot(pbmc)


pbmc <- FindClusters(pbmc,
                     reduction.type = "pca",
                     dims.use = 1:5,
                     resolution = 0.6
                     )

head(pbmc@ident,10)

pbmc <- RunTSNE(pbmc,
                dims.use = 1:5,
                do.fast=TRUE)

TSNEPlot(pbmc, do.label = TRUE)

#Marker genes
pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = TRUE,
                               min.pct =0.25,
                               logfc.threshold = 0.25)

top5.markers <- pbmc.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

best.markers <- pbmc.markers %>% group_by(cluster) %>% top_n(1, avg_logFC)

FeaturePlot(pbmc,
            features.plot = best.markers$gene,
            cols.use = c("grey", "blue"),
            reduction.use = "tsne")

