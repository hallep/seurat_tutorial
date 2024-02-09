library(dplyr)
library(Seurat)
library(patchwork)

# ========== LOAD DATA ========== #

# Read10X() - reads matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv files
# pbmc.data = unique molecular identified (UMI) count matrix
pbmc.data <- Read10X(data.dir = "~/seurat_tutorial/data/filtered_gene_bc_matrices/hg19/")

# CreateSeuratObject() - creates Seurat object (contains data and analysis for single-cell data)
# pbmc[["RNA"]]$counts = count matrix
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# displays number of molecules of CD3D, TCL1A, and MS4A1 transcripts in the
#     first 30 cells (prints to terminal)
# pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# ========== QC METRICS ========== #

# PercentageFeatureSet() - creates QC feature
#     percentage of mtDNA
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# ========== VISUALIZE DATA ========== #

# head() - displays QC metrics for first 5 cells (prints to terminal)
# head(pbmc@meta.data, 5)

# VlnPlot() - visualizes QC metrics as a violin plot (saved as Rplots.pdf)
# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter() - visualizes QC metrics as a scatter plot (saved as Rplots.pdf)
# plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

# ========== FILTER & NORMALIZE DATA ========== #

# subset() - filters data
#     cells with 200 to 2500 features and <5% mt count
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

# NormalizeData() - normalizes individual cell feature expression by total expression,
#     scale.factor = X (default = 10000) - multiplies values by scale factor X
#     normalization.method = "X" (default = "LogNormalize") - transforms by method X
# pbmc[["RNA"]]$data = normalized data
pbmc <- NormalizeData(pbmc)

# ========== VARIABLE FEATURES ========== #

# FindVariableFeatures() - finds the features with the highest cell-to-cell variation
#     selection.method = "X" (typical = "vst") - chooses features via method X
#     nfeatures = X - returns X features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# top10 = the 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# VariableFeaturePlot: visualizes standard variance vs. average expression
#     of variable features (saved to Rplots.pdf)
# VariableFeaturePlot(pbmc)
# LabelPoints() - adds labels for top [points] most variable genes
# LabelPoints(plot = VariableFeaturePlot(pbmc), points = top10, repel = TRUE)

# ========== SCALING ========== #

# ScaleData() - shifts expression of each genes so that mean is 0
#               scales expression of each gene so that variance is 1
#     features = X (default: only variable features) - features to scale
#     vars.to.regress = "X" - remove heterogeneity associated with feature X
# pbmc[["RNA"]]$scale.data = scaled data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ========== LINEAR DIMENSIONAL REDUCTION ========== #

# RunPCA() - returns genes with the most positive  or negative loadings
#             representing genes with correlation or anticorrelation across cells
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# returns first [dims] principle components and
#     [nfeatures] genes with the most positive and negative loadings
# print(pbmc[["pca"]], dims = 1:5, nfeatures = 10)

# VizDimLoadings() - shows gene vs. PC_1 and gene vs. PC_2 as separate graphs
# VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# DimPlot() - shows genes as dots plotted on PC_2 vs. PC_1
# DimPlot(pbmc, reduction = "pca") + NoLegend()

# DimHeatmap() - represents PC_1 values for 500 cells by color
# DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# ========== DIMENSIONALITY ========== #
# ElbowPlot() - ranks principle components based on percentage of variance
#     that can be explained by each component
# ElbowPlot(pbmc)

# ========== CLUSTERING CELLS ========== #

# FindNeighbors() - K-nearest neighbor (KNN)
#     calculates KNN based on euclidean distance in PCA space
#     redefines edge weights between any two sells based on shared overlap
#       in their local neighborhoods (Jaccard similarity)
#     returns number of nodes and edegs
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# FindClusters() - clusters the cells to optimize the standard modularity function
#     resolution: determines granularity of downstream clustering
#         higher resolution, more clusters
#         0.4 - 1.2 is best for single-cell data of ~3k cells
#     returns number of communities
pbmc <- FindClusters(pbmc, resolution = 0.5)

# returns cluster IDs of first 5 cells
head(Idents(pbmc), 5)

# ========== NON-LINEAR DIMENSIONAL REDUCTION ========== #
# RunUMAP() - Uniform Manifold Approximation and Projection dimensional reduction
# need to install
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

# saves current Seurat object
saveRDS(pbmc, file = "pbmc_tutorial.rds")

# ========== DIFFERENTIALLY EXPRESSED FEATURES (BIOMARKER CLUSTERS) ========== #
# finds differential expressed genes (markers) for identity classes
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0,3))
head(cluster5.markers, n = 5)

# ========== ASSIGNING CELL TYPE IDENTITY ========== #
