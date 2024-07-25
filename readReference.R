#######
# Script used for processing reference single cell data
#######


set.seed(1234)
# Libraries ---------------------------------------------------------------
library("SingleCellExperiment")
library("scater")
library("scran")



# Read Data ---------------------------------------------------------------
# Read10X
uw1 <- Read10X("./data/public/GSE142471/uW1_GSM4230076/")
uw2 <- Read10X("./data/public/GSE142471/uW2_GSM4230077/")
wo1 <- Read10X("./data/public/GSE142471/W1_GSM4230078/")
wo2 <- Read10X("./data/public/GSE142471/W2_GSM4230079/")
wo3 <- Read10X("./data/public/GSE142471/W3_GSM4230080/")
# Create seurat
uw1 <- CreateSeuratObject(counts = uw1, project = "skin_public")
uw2 <- CreateSeuratObject(counts = uw2, project = "skin_public")
wo1 <- CreateSeuratObject(counts = wo1, project = "skin_public")
wo2 <- CreateSeuratObject(counts = wo2, project = "skin_public")
wo3 <- CreateSeuratObject(counts = wo3, project = "skin_public")
# Metadata
meta <- read.csv(
  file = "./data/public/GSE142471/meta_Total_Cells_integration_UW_WO_CellRep2020_FigS2E.txt",
  sep = "\t")



# Apply Metadata ----------------------------------------------------------
## Refine barcodes ----
# Set row.names to barcodes
row.names(meta) <- meta$Row


## Define sample relation ----
# bs: unwounded / sw: wounded
### Subset ----
meta.bs1 <- meta[meta$Run == "bs_1", ]
meta.bs2 <- meta[meta$Run == "bs_2", ]
meta.sw1 <- meta[meta$Run == "sw_1", ]
meta.sw2 <- meta[meta$Run == "sw_2", ]
meta.sw3 <- meta[meta$Run == "sw_3", ]

### Change barcodes ----
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4230076
# The order of aggregation is
# uw1, wo1, uw2, wo2, wo3
# thus, the number at the end of the cell barcodes are respectively 1~5
# metadata follows this, cell barcodes in seurat don't
Cells(uw1)  # no need to change
colnames(wo1) <- gsub("1", "2", colnames(wo1)); Cells(wo1)
colnames(uw2) <- gsub("1", "3", colnames(uw2)); Cells(uw2)
colnames(wo2) <- gsub("1", "4", colnames(wo2)); Cells(wo2)
colnames(wo3) <- gsub("1", "5", colnames(wo3)); Cells(wo3)

### Check overlap ----
# TRUE / FALSE
table(Cells(uw1) %in% row.names(meta.bs1))  # 5260 / 112
table(Cells(uw2) %in% row.names(meta.bs2))  # 5355 / 162
table(Cells(wo1) %in% row.names(meta.sw1))  # 7501 / 77
table(Cells(wo2) %in% row.names(meta.sw2))  # 4557 / 141
table(Cells(wo3) %in% row.names(meta.sw3))  # 4106 / 46



# Subset Based on Metadata ------------------------------------------------
uw1 <- uw1[, colnames(uw1) %in% row.names(meta.bs1)]
uw2 <- uw2[, colnames(uw2) %in% row.names(meta.bs2)]
wo1 <- wo1[, colnames(wo1) %in% row.names(meta.sw1)]
wo2 <- wo2[, colnames(wo2) %in% row.names(meta.sw2)]
wo3 <- wo3[, colnames(wo3) %in% row.names(meta.sw3)]



# Add Metadata ------------------------------------------------------------
### Add run metadata ----
uw1$run <- "UW1"
uw2$run <- "UW2"
wo1$run <- "WO1"
wo2$run <- "WO2"
wo3$run <- "WO3"

### Add condition metadata ----
uw1$condition <- "un-wounded"
uw2$condition <- "un-wounded"
wo1$condition <- "wounded"
wo2$condition <- "wounded"
wo3$condition <- "wounded"

### Add cell type metadata ----
uw1$cellType <- meta.bs1[colnames(uw1), "final_labels"]
uw2$cellType <- meta.bs2[colnames(uw2), "final_labels"]
wo1$cellType <- meta.sw1[colnames(wo1), "final_labels"]
wo2$cellType <- meta.sw2[colnames(wo2), "final_labels"]
wo3$cellType <- meta.sw3[colnames(wo3), "final_labels"]



# Merge -------------------------------------------------------------------
## Merge ----
sce.daniel <- cbind(
  SingleCellExperiment(assays = list(counts = uw1@assays$RNA$counts)),
  SingleCellExperiment(assays = list(counts = uw2@assays$RNA$counts)),
  SingleCellExperiment(assays = list(counts = wo1@assays$RNA$counts)),
  SingleCellExperiment(assays = list(counts = wo2@assays$RNA$counts)),
  SingleCellExperiment(assays = list(counts = wo3@assays$RNA$counts))
)


## Metadata ----
# run
sce.daniel$run <- c(
  rep("UW1", times = ncol(uw1)),
  rep("UW2", times = ncol(uw2)),
  rep("WO1", times = ncol(wo1)),
  rep("WO2", times = ncol(wo2)),
  rep("WO3", times = ncol(wo3))
)
sce.daniel$run %>% is.na() %>% sum()
# condition
sce.daniel$condition <- c(
  rep("un-wounded", times = (ncol(uw1) + ncol(uw2))),
  rep("wounded", times = (ncol(wo1) + ncol(wo2) + ncol(wo3)))
)
sce.daniel$condition %>% is.na() %>% sum()
# cellType
sce.daniel$cellType <- NA
names(sce.daniel$cellType) <- colnames(sce.daniel)
sce.daniel$cellType[colnames(uw1)] <- meta.bs1[colnames(uw1), "final_labels"]
sce.daniel$cellType[colnames(uw2)] <- meta.bs2[colnames(uw2), "final_labels"]
sce.daniel$cellType[colnames(wo1)] <- meta.sw1[colnames(wo1), "final_labels"]
sce.daniel$cellType[colnames(wo2)] <- meta.sw2[colnames(wo2), "final_labels"]
sce.daniel$cellType[colnames(wo3)] <- meta.sw3[colnames(wo3), "final_labels"]
sce.daniel$cellType %>% is.na() %>% sum()



# Preprocess --------------------------------------------------------------
## Feature selection ----
# quickCluster
clusters.daniel <- quickCluster(sce.daniel)
# computeSumFactors
sce.daniel <- computeSumFactors(sce.daniel, clusters = clusters.daniel)
# log-Normalize
norm.daniel <- logNormCounts(sce.daniel, pseudo_count = 1)


## Feature selection ----
# modelGeneVar
var.daniel <- modelGeneVar(norm.daniel)
# plot
plot(var.daniel$mean, var.daniel$total,
     xlab = "mean_log-expression", ylab = "variance", ylim = c(0, 10))
curve(metadata(var.daniel)$trend(x), col = "blue", add = TRUE)
# HVGs
hvg.daniel <- getTopHVGs(var.daniel, fdr.threshold = 0.05)
hvg.daniel %>% length()  # 936 / 1540



# Make Seurat -------------------------------------------------------------
daniel.merge <- as.Seurat(norm.daniel, counts = "counts", data = "logcounts")
VariableFeatures(daniel.merge) <- hvg.daniel



# Dimension reduction -----------------------------------------------------
# scale data and PCA
daniel.merge <- daniel.merge %>% 
  ScaleData() %>% 
  RunPCA()
# ElbowPlot
ElbowPlot(daniel.merge, ndims = 50) +
  geom_hline(yintercept = c(3, 5)) +
  geom_vline(xintercept = c(10, 15))
PCs <- 10



# Clustering --------------------------------------------------------------
daniel.merge <- daniel.merge %>% 
  RunPCA(npcs = PCs) %>% 
  FindNeighbors(dims = 1:PCs)
daniel.merge <- daniel.merge %>% FindClusters(resolution = 0.8)
# Reduction
daniel.merge <- daniel.merge %>% 
  RunUMAP(dims = 1:PCs) %>% 
  RunTSNE(dims = 1:PCs)



# Refine Metadata ---------------------------------------------------------
### Run ----
daniel.merge$run <- factor(daniel.merge$run,
                           levels = c("UW1", "UW2", "WO1", "WO2", "WO3"))

### Condition ----
daniel.merge$condition <- factor(daniel.merge$condition, levels = c("un-wounded", "wounded"))

### cellType.simple ----
daniel.merge$cellType.simple <- NA
daniel.merge$cellType.simple[grep("^Basal", daniel.merge$cellType, ignore.case = T)] <- "Basal"
daniel.merge$cellType.simple[daniel.merge$cellType == "Prolif. Basal"] <- "Prolif.Basal"
daniel.merge$cellType.simple[daniel.merge$cellType == "Spinous"] <- "Spinous"
daniel.merge$cellType.simple[daniel.merge$cellType == "HFSC"] <- "HFSC"
daniel.merge$cellType.simple[grep("^HF ", daniel.merge$cellType)] <- "HF"
daniel.merge$cellType.simple[grep("^Fibro", daniel.merge$cellType)] <- "Fibroblast"
daniel.merge$cellType.simple[grep("^Myo", daniel.merge$cellType)] <- "Myofibroblast"
daniel.merge$cellType.simple[grep("^Endo", daniel.merge$cellType)] <- "Endothelial"
daniel.merge$cellType.simple[grep("^Skel", daniel.merge$cellType)] <- "Skeletal muscle"
daniel.merge$cellType.simple[grep("^T cell ", daniel.merge$cellType)] <- "T cell"
daniel.merge$cellType.simple[grep("^Macro", daniel.merge$cellType)] <- "Macrophage"
daniel.merge$cellType.simple[grep("^Lang", daniel.merge$cellType)] <- "Langerhans"
daniel.merge$cellType.simple[grep("^Dendrit", daniel.merge$cellType)] <- "Dendritic cell"
daniel.merge$cellType.simple %>% is.na() %>% sum()
daniel.merge$cellType.simple <- factor(
  daniel.merge$cellType.simple,
  levels = c("Basal", "Prolif.Basal", "Spinous", "HFSC", "HF",
             "Fibroblast", "Myofibroblast", "Endothelial", "Skeletal muscle",
             "T cell", "Macrophage", "Langerhans", "Dendritic cell"))



# Plot Data ---------------------------------------------------------------
# cellType DimPlot
p <- DimPlot(daniel.merge, group.by = "cellType") +
  NoLegend() + NoAxes() +
  theme(plot.title = element_blank())
LabelClusters(plot = p, id = "cellType", repel = T,
              size = 4, fontface = "bold")
# condition DimPlot
DimPlot(daniel.merge, group.by = "condition") +
  NoLegend() + NoAxes() +
  theme(plot.title = element_blank())



# Subset IFE --------------------------------------------------------------
## Subset ----
daniel.merge.ife <- daniel.merge[, daniel.merge$cellType.simple %in% 
                                   c("Basal", "Prolif.Basal", "Spinous")]


## Make SCE ----
sce.daniel.ife <- as.SingleCellExperiment(daniel.merge.ife)


## Normalization ----
# quickCluster
cluster.daniel <- quickCluster(sce.daniel.ife)
# computeSumFactors
sce.daniel.ife <- computeSumFactors(sce.daniel.ife, clusters = cluster.daniel)
# log-Normalize
norm.daniel <- logNormCounts(sce.daniel.ife, pseudo_count = 1)


## Feature selection ----
# modelGeneVar
var.daniel <- modelGeneVar(norm.daniel)
# plot
plot(var.daniel$mean, var.daniel$total,
     xlab = "mean_log-expression", ylab = "variance", ylim = c(0, 10))
curve(metadata(var.daniel)$trend(x), col = "blue", add = TRUE)
# HVG
hvg.daniel <- getTopHVGs(var.daniel, fdr.threshold = 0.05)
hvg.daniel %>% length()  # 250


## Make Seurat ----
daniel.merge.ife <- as.Seurat(norm.daniel, counts = "counts", data = "logcounts")
VariableFeatures(daniel.merge.ife) <- hvg.daniel


## Dimension reduction ----
# remove existing information
daniel.merge.ife$ident <- NULL
daniel.merge.ife@reductions$PCA <- NULL
daniel.merge.ife@reductions$UMAP <- NULL
daniel.merge.ife@reductions$TSNE <- NULL
# scale data and PCA
daniel.merge.ife <- daniel.merge.ife %>% 
  ScaleData() %>% RunPCA()
# check dimension
ElbowPlot(daniel.merge.ife, ndims = 50) +
  geom_hline(yintercept = c(3, 5)) +
  geom_vline(xintercept = c(10, 15))
PCs <- 5


## Clustering ----
daniel.merge.ife <- daniel.merge.ife %>% 
  RunPCA(npcs = PCs) %>% 
  FindNeighbors(dims = 1:PCs)
daniel.merge.ife <- daniel.merge.ife %>% FindClusters(resolution = 0.8)
daniel.merge.ife$cluster.ife <- daniel.merge.ife$seurat_clusters
daniel.merge.ife$seurat_clusters <- NULL
daniel.merge.ife$originalexp_snn_res.0.8 <- NULL


## Visualization ----
daniel.merge.ife <- daniel.merge.ife %>% 
  RunUMAP(dims = 1:PCs) %>% 
  RunTSNE(dims = 1:PCs)


## Plot Seurat ----
DimPlot(daniel.merge.ife, label = T, label.size = 4, repel = T) + NoLegend()
DimPlot(daniel.merge.ife, split.by = "condition")
DimPlot(daniel.merge.ife, group.by = "cellType", split.by = "condition")



# Check  ------------------------------------------------------------------
## Plot by clusters and cellType ----
DimPlot(daniel.merge.ife,
             group.by = "cluster.ife"
             # group.by = "cellType"
        ) +
  theme_classic() + NoAxes() + NoLegend() + ggtitle("Total") +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  scale_x_continuous(limits = c(-12.5, 11), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-12, 10), expand = c(0, 0))


## Plot density ----
# get embedding information
plt.umap <- Embeddings(daniel.merge.ife, reduction = "umap") %>% as.data.frame()
plt.meta <- daniel.merge.ife@meta.data
plt.umap <- cbind(plt.umap, plt.meta); rm(plt.meta)
plt.umap %>% head()

### un-wounded ----
plt.umap %>% 
  filter(condition == "un-wounded") %>%
  ggplot() +
  # plot cells
  geom_point(
    aes(x = umap_1, y = umap_2, 
        # color = cluster.ife
        color = cellType
        ),
    size = 1.5, alpha = 0.2) +
  # plot density
  geom_density_2d(
    aes(x = umap_1, y = umap_2, alpha = after_stat(level)),
    bins = 10, linewidth = 0.6, color = "red") +
  theme_classic() + NoAxes() + NoLegend() + 
  ggtitle("un-wounded") +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  scale_x_continuous(limits = c(-12.5, 11), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-12, 10), expand = c(0, 0))

### wounded ----
plt.umap %>% 
  filter(condition == "wounded") %>%
  ggplot() +
  # plot cells
  geom_point(
    aes(x = umap_1, y = umap_2, 
        color = cluster.ife
        # color = cellType
    ),
    size = 1.5, alpha = 0.2) +
  # plot density
  geom_density_2d(
    aes(x = umap_1, y = umap_2, alpha = after_stat(level)),
    bins = 10, linewidth = 0.6, color = "red") +
  theme_classic() + NoAxes() + NoLegend() + 
  ggtitle("wounded") +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  scale_x_continuous(limits = c(-12.5, 11), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-12, 10), expand = c(0, 0))



#####