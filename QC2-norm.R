#######
# Script used for the QC process of single-cell RNA-sequencing
# from mCherry-sorted epidermis (second dataset)
#######


set.seed(1234)
# Libraries ---------------------------------------------------------------
library("SingleCellExperiment")
library("scater")
library("scran")



# Read and load data ------------------------------------------------------
# read raw as SCE
data.raw <- Read10X("./data/240123_second/cellranger/multi/raw_feature_bc_matrix/")
# convert to SingleCellExperiment
sce.raw <- SingleCellExperiment(assays = list(counts = data.raw$`Gene Expression`))
# add CMO counts to colData
colData(sce.raw) <- cbind(colData(sce.raw), t(data.raw$`Multiplexing Capture`))
sce.raw %>% colData() %>% colnames()



# Filter by CMO -----------------------------------------------------------
# read CMO assignments
cmo <- read.csv("./data/240123_second/cellranger/multi/multiplexing_analysis/assignment_confidence_table.csv")
# change column names
colnames(cmo) <- gsub(pattern = "CMO", "assign.CMO", x = colnames(cmo))
# subset with cells with info in CMO assignment file
table(colnames(sce.raw) %in% cmo$Barcode)  # T: 20886 cells, F: 2623810
sce.raw <- sce.raw[, colnames(sce.raw) %in% cmo$Barcode]
# add CMO assignments to colData
row.names(cmo) <- cmo$Barcode
colData(sce.raw) <- cbind(colData(sce.raw),
                          cmo[colnames(sce.raw),
                              c("assign.CMO301", "assign.CMO302", "assign.CMO303",
                                "Multiplet", "Blank",
                                "Assignment", "Assignment_Probability")])



# Calculate QC metrics ----------------------------------------------------
# get mitochondrial genes
is.mito <- grep("^mt-", row.names(sce.raw))
# calculate QC metrics
sce.raw <- addPerCellQC(
  sce.raw, subsets = list(MT = is.mito),
  percent.top = c(50, 100, 200, 500),
  threshold = 0
)
# log10 sum, detected
sce.raw$log10_sum <- log10(sce.raw$sum)
sce.raw$log10_detected <- log10(sce.raw$detected)



# Check distribution ------------------------------------------------------
### filter thresholds ----
thr.umi <- 500
thr.mt <- 20  # percent

### Total UMI counts ----
hist(sce.raw$log10_sum, breaks = 50, xlim = c(2.5, 5),
     main = "Log10(counts)", xlab = NULL,
     cex.main = 2, cex.lab = 1.5, cex.axis = 1.5) %>% 
  abline(v = log10(thr.umi), col = "red")

### MT percent ----
hist(sce.raw$subsets_MT_percent,
     breaks = 50, xlim = c(0, 100),
     main = "Mito%", xlab = NULL,
     cex.main = 2, cex.lab = 1.5, cex.axis = 1.5) %>%
  abline(v = thr.mt, col = "red")



# Filter results ----------------------------------------------------------
keep.umi <- sce.raw$log10_sum > log10(thr.umi)
keep.mt <- sce.raw$subsets_MT_percent < thr.mt
sce.raw$keep <- keep.umi & keep.mt
# filter
sce.raw <- sce.raw[, keep]



# Subset by CMO Assignment ------------------------------------------------
# CMO301: control, CMO302: mCherry-pos, CMO303: mCherry-neg
sce.cmo <- sce.raw[, sce.raw$Assignment %in% c("CMO301", "CMO302", "CMO303")]



# Normalization -----------------------------------------------------------
# quickCluster
clusters <- quickCluster(sce.cmo)
# computeCumFactors
sce.cmo <- computeSumFactors(sce.cmo, clusters = clusters)
# log-normalize
norm <- logNormCounts(sce.cmo, pseudo_count = 1)



# Feature selection -------------------------------------------------------
# modelGeneVar
var <- modelGeneVar(norm)
# plot
plot(var$mean, var$total,
     xlab = "mean_log-expression", ylab = "variance", ylim = c(0, 10))
curve(metadata(var)$trend(x), col = "blue", add = TRUE)
# HVG
hvg <- getTopHVGs(var, fdr.threshold = 0.001)
hvg %>% length()  # 2037



# Make Seurat object ------------------------------------------------------
skin <- as.Seurat(norm, counts = "counts", data = "logcounts")
VariableFeatures(skin) <- hvg



# Dimension reduction -----------------------------------------------------
# scale data and PCA
skin <- skin %>% 
  ScaleData() %>% 
  RunPCA()
# check dimension
ElbowPlot(skin, ndims = 50) +
  geom_hline(yintercept = 1:3, col = "red") + 
  geom_vline(xintercept = c(10, 15), col = "red")
PCs <- 15



# Clustering --------------------------------------------------------------
skin <- skin %>% 
  RunPCA(npcs = PCs) %>% 
  FindNeighbors(dims = 1:PCs)
skin <- FindClusters(skin, resolution = 0.8)
skin$cluster.total <- skin$seurat_clusters
skin$seurat_clusters <- NULL
skin$originalexp_snn_res.0.8 <- NULL



# Visualization -----------------------------------------------------------
skin <- skin %>% 
  RunUMAP(
    reduction = "pca",
    dims = 1:PCs,
    reduction.name = "umap",
    reduction.key = "UMAP_"
  ) %>% 
  RunTSNE(
    reduction = "pca",
    dims = 1:PCs,
    reduction.name = "tsne",
    reduction.key = "TSNE_"
  )



# Check -------------------------------------------------------------------
DimPlot(skin, reduction = "umap", group.by = "cluster.total",
        label = T, label.size = 7) + NoLegend()
DimPlot(skin, reduction = "umap", split.by = "Assignment")
DimPlot(skin, reduction = "tsne", label = T, label.size = 7) + NoLegend()




# Refine Meta.data --------------------------------------------------------
### CMO count matrix ----
# Create another assay
skin[["CMO"]] <- CreateAssayObject(
  counts = t(skin@meta.data[, paste0("CMO", 301:312)]),
  min.cells = 0, min.features = 0)
# erase metadata
skin@meta.data[, paste0("CMO", 301:312)] <- NULL

### remove some meta.data ----
skin@meta.data %>% colnames()
skin$orig.ident <- NULL
skin@meta.data[, paste0("assign.CMO", 301:303)] <- NULL
skin$Multiplet <- NULL
skin$Blank <- NULL
skin@meta.data[, paste0("percent.top_", c("100", "200", "500"))] <- NULL
skin@meta.data[, paste0("subsets_MT_", c("sum", "detected"))] <- NULL
skin$total <- NULL
skin$keep <- NULL
skin$sizeFactor <- NULL
skin@meta.data %>% colnames()

### Add condition metadata ----
# CMO301: control, CMO302: mCherry-pos, CMO303: mCherry-neg
skin$condition <- as.character(skin$Assignment)
skin$condition[skin$condition == "CMO301"] <- "Control"
skin$condition[skin$condition == "CMO302"] <- "OSKM_mCherry-pos"
skin$condition[skin$condition == "CMO303"] <- "OSKM_mCherry-neg"
skin$condition %>% table()
skin$condition <- factor(
  skin$condition, levels = c("OSKM_mCherry-pos", "OSKM_mCherry-neg", "Control")
)

### Reorder metadata ----
skin@meta.data %>% colnames()
skin@meta.data <- skin@meta.data[, c(1, 2, 12, 13,
                                     5, 6, 9, 10, 8, 7,
                                     3, 4, 14,
                                     11)]
skin@meta.data %>% colnames()



#####