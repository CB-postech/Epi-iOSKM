#######
# Script used for the QC process of single-cell RNA sequencing data
# from Total dorsal epidermis (first dataset)
#######


set.seed(1234)
# Libraries ---------------------------------------------------------------
library("SingleCellExperiment")
library("DropletUtils")
library("scater")
library("scran")



# Read and load data ------------------------------------------------------
back.oskm <- read10xCounts("./data/back_oskm/count/sample_filtered_feature_bc_matrix/")
back.ctl <- read10xCounts("./data/back_control/count/sample_filtered_feature_bc_matrix/")



# Remove CMO columns ------------------------------------------------------
back.oskm <- back.oskm[row.names(back.oskm[1:(nrow(back.oskm) - 12), ]), ]
back.ctl <- back.ctl[row.names(back.ctl[1:(nrow(back.ctl) - 12), ]), ]



# Set col and row names ---------------------------------------------------
# colnames, cell names ----
colnames(back.oskm) <- colData(back.oskm)$Barcode
colnames(back.ctl) <- colData(back.ctl)$Barcode
# rownames 1, unquify gene symbols
row.names(back.oskm) <-
  uniquifyFeatureNames(rowData(back.oskm)$ID, rowData(back.oskm)$Symbol)
row.names(back.ctl) <-
  uniquifyFeatureNames(rowData(back.ctl)$ID, rowData(back.ctl)$Symbol)
# rownames 2, replace underscores with dashes for duplicates
row.names(back.oskm) <- gsub("_", "-", row.names(back.oskm))
row.names(back.ctl) <- gsub("_", "-", row.names(back.ctl))



# Calculate QC metrics ----------------------------------------------------
# get mitochondrial genes
is.mito <- grep("^mt-", row.names(back.oskm))  # same for both objects
# calculate QC metrics
back.oskm <- addPerCellQC(
  back.oskm, subsets = list(MT = is.mito),
  percent.top = c(50, 100, 200, 500),
  threshold = 0
)
back.ctl <- addPerCellQC(
  back.ctl, subsets = list(MT = is.mito),
  percent.top = c(50, 100, 200, 500),
  threshold = 0
)
# log10 sum / detected
back.oskm$log10_sum <- log10(back.oskm$sum)
back.oskm$log10_detected <- log10(back.oskm$detected)
back.ctl$log10_sum <- log10(back.ctl$sum)
back.ctl$log10_detected <- log10(back.ctl$detected)



# Check distribution ------------------------------------------------------
### Filter thresholds ----
thr.umi <- 500
thr.mt <- 10  # percent

### Total UMI counts ----
hist(back.oskm$log10_sum, breaks = 100, xlim = c(2.5, 5.5),
     main = "Log10(Sum) Histogram of Back OSKM", xlab = NULL) %>%
  abline(v = log10(thr.umi), col = "red")
hist(back.ctl$log10_sum, breaks = 100, xlim = c(2.5, 5.5),
     main = "Log10(Sum) Histogram of Back Control", xlab = NULL) %>%
  abline(v = log10(thr.umi), col = "red")

### MT percent ----
hist(back.oskm$subsets_MT_percent, breaks = 100, xlim = c(0, 100),
     main = "Mito.percent Histogram of Back OSKM", xlab = NULL) %>%
  abline(v = thr.mt, col = "red")
hist(back.ctl$subsets_MT_percent, breaks = 100, xlim = c(0, 100),
     main = "Mito.percent Histogram of Back Control", xlab = NULL) %>%
  abline(v = thr.mt, col = "red")



# Filter results ----------------------------------------------------------
### back.oskm ----
keep.umi <- back.oskm$log10_sum > log10(thr.umi)
keep.mt <- back.oskm$subsets_MT_percent < thr.mt
back.oskm$keep <- keep.umi & keep.mt

### back.control ----
keep.umi <- back.ctl$log10_sum > log10(thr.umi)
keep.mt <- back.ctl$subsets_MT_percent < thr.mt
back.ctl$keep <- keep.umi & keep.mt

### Filter ----
back.oskm <- back.oskm[, back.oskm$keep]
back_ctl <- back.ctl[, back.ctl$keep]



# Merge objects -----------------------------------------------------------
### Put condition data in metadata ----
# full
back.oskm$condition <- "back_oskm"
back.ctl$condition <- "back_control"
# OSKM
back.oskm$oskm <- "oskm"
back.ctl$oskm <- "control"
# modify cell barcode according to condition
colnames(back.oskm) <- paste0(colnames(back.oskm), "-BO")
colnames(back.ctl) <- paste0(colnames(back.ctl), "-BC")

### Merge ----
back <- cbind(back.oskm, back.ctl)
# factorize metadata
back$condition <- factor(back$condition, levels = c("back_oskm", "back_control"))
back$oskm <- factor(back$oskm, levels = c("control", "oskm"))



# Normalization -----------------------------------------------------------
# quickCluster
clusters <- quickCluster(back)
# computeSumFactors
back <- computeSumFactors(back, clusters = clusters)
# log-normalize
norm <- logNormCounts(back, pseudo_count = 1)



# Feature selection -------------------------------------------------------
# modelGeneVar
var <- modelGeneVar(norm)
# plot
plot(var$mean, var$total,
     xlab = "mean_log-expression", ylab = "variance", ylim = c(0, 10))
curve(metadata(var)$trend(x), col = "blue", add = TRUE)
# HVGs
hvg <- getTopHVGs(var, fdr.threshold = 0.05)
length(hvg)



# Make Seurat object ------------------------------------------------------
back <- as.Seurat(norm, counts = "counts", data = "logcounts")
VariableFeatures(back) <- hvg



# Dimension reduction -----------------------------------------------------
# scale data and PCA
back <- back %>% 
  ScaleData() %>% 
  RunPCA()
# check dimension
ElbowPlot(back, ndims = 50) +
  geom_hline(yintercept = 3) + 
  geom_vline(xintercept = 10)
PCs <- 20



# Clustering --------------------------------------------------------------
back <- back %>% 
  RunPCA(npcs = PCs) %>% 
  FindNeighbors(dims = 1:PCs)
back <- FindClusters(back, resolution = 0.8)



# Visualization -----------------------------------------------------------
# UMAP
back <- RunUMAP(back, dims = 1:PCs, seed.use = 1234)



# Check -------------------------------------------------------------------
DimPlot(back, reduction = "umap", label = TRUE, label.size = 5) + 
  theme(legend.position = "none", axis.title = element_blank())



#####
