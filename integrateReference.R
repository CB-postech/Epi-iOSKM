#######
# Script used for integrating reference data with 
# mCherry-sorted dorsal epidermis single cell data (second dataset)
#######


set.seed(1234)
# Load data ---------------------------------------------------------------
daniel.ife <- readRDS("./rds/2024/0311/daniel_ife.rds")
skin2.ife <- readRDS("./rds/2024/0325/skin2_ife-oskm.rds")



# Prepare Data for Integration --------------------------------------------
### Cell names ----
colnames(daniel.ife) <- paste0("daniel_", colnames(daniel.ife))
colnames(skin2.ife) <- paste0("oskm_", colnames(skin2.ife))

### Data Source meta.data ----
daniel.ife$data <- "daniel"
skin2.ife$data <- "oskm"

### Assay ----
names(daniel.ife@assays) <- "RNA"
colnames(daniel.ife@meta.data)[2:3] <- paste0(c("nCount", "nFeature"), "_RNA")
colnames(daniel.ife@meta.data)
DefaultAssay(daniel.ife) <- "RNA"

### data.list ----
data.list <- list(daniel.ife, skin2.ife)



# Select Features ---------------------------------------------------------
features <- SelectIntegrationFeatures(
  object.list = data.list,
  nfeatures = 2000
)
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})



# Integration -------------------------------------------------------------
# get anchors
data.anchors <- FindIntegrationAnchors(
  object.list = data.list,
  anchor.features = features,
  reduction = "rpca",
)
# integrate
skin <- IntegrateData(anchorset = data.anchors)



# Processing --------------------------------------------------------------
### Dimension reduction ----
skin <- skin %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50)
ElbowPlot(skin, ndims = 50) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  geom_hline(yintercept = c(2, 4, 6), color = "red") +
  geom_vline(xintercept = c(5, 10, 15), color = "green")
PCs <- 10

### Visualization ----
skin <- skin %>% 
  RunPCA(npcs = PCs) %>%
  RunUMAP(
    dims = 1:PCs,
    reduction.name = "umap.integrated",
    reduction.key = "integratedUMAP_"
  ) %>% 
  FindNeighbors(reduction = "pca", dims = 1:PCs)
skin <- FindClusters(skin, resolution = 0.5)
skin$cluster.integrated <- skin$seurat_clusters
skin$seurat_clusters <- NULL
skin$integrated_snn_res.0.5 <- NULL

### Check ----
DimPlot(skin, reduction = "umap.integrated", 
        group.by = "cluster.integrated", pt.size = 0.7)
DimPlot(skin, reduction = "umap.integrated",
        group.by = "cluster.integrated", split.by = "data",
        pt.size = 0.7)



# Process RNA Data --------------------------------------------------------
DefaultAssay(skin) <- "RNA"
skin <- skin %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()



#####