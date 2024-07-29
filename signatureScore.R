#######
# Script used to perform module score calculation and visualization
#######


set.seed(1234)
# Libraries ---------------------------------------------------------------
library("msigdbr")



# Load data ---------------------------------------------------------------
skin2.ife <- readRDS("./rds/skin2_ife.rds")



# GO terms ----------------------------------------------------------------
### Get gene sets from msigdbr ----
gs.all <- msigdbr(species = "Mus musculus")
gs.gobp <- gs.all %>%
  filter(gs_subcat == "GO:BP") %>% 
  split(x = .$gene_symbol, f = .$gs_name)
# select terms
terms_ <- c(
  "GOBP_DEFENSE_RESPONSE",
  "GOBP_DEFENSE_RESPONSE_TO_OTHER_ORGANISM",
  "GOBP_BIOLOGICAL_PROCESS_INVOLVED_IN_INTERSPECIES_INTERACTION_BETWEEN_ORGANISMS",
  "GOBP_IMMUNE_RESPONSE",
  "GOBP_DEFENSE_RESPONSE_TO_SYMBIONT",
  "GOBP_INNATE_IMMUNE_RESPONSE",
  "GOBP_RESPONSE_TO_VIRUS",
  "GOBP_RESPONSE_TO_BACTERIUM",
  "GOBP_RESPONSE_TO_CYTOKINE",
  "GOBP_DEFENSE_RESPONSE_TO_BACTERIUM",
  "GOBP_REGULATION_OF_RESPONSE_TO_EXTERNAL_STIMULUS",
  "GOBP_RESPONSE_TO_INTERFERON_GAMMA",
  "GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS",
  "GOBP_RIBOSOME_BIOGENESIS",
  "GOBP_PROGRAMMED_CELL_DEATH",
  "GOBP_TISSUE_DEVELOPMENT")
gs.gobp <- gs.gobp[terms_]

### Calculate scores ----
skin2.ife <- AddModuleScore(skin2.ife, features = gs.gobp, name = "gobp.")

### score matrix ----
# score
skin2.ife$condition <- factor(
  skin2.ife$condition,
  levels = rev(c("Control", "OSKM_mCherry-neg", "OSKM_mCherry-pos")))
score.ife.b <- skin2.ife@meta.data[skin2.ife$cellType.epi == "IFE.B", ] %>% 
  select(cellType.epi, condition, starts_with("gobp.")) %>% 
  group_by(cellType.epi, condition) %>%
  summarise(across(starts_with("gobp."), mean)) %>%
  ungroup() %>% as.matrix()
score.ife.sb <- skin2.ife@meta.data[skin2.ife$cellType.epi == "IFE.SB", ] %>% 
  select(cellType.epi, condition, starts_with("gobp.")) %>% 
  group_by(cellType.epi, condition) %>%
  summarise(across(starts_with("gobp."), mean)) %>%
  ungroup() %>% as.matrix()
score.ife.bc <- skin2.ife@meta.data[skin2.ife$cellType.epi == "IFE.B.C", ] %>% 
  select(cellType.epi, condition, starts_with("gobp.")) %>% 
  group_by(cellType.epi, condition) %>%
  summarise(across(starts_with("gobp."), mean)) %>%
  ungroup() %>% as.matrix()
# metadata and matrix
meta.ife.b <- as.data.frame(score.ife.b[, c(2, 1)])
meta.ife.sb <- as.data.frame(score.ife.sb[, c(2, 1)])
meta.ife.bc <- as.data.frame(score.ife.bc[, c(2, 1)])
colnames(meta.ife.b) <- c("condition", "cellType")
colnames(meta.ife.sb) <- c("condition", "cellType")
colnames(meta.ife.bc) <- c("condition", "cellType")
row.names(meta.ife.b) <- paste0("C", 1:3)
row.names(meta.ife.sb) <- paste0("C", 1:3)
row.names(meta.ife.bc) <- paste0("C", 1:3)
score.ife.b <- t(apply(score.ife.b[, -c(1, 2)], 2, as.numeric))
score.ife.sb <- t(apply(score.ife.sb[, -c(1, 2)], 2, as.numeric))
score.ife.bc <- t(apply(score.ife.bc[, -c(1, 2)], 2, as.numeric))
colnames(score.ife.b) <- paste0("C", 1:3)
colnames(score.ife.sb) <- paste0("C", 1:3)
colnames(score.ife.bc) <- paste0("C", 1:3)

### Colors ----
annot_col2 <- list(
  condition = c("Control" = "#a6cde1",
                "OSKM_mCherry-neg" = "#008080",
                "OSKM_mCherry-pos" = "#ff7f50"),
  cellType = c("IFE.B" = "#f8766d", "IFE.SB" = "#619cff", "IFE.B.C" = "#00ba38")
)

### Plot ----
pheatmap(score.ife.b, scale = "row",
         cluster_rows = T, cluster_cols = F,
         annotation_col = meta.ife.b, labels_row = terms_,
         show_colnames = F, annotation_colors = annot_col2,
         color = colorRampPalette(c("royalblue", "white", "firebrick2"))(100))
pheatmap(score.ife.sb, scale = "row",
         cluster_rows = T, cluster_cols = F,
         annotation_col = meta.ife.sb, labels_row = terms_,
         show_colnames = F, annotation_colors = annot_col2,
         color = colorRampPalette(c("royalblue", "white", "firebrick2"))(100))
pheatmap(score.ife.bc, scale = "row",
         cluster_rows = T, cluster_cols = F,
         annotation_col = meta.ife.bc, labels_row = terms_,
         show_colnames = F, annotation_colors = annot_col2,
         color = colorRampPalette(c("royalblue", "white", "firebrick2"))(100))



# REACTOME Terms ----------------------------------------------------------
### Load data ----
skin2.ife <- readRDS("./rds/2024/0226/skin2_epi_ife.rds")
skin2.ife <- readRDS("./rds/2024/0415/skin2_ife.rds")
skin2.ife@meta.data %>% colnames()

### Get gene sets from msigdb ----
reac.cytokine <- fromJSON("./data/geneset/reactome/REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM.v2023.2.Mm.json")
reac.cytokine <- reac.cytokine$REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM$geneSymbols
reac.hypoxia <- fromJSON("./data/geneset/reactome/REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA.v2023.2.Mm.json")
reac.hypoxia <- reac.hypoxia$REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA$geneSymbols

### Calculate scores ----
skin2.ife <- skin2.ife %>% 
  AddModuleScore(features = list(reac.cytokine), name = "C") %>% 
  AddModuleScore(features = list(reac.hypoxia), name = "H")

### Plot ----
# Cytokine signaling in immune system
VlnPlot(skin2.ife, features = "C1", group.by = "condition", pt.size = 0) +
  geom_boxplot(width = 0.3) +
  theme(axis.title = element_blank()) + ggtitle(NULL)
# Cellular resonse to hypoxia
VlnPlot(skin2.ife, features = "H1", group.by = "condition", pt.size = 0) +
  geom_boxplot(width = 0.3) +
  theme(axis.title = element_blank()) + ggtitle(NULL)



#####
