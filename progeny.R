#######
# Script used for pathway activation score inference using PROGENy
#######


set.seed(1234)
# Libraries ---------------------------------------------------------------
library("progeny")



# Load data ----------------------------------------------------------------
skin2.ife <- readRDS("./rds/skin2_ife.rds")



# Progeny -----------------------------------------------------------------
skin2.ife <- progeny(skin2.ife, scale = F, organism = "Mouse",
                     top = 500, perm = 1, return_assay = T,
                     assay_name = "RNA")
skin2.ife <- ScaleData(skin2.ife, assay = "progeny")



# Summarize progeny scores ------------------------------------------------
progeny.df <- skin2.ife %>% 
  GetAssayData(assay = "progeny", layer = "scale.data") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)
progeny.df <- progeny.df %>% 
  inner_join(., 
             data.frame(Condition = skin2.ife$condition) %>% 
               rownames_to_column("Cell"),
             by = "Cell")
progeny.df.sum <- progeny.df %>% 
  group_by(Pathway, Condition) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
progeny.df.sum <- progeny.df.sum %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)



# Plot --------------------------------------------------------------------
# color palette
paletteLength <- 100
myColor = colorRampPalette(c("Darkblue", "white", "red"))(paletteLength)
progenyBreaks <- c(seq(min(progeny.df.sum), 0, 
                      length.out = ceiling(paletteLength / 2) + 1),
                  seq(max(progeny.df.sum) / paletteLength, 
                      max(progeny.df.sum), 
                      length.out = floor(paletteLength / 2)))
progeny.hmap <- pheatmap(t(progeny.df.sum), fontsize=14, 
                        fontsize_row = 10, scale = "none",
                        color = myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0, treeheight_row = 75,
                        border_color = NA,
                        cluster_cols = F, cluster_rows = T)



#####
