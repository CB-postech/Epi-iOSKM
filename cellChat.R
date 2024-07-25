#######
# Script used to perform cell-cell interaction comparison 
# across conditions using the first dataset
#######


set.seed(1234)
# Libraries ---------------------------------------------------------------
library("CellChat")



# Read data ---------------------------------------------------------------
skin1.ife <- readRDS("./rds/2024/0415/skin1_ife.rds")



# Split Data by Condition -------------------------------------------------
skin1.ife.ctrl <- skin1.ife[, skin1.ife$oskm == "control"] %>% 
  NormalizeData()
skin1.ife.oskm <- skin1.ife[, skin1.ife$oskm == "oskm"] %>% 
  NormalizeData()



# Configure CellChat DB to use --------------------------------------------
ccDB <- CellChatDB.mouse
ccDB.use <- subsetDB(ccDB)



# Process CCC for Each Data -----------------------------------------------
### Configure CellChat DB to use ----
### function ----
process_ccc <- function(seurat, assay = "RNA", group.by = NULL) {
  if (is.null(group.by)) {
    seurat$cluster <- paste0("g", Idents(seurat))
    group.by <- "cluster"
  }
  # create cc object
  cc <- createCellChat(object = seurat, group.by = group.by, assay = assay)
  cc@DB <- ccDB.use
  # pre-processing
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, group.by = group.by)
  cc <- identifyOverExpressedInteractions(cc)
  # inference
  cc <- computeCommunProb(cc, type = "triMean")
  # infer ccc at a signaling pathway level
  cc <- computeCommunProbPathway(cc)
  # return
  return(cc)
}
cc.ctrl <- process_ccc(skin1.ife.ctrl, assay = "originalexp")
cc.oskm <- process_ccc(skin1.ife.oskm, assay = "originalexp")



# Merge CC Objects --------------------------------------------------------
cc <- mergeCellChat(object = list(CTRL = cc.ctrl, OSKM = cc.oskm),
                    add.names = c("CTRL", "OSKM"))



# Identify altered signaling with distinct interaction strength -----------
### Compare overall information flow of each sig/LP pair ----
rankNet(cc, mode = "comparison", measure = "weight", stacked = T, do.stat = T,
        color.use = c("cornflowerblue", "indianred2"))



#####