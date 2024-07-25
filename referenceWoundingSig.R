#######
# Scripts used for identifying DEGs that represent wounding signature
# using the reference data
#######


set.seed(1234)
# Load data ---------------------------------------------------------------
daniel.ife <- readRDS("./rds/2024/0610/daniel_ife.rds")



# DEGs --------------------------------------------------------------------
df.plt <- FindMarkers(
  daniel.ife, group.by = "condition", 
  ident.1 = "wounded", ident.2 = "un-wounded",
  logfc.threshold = 0, min.pct = 0)



# Volcano plot ------------------------------------------------------------
### Prepare data.frame for plotting ----
df.plt$gene <- row.names(df.plt)
min(df.plt$p_val_adj[df.plt$p_val_adj != 0])  # 3.08e-282
df.plt$p_val_adj[df.plt$p_val_adj == 0] <- min(df.plt$p_val_adj[df.plt$p_val_adj != 0])

### Select genes of interest ----
c.ccl <- c(
  "Ccl3", "Ccl4", "Ccl5", "Ccl6", "Ccl7", "Ccl8", "Ccl9", "Ccl11",
  "Ccl20", "Ccl22", "Cxcl3", "Cxcl5", "Cxcl10", "Cxcl16", 
  "Il1a", "Il1b", "Il6", "Tnf"
)


### Set threshold for significance ----
thr.p <- 0.05
thr.lfc <- 1.2

### Plot ----
ggplot() +
  # plot background genes
  geom_point(
    data = df.plt %>% filter(p_val_adj > thr.p | abs(avg_log2FC) < thr.lfc),
    aes(
      x = avg_log2FC, y = -log10(p_val_adj)
    ),
    color = "grey"
  ) +
  # plot Up-regulated genes
  geom_point(
    data = df.plt[df.plt$p_val_adj < thr.p & df.plt$avg_log2FC > thr.lfc, ],
    aes(
      x = avg_log2FC, y = -log10(p_val_adj)
    ),
    color = "firebrick2", alpha = 0.5
  ) +
  # plot Down-regulated genes
  geom_point(
    data = df.plt[df.plt$p_val_adj < thr.p & df.plt$avg_log2FC < -thr.lfc, ],
    aes(
      x = avg_log2FC, y = -log10(p_val_adj)
    ),
    color = "royalblue", alpha = 0.5
  ) +
  # annotate Ccl genes
  geom_text_repel(
    data = df.plt[df.plt$avg_log2FC > thr.lfc & df.plt$gene %in% c.ccl, ],
    aes(
      x = avg_log2FC, y = -log10(p_val_adj), label = gene
    ),
    segment.colour = "black", min.segment.length = 0,
    color = "black", size = 7, force = 30,
    max.overlaps = Inf
  ) +
  theme_classic() + labs(
    x = "log2FC (Wounded vs Unwounded)", y = "-log10(adj.p)"
  ) +
  theme(
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 14),
    panel.background = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  scale_x_continuous(limits = c(-10, 20), expand = c(0, 0))



#####