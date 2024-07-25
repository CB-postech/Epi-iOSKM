#######
# Scripts used for the enrichment analysis using either REACTOME or SMART DB
#######


set.seed(1234)
# Libraries ---------------------------------------------------------------
library("gprofiler2")
library("STRINGdb")



# Load data ---------------------------------------------------------------
skin2.ife <- readRDS("./rds/skin2_ife.rds")
gocc.exs <- fromJSON("./data/geneSet/GOCC_EXTRACELLULAR_SPACE.v2023.2.Mm.json")
gocc.exs <- gocc.exs$GOCC_EXTRACELLULAR_SPACE$geneSymbols



# Get extracellular genes from DEGs ---------------------------------------
# get differential expression
f <- FindMarkers(skin2.ife, group.by = "condition", 
                 ident.1 = "OSKM_mCherry-pos", ident.2 = "Control",
                 min.pct = 0, logfc.threshold = 0)
# subset extracellular genes from DEGs
f.up <- f %>% 
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 1.2,
    row.names(.) %in% gocc.exs
  ) %>% 
  arrange(p_val_adj, desc(avg_log2FC))



# GO enrichment -----------------------------------------------------------
# enrichment
go <- gprofiler2::gost(
  query = row.names(f.up), organism = "mmusculus", 
  correction_method = "bonferroni", sources = "REAC")
go.res <- go$result; go.meta <- go$meta
# plot
go.res %>%
  arrange(p_value) %>% 
  select(term_name, p_value) %>%
  mutate(term_name = factor(term_name, levels = rev(term_name))) %>% 
  ggplot() +
  geom_col(
    aes(x = term_name, y = -log10(p_value)),
    fill = "firebrick3"
  ) +
  theme_classic() +
  labs(x = "", y = "-log10(p.value)") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  theme(
    axis.title = element_text(size = 17),
    axis.text.x = element_text(),
    axis.text = element_text(size = 14),
    panel.background = element_blank(),
    plot.margin = margin(.5, .5, .5, .5, "cm")
  )



# SMART enrichment --------------------------------------------------------
# initialize STRINGdb object
string.db <- STRINGdb::STRINGdb$new(
  version = "11.5", species = 10090, score_threshold = 400, network_type = "full",
  input_directory = "", protocol = "http")
# retrieve domain annotations
domain.info <- string.db$get_enrichment(row.names(f.up))
# plot
domain.info %>%
  filter(category == "SMART") %>% 
  mutate(value = -log10(fdr)) %>%
  arrange(desc(value)) %>%
  mutate(description = gsub("\\.", "", description)) %>%
  mutate(description = factor(description, levels = rev(description))) %>%
  ggplot() +
  geom_col(
    aes(x = description, y = value),
    fill = "#cdb2cf", color = "black", lwd = 0.5
  ) +
  theme_classic() +
  labs(x = "", y = "-log10(adj.p)",
       title = "Protein domains (SMART)") +
  scale_y_continuous(limits = c(0, 7.2), expand = c(0, 0)) +
  coord_flip() +
  theme(
    plot.title = element_text(size = 17),
    axis.title = element_text(size = 17),
    axis.text.x = element_text(),
    axis.text = element_text(size = 14),
    panel.background = element_blank(),
    plot.margin = margin(.5, .5, .5, .5, "cm")
  )



#####
