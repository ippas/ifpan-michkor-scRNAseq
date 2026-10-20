library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(forcats)
library(ggtext)

cutoff <- -log10(0.05)

# ============================================
# 1) Lista neuronalnych terminów
# ============================================

neural_terms <- c(
  "Neuron Differentiation (GO:0030182)",
  "Positive Regulation of Neuron Projection Development (GO:0010976)",
  "Axon Guidance (GO:0007411)",
  "Nervous System Development (GO:0007399)",
  "Neuron Projection Guidance (GO:0097485)",
  "Glutamatergic Neuron Differentiation (GO:1905962)",
  "Cerebellar Granule Cell Differentiation (GO:0021707)",
  "Regulation of Synapse Assembly (GO:0051963)",
  "Positive Regulation of Synapse Assembly (GO:0051965)",
  "Axonogenesis (GO:0007409)"
)

# ============================================
# 2) Wyciągamy GO BP i wybieramy TOP 5 per cluster
#    + ODWRÓCONE SORTOWANIE (najmniej istotne u góry)
# ============================================

df_go <- enrichr_results$GO_Biological_Process_2025 %>%
  bind_rows(.id = "cluster") %>%
  filter(n_genes > 0) %>%
  group_by(cluster) %>%
  arrange(pvalue) %>%                 # sortujemy od najlepszych
  slice_head(n = 5) %>%
  mutate(
    neglog10_p = -log10(pvalue),
    cluster = factor(cluster, levels = paste0("cluster_", 1:5)),
    term = as.character(term),
    
    term_type = ifelse(term %in% neural_terms, "neural", "other"),
    
    term_label_raw = case_when(
      term_type == "neural" ~ paste0("<b><span style='color:black;'>", term, "</span></b>"),
      TRUE                  ~ paste0("<span style='color:#555555;'>", term, "</span>")
    )
  ) %>%
  # <<< KLUCZOWA ZMIANA → ODWRÓCENIE kolejności LEVELS >>>
  mutate(term_label = factor(term_label_raw, levels = rev(term_label_raw))) %>%
  ungroup()

# ============================================
# 3) Paleta kolorów per cluster
# ============================================

cluster_colors <- c(
  "cluster_1" = "#6BAED6",
  "cluster_2" = "#3182BD",
  "cluster_3" = "#08519C",
  "cluster_4" = "#CB181D",
  "cluster_5" = "#FB6A4A"
)


cluster_colors <- c(
  "cluster_1" = "#FB6A4A",
  "cluster_2" = "#FB6A4A",
  "cluster_3" = "#FB6A4A",
  "cluster_4" = "#FB6A4A",
  "cluster_5" = "#FB6A4A"
)

cluster_colors <- c(
  "cluster_1" = "#6BAED6",
  "cluster_2" = "#6BAED6",
  "cluster_3" = "#6BAED6",
  "cluster_4" = "#6BAED6",
  "cluster_5" = "#6BAED6"
)

# ============================================
# 4) Wykres — termy posortowane od NAJMNIEJ do NAJBARDZIEJ istotnych
# ============================================

ggplot(df_go, aes(x = term_label, y = neglog10_p, fill = cluster)) +
  geom_col(color = "black", width = 0.75) +
  
  geom_text(
    aes(label = n_genes),
    hjust = -0.2,
    size = 3.5,
    color = "black"
  ) +
  
  facet_grid(
    rows = vars(cluster),
    switch = "y",
    scales = "free_y"
  ) +
  
  scale_fill_manual(values = cluster_colors) +
  
  labs(
    x = "",
    y = "-log10(p-value)",
    title = "Top 5 enriched GO Biological Processes per GR-dependent gene cluster"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    text = element_text(color = "black"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(size = 12, color = "black", angle = 90),
    axis.text.y = element_markdown(size = 12),
    legend.position = "none"
  ) +
  
  expand_limits(y = max(df_go$neglog10_p) + 0.3) +
  coord_flip() -> p



svg("results/figures/barplot_functionalAnalysisTest_08.12.2025/top5_GObiologicalProces_clusters_v3_09.12.2025.svg", width = 10, height = 10)
print(p)       # albo print(last_plot())
dev.off()



library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(forcats)

cutoff <- -log10(0.05)

# ============================================
# GO BP — TOP 5 per GROUP (cluster_2_3, cluster_4_5)
# ============================================

df_go <- enrichr_results$GO_Biological_Process_2025 %>%
  bind_rows(.id = "cluster_group") %>%
  filter(n_genes > 0) %>%
  group_by(cluster_group) %>%
  arrange(pvalue) %>%
  slice_head(n = 5) %>%
  mutate(
    neglog10_p = -log10(pvalue),
    cluster_group = factor(
      cluster_group,
      levels = c("cluster_2_3", "cluster_4_5")
    ),
    term_label = factor(term, levels = rev(term))
  ) %>%
  ungroup()

# ============================================
# Kolory per GROUP
# ============================================

cluster_colors <- c(
  "cluster_2_3" = "#6BAED6",
  "cluster_4_5" = "#6BAED6"
)

# ============================================
# Wykres
# ============================================

ggplot(df_go, aes(x = term_label, y = neglog10_p, fill = cluster_group)) +
  geom_col(color = "black", width = 0.75) +
  
  geom_text(
    aes(label = n_genes),
    hjust = -0.2,
    size = 3.5,
    color = "black"
  ) +
  
  facet_grid(
    rows = vars(cluster_group),
    switch = "y",
    scales = "free_y"
  ) +
  
  scale_fill_manual(values = cluster_colors) +
  
  labs(
    x = "",
    y = "-log10(p-value)",
    title = "Top 5 enriched GO Biological Processes per GR-dependent cluster group"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(size = 12, color = "black", angle = 90),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  ) +
  
  expand_limits(y = max(df_go$neglog10_p) + 0.3) +
  coord_flip() -> p

svg(
  "results/figures/barplot_functionalAnalysisTest_08.12.2025/top5_GObiologicalProcess_clusterGroups_2-3_vs_4-5_simple.svg",
  width = 10,
  height = 8
)
print(p)
dev.off()
