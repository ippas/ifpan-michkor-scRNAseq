library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(forcats)
library(ggtext)

# ===============================
# 1. Lista neuronalnych terminów
# ===============================

neural_terms <- c(
  "Neuron Differentiation (GO:0030182)",
  "Positive Regulation of Neuron Projection Development (GO:0010976)",
  "Axon Guidance (GO:0007411)",
  "Nervous System Development (GO:0007399)",
  "Neuron Projection Guidance (GO:0097485)",
  "Glutamatergic Neuron Differentiation (GO:1905962)",
  "Cerebellar Granule Cell Differentiation (GO:0021707)",
  "Regulation of Synapse Assembly (GO:0051963)",
  "Positive Regulation of Synapse Assembly (GO:0051965)"
)

# ===============================
# 2. Przetwarzanie GO
# ===============================

df_go <- enrichr_results$GO_Biological_Process_2025 %>%
  bind_rows(.id = "cluster") %>%
  filter(n_genes > 0) %>%
  group_by(cluster) %>%
  arrange(pvalue) %>%
  slice_head(n = 5) %>%
  mutate(
    neglog10_p = -log10(pvalue),
    cluster = factor(cluster, levels = paste0("Cluster_", 1:5)),
    term = as.character(term),
    term_type = ifelse(term %in% neural_terms, "neural", "other"),
    term_label_raw = case_when(
      term_type == "neural" ~ paste0("<b><span style='color:black;'>", term, "</span></b>"),
      TRUE                  ~ paste0("<span style='color:#555555;'>", term, "</span>")
    ),
    term_label = factor(term_label_raw, levels = rev(term_label_raw))
  ) %>%
  ungroup()

# ===============================
# 3. Kolory
# ===============================

cluster_colors <- c(
  "Cluster_1" = "#6BAED6",
  "Cluster_2" = "#3182BD",
  "Cluster_3" = "#08519C",
  "Cluster_4" = "#CB181D",
  "Cluster_5" = "#FB6A4A"
)

# ===============================
# 4. WYKRES – ręczna oś X
# ===============================

# maksymalna wartość p-value na potrzeby osi X
x_breaks <- c(0, 2, 4, 6)
x_labels <- x_breaks

p <- ggplot(df_go, aes(x = term_label, y = neglog10_p, fill = cluster)) +
  geom_col(color = "black", width = 0.75) +
  geom_text(aes(label = n_genes), hjust = -0.2, size = 3.5) +
  
  facet_grid(
    rows = vars(cluster),
    switch = "y",
    scales = "free_y",
    labeller = labeller(cluster = function(x) x)
  ) +
  
  scale_fill_manual(values = cluster_colors) +
  
  labs(
    x = NULL,
    y = "GO Biological Process terms",
    title = "Top 5 enriched GO Biological Processes per GR-dependent gene cluster"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(size = 12, color = "black", angle = 0),
    axis.text.y = element_markdown(size = 12),
    axis.text.x = element_text(size = 11),
    legend.position = "none"
  ) +
  
  coord_flip() +
  
  # <<< RĘCZNE DODANIE OSI X DO KAŻDEGO PANELU >>>
  annotate("segment",
           x = Inf, xend = Inf,
           y = x_breaks, yend = x_breaks,
           colour = "black") +
  annotate("text",
           x = Inf,
           y = x_breaks,
           label = x_labels,
           hjust = -0.3,
           size = 3.2)

p

