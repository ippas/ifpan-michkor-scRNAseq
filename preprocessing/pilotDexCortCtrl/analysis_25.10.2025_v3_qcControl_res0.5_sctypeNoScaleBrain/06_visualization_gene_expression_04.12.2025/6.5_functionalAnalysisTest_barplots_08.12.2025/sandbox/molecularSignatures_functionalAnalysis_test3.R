library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(forcats)
library(ggtext)   # <<< potrzebne do HTML w etykietach osi

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
  "Positive Regulation of Synapse Assembly (GO:0051965)"
)

# ============================================
# 2) Wyciągamy GO BP i wybieramy TOP 5 per cluster
# ============================================

df_go <- enrichr_results$GO_Biological_Process_2025 %>%
  bind_rows(.id = "cluster") %>%
  filter(n_genes > 0) %>%
  group_by(cluster) %>%
  arrange(pvalue) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    neglog10_p = -log10(pvalue),
    cluster = factor(cluster, levels = paste0("Cluster_", 1:5)),
    term = as.character(term),
    
    # neuronalne vs pozostałe
    term_type = ifelse(term %in% neural_terms, "neural", "other"),
    
    # HTML dla ggtext: neuronalne czarne + bold, pozostałe ciemny szary (#555555)
    term_label = case_when(
      term_type == "neural" ~ paste0("<b><span style='color:black;'>", term, "</span></b>"),
      TRUE ~ paste0("<span style='color:#555555;'>", term, "</span>")
    ),
    
    term_label = factor(term_label, levels = unique(term_label))
  )

# ============================================
# 3) Paleta kolorów per cluster
# ============================================

cluster_colors <- c(
  "Cluster_1" = "#6BAED6",
  "Cluster_2" = "#3182BD",
  "Cluster_3" = "#08519C",
  "Cluster_4" = "#CB181D",
  "Cluster_5" = "#FB6A4A"
)

# ============================================
# 4) Wykres z HTML w etykietach osi X (neuronalne czarne, inne ciemne szare)
# ============================================

ggplot(df_go, aes(x = term_label, y = neglog10_p, fill = cluster)) +
  geom_col(color = "black", width = 0.75) +
  
  geom_text(aes(label = n_genes),
            vjust = -0.3,
            size = 3.5,
            color = "black") +
  
  
  facet_wrap(vars(cluster), nrow = 1, scales = "free_x") +
  scale_fill_manual(values = cluster_colors) +
  
  labs(
    x = "",
    y = "-log10(p-value)",
    title = "Top 5 enriched GO Biological Processes per GR-dependent gene cluster"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    text = element_text(color = "black"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 12, color = "black"),
    
    # obsługa HTML w etykietach osi X z kolorami i bold
    axis.text.x = element_markdown(angle = 45, hjust = 1, size = 9),
    
    legend.position = "none"
  ) +
  
  expand_limits(y = max(df_go$neglog10_p, cutoff) + 0.3)

