library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(forcats)

cutoff <- -log10(0.05)

# 1) WYCIĄGNIĘCIE DANYCH DO JEDNEJ TABELI
df_plot <- enrichr_results$DisGeNET %>%
  bind_rows(.id = "cluster") %>%
  filter(term %in% c("Bipolar Disorder", "Unipolar Depression", "Schizophrenia"),
         n_genes > 0) %>%
  mutate(
    neglog10_p = -log10(pvalue),
    term = factor(term, levels = c("Schizophrenia", "Bipolar Disorder", "Unipolar Depression")),
    cluster = factor(cluster, levels = paste0("Cluster_", 1:5))
  )

ggplot(df_plot, aes(x = term, y = neglog10_p, fill = term)) +
  geom_col(color = "black", width = 0.75) +
  
  # ---- N GENÓW NAD SŁUPKAMI ----
geom_text(aes(label = n_genes),
          vjust = -0.3,
          size = 4,
          color = "black") +  # <<<<<< CZARNY TEKST
  
  # ---- CZERWONA LINIA –log10(0.05) ----
geom_hline(yintercept = cutoff,
           linetype = "dashed",
           color = "red",
           linewidth = 0.7) +
  
  facet_wrap(vars(cluster), nrow = 1) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "",
    y = "-log10(p-value)",
    title = "DisGeNET enrichments for GR-dependent clusters"
  ) +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(color = "black"),           # <<<<<< CZARNY TEKST GLOBALNIE
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    plot.title  = element_text(color = "black"),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 12, color = "black")  # <<<<<< CZARNY
  ) +
  expand_limits(y = max(df_plot$neglog10_p, cutoff) + 0.5)