library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

df_plot_fc <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain[c(
  "Radial glial cells",
  "Astrocytes",
  "OPCs",
  "Oligodendrocytes",
  "Immature neurons",
  "Mature neurons"
)] %>%
  purrr::map(function(x) {
    x %>%
      filter(gene_symbol %in% res_heatmap$genes_per_cluster$cluster_1) %>%
      select(
        gene_symbol,
        log2FC_CORT_CPM,
        log2FC_DEX_CPM
      ) %>%
      pivot_longer(
        cols      = starts_with("log2FC_"),
        names_to  = "treatment",
        values_to = "log2FC"
      ) %>%
      mutate(
        treatment = case_when(
          treatment == "log2FC_CORT_CPM" ~ "CORT",
          treatment == "log2FC_DEX_CPM"  ~ "DEX"
        )
      )
  }) %>%
  bind_rows(.id = "celltype") %>%
  mutate(
    treatment = factor(treatment, levels = c("CORT", "DEX")),
    celltype  = factor(
      celltype,
      levels = c(
        "Radial glial cells",
        "Astrocytes",
        "OPCs",
        "Oligodendrocytes",
        "Immature neurons",
        "Mature neurons"
      )
    )
  )


ggplot(df_plot_fc) +
  geom_boxplot(
    aes(
      x = celltype,
      y = log2FC,
      color = treatment,
      group = interaction(celltype, treatment)
    ),
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    width = 0.6,
    linewidth = 0.6
  ) +
  geom_jitter(
    aes(
      x = celltype,
      y = log2FC,
      color = treatment
    ),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.75
    ),
    alpha = 0.4,
    size = 0.8
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  scale_color_manual(
    values = c(
      "CORT" = "firebrick2",
      "DEX"  = "steelblue3"
    )
  ) +
  theme_classic() +
  labs(
    x = "Cell type",
    y = "log2 fold change vs CTRL",
    color = "Treatment"
  )


# ##############################################################################

df_counts <- df_plot_fc %>%
  group_by(celltype, treatment) %>%
  nest() %>%
  mutate(
    n_up = map_int(data, ~ sum(.x$log2FC >=  1, na.rm = TRUE)),
    n_down = map_int(data, ~ sum(.x$log2FC <= -1, na.rm = TRUE))
  ) %>%
  select(celltype, treatment, n_up, n_down) %>%
  ungroup()

df_labels <- bind_rows(
  df_counts %>%
    mutate(
      y =  2.5,
      label = paste0("n = ", n_up)
    ),
  df_counts %>%
    mutate(
      y = -2.5,
      label = paste0("n = ", n_down)
    )
)


ggplot(df_plot_fc) +
  ## BOXPLOT
  geom_boxplot(
    aes(
      x = celltype,
      y = log2FC,
      color = treatment,
      group = interaction(celltype, treatment)
    ),
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    width = 0.6,
    linewidth = 0.6
  ) +
  
  ## PUNKTY
  geom_jitter(
    aes(
      x = celltype,
      y = log2FC,
      color = treatment
    ),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.75
    ),
    alpha = 0.35,
    size = 0.8
  ) +
  
  ## PROGI ±1
  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.6
  ) +
  
  ## LINIA BAZOWA
  geom_hline(
    yintercept = 0,
    color = "grey60",
    linewidth = 0.5
  ) +
  
  ## 🔢 LICZBY GENÓW NAD / POD PROGIEM
  geom_text(
    data = df_labels,
    aes(
      x = celltype,
      y = y,
      label = label,
      group = interaction(celltype, treatment)
    ),
    color = "black",
    position = position_dodge(width = 0.75),
    size = 4,
    vjust = ifelse(df_labels$y > 0, 0, 1)
  ) +
  
  scale_color_manual(
    values = c(
      "CORT" = "firebrick2",
      "DEX"  = "steelblue3"
    )
  ) +
  
  theme_classic() +
  
  ## 🔒 WYMUSZENIE CZARNEGO TEKSTU
  theme(
    text              = element_text(color = "black"),
    axis.text.x       = element_text(color = "black"),
    axis.text.y       = element_text(color = "black"),
    axis.title.x      = element_text(color = "black"),
    axis.title.y      = element_text(color = "black"),
    legend.text       = element_text(color = "black"),
    legend.title      = element_text(color = "black")
  ) +
  
  labs(
    x = "Cell type",
    y = "log2 fold change vs CTRL",
    color = "Treatment"
  )

