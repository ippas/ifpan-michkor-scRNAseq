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
      filter(gene_symbol %in% res_heatmap$genes_per_cluster$cluster_4) %>%
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


# ============================================================================
# 0. PARAMETRY GLOBALNE
# ============================================================================

y_max <- max(df_plot_fc$log2FC, na.rm = TRUE)
y_min <- min(df_plot_fc$log2FC, na.rm = TRUE)

n_genes_total <- df_plot_fc %>%
  distinct(gene_symbol) %>%
  nrow()

# ============================================================================
# 1. n_up / n_down (|log2FC| >= 1) ---- NEST
#     -> DALEJ od linii
# ============================================================================

df_counts_thresh <- df_plot_fc %>%
  group_by(celltype, treatment) %>%
  nest() %>%
  mutate(
    n_up   = map_int(data, ~ sum(.x$log2FC >=  1, na.rm = TRUE)),
    n_down = map_int(data, ~ sum(.x$log2FC <= -1, na.rm = TRUE))
  ) %>%
  select(celltype, treatment, n_up, n_down) %>%
  ungroup()

df_labels_thresh <- bind_rows(
  df_counts_thresh %>%
    mutate(
      y =  1.85,                 # dalej nad +1
      label = paste0("n = ", n_up)
    ),
  df_counts_thresh %>%
    mutate(
      y = -1.85,                 # dalej pod −1
      label = paste0("n = ", n_down)
    )
)

# ============================================================================
# 2. ALL UP / ALL DOWN (log2FC > 0 / < 0) ---- NEST
#     -> BLIŻEJ linii przerywanej
# ============================================================================

df_counts_all <- df_plot_fc %>%
  group_by(celltype, treatment) %>%
  nest() %>%
  mutate(
    n_all_up   = map_int(data, ~ sum(.x$log2FC >  0, na.rm = TRUE)),
    n_all_down = map_int(data, ~ sum(.x$log2FC <  0, na.rm = TRUE))
  ) %>%
  select(celltype, treatment, n_all_up, n_all_down) %>%
  ungroup()

df_labels_all <- bind_rows(
  df_counts_all %>%
    mutate(
      y =  1.45,                 # bliżej +1
      label = paste0("all↑ = ", n_all_up)
    ),
  df_counts_all %>%
    mutate(
      y = -1.45,                 # bliżej −1
      label = paste0("all↓ = ", n_all_down)
    )
)

# ============================================================================
# 3. ŚREDNIA (x̄) ---- NEST
# ============================================================================

df_stats_mean <- df_plot_fc %>%
  group_by(celltype, treatment) %>%
  nest() %>%
  mutate(
    mean_log2FC = map_dbl(data, ~ mean(.x$log2FC, na.rm = TRUE))
  ) %>%
  select(celltype, treatment, mean_log2FC) %>%
  ungroup() %>%
  mutate(
    y = 3,
    label = paste0("x̄ = ", round(mean_log2FC, 2))
  )

# ============================================================================
# 4. WYKRES
# ============================================================================

ggplot(df_plot_fc, aes(x = celltype, y = log2FC, color = treatment)) +
  
  ## BOXPLOT
  geom_boxplot(
    aes(group = interaction(celltype, treatment)),
    position = position_dodge(width = 0.75),
    width = 0.6,
    linewidth = 0.6,
    outlier.shape = NA
  ) +
  
  ## JITTER
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.75
    ),
    alpha = 0.35,
    size = 0.8
  ) +
  
  ## LINIE REFERENCYJNE
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.6
  ) +
  
  ## ALL UP / ALL DOWN (BLIŻEJ LINII)
  geom_text(
    data = df_labels_all,
    aes(
      x = celltype,
      y = y,
      label = label,
      group = interaction(celltype, treatment)
    ),
    position = position_dodge(width = 0.75),
    size = 4,
    color = "black",
    vjust = ifelse(df_labels_all$y > 0, 0, 1)
  ) +
  
  ## n_up / n_down (DALEJ)
  geom_text(
    data = df_labels_thresh,
    aes(
      x = celltype,
      y = y,
      label = label,
      group = interaction(celltype, treatment)
    ),
    position = position_dodge(width = 0.75),
    size = 4,
    color = "black",
    vjust = ifelse(df_labels_thresh$y > 0, 0, 1)
  ) +
  
  ## ŚREDNIA
  geom_text(
    data = df_stats_mean,
    aes(
      x = celltype,
      y = y,
      label = label,
      group = interaction(celltype, treatment)
    ),
    position = position_dodge(width = 0.75),
    size = 4,
    color = "black",
    vjust = 0
  ) +
  
  ## KOLORY
  scale_color_manual(
    values = c(
      "CORT" = "firebrick2",
      "DEX"  = "steelblue3"
    )
  ) +
  
  ## ZAKRES Y
  scale_y_continuous(
    limits = c(
      y_min - 0.4,
      y_max + 1.0
    )
  ) +
  
  ## THEME
  theme_classic() +
  theme(
    text         = element_text(color = "black"),
    axis.text.x  = element_text(color = "black"),
    axis.text.y  = element_text(color = "black"),
    axis.title   = element_text(color = "black"),
    legend.text  = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  
  ## OPISY
  labs(
    title = paste0(
      "Distribution of log2FC values (n = ",
      n_genes_total,
      " genes)"
    ),
    x = "Cell type",
    y = "log2 fold change vs CTRL",
    color = "Treatment"
  )
