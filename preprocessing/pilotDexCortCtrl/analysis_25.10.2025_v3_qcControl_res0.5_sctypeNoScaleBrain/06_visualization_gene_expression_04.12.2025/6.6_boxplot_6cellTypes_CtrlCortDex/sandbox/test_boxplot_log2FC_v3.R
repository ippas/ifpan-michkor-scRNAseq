library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# ============================================================================
# 1. LICZBA GENÓW > +1 i < −1  (Z NEST)
# ============================================================================

df_counts <- df_plot_fc %>%
  group_by(celltype, treatment) %>%
  nest() %>%
  mutate(
    n_up   = map_int(data, ~ sum(.x$log2FC >=  1, na.rm = TRUE)),
    n_down = map_int(data, ~ sum(.x$log2FC <= -1, na.rm = TRUE))
  ) %>%
  select(celltype, treatment, n_up, n_down) %>%
  ungroup()

df_labels_n <- bind_rows(
  df_counts %>%
    mutate(y =  2.5, label = paste0("n = ", n_up)),
  df_counts %>%
    mutate(y = -2.5, label = paste0("n = ", n_down))
)

# ============================================================================
# 2. STATYSTYKA BOXPLOTU: ŚREDNIA (MEAN)  (Z NEST)
# ============================================================================

df_stats <- df_plot_fc %>%
  group_by(celltype, treatment) %>%
  nest() %>%
  mutate(
    mean = map_dbl(data, ~ mean(.x$log2FC, na.rm = TRUE))
  ) %>%
  select(celltype, treatment, mean) %>%
  ungroup()

# ============================================================================
# 3. POZYCJA TEKSTU (STAŁA)
# ============================================================================

y_max <- max(df_plot_fc$log2FC, na.rm = TRUE)

df_labels_stats <- df_stats %>%
  mutate(
    y = 3,
    label = paste0("x̄ = ", round(mean, 2))
  )

# ============================================================================
# 4. FINALNY WYKRES
# ============================================================================

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
  
  ## LINIE REFERENCYJNE
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.6
  ) +
  
  ## n_up / n_down
  geom_text(
    data = df_labels_n,
    aes(
      x = celltype,
      y = y,
      label = label,
      group = interaction(celltype, treatment)
    ),
    color = "black",
    position = position_dodge(width = 0.75),
    size = 4,
    vjust = ifelse(df_labels_n$y > 0, 0, 1)
  ) +
  
  ## ŚREDNIA (JEDNA LINIA)
  geom_text(
    data = df_labels_stats,
    aes(
      x = celltype,
      y = y,
      label = label,
      group = interaction(celltype, treatment)
    ),
    color = "black",
    position = position_dodge(width = 0.75),
    size = 4,
    vjust = 0
  ) +
  
  ## KOLORY
  scale_color_manual(
    values = c(
      "CORT" = "firebrick2",
      "DEX"  = "steelblue3"
    )
  ) +
  
  ## ZAKRES OSI Y
  scale_y_continuous(
    limits = c(
      min(df_plot_fc$log2FC, na.rm = TRUE) - 0.3,
      y_max + 0.8
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
    legend.title = element_text(color = "black")
  ) +
  
  labs(
    x = "Cell type",
    y = "log2 fold change vs CTRL",
    color = "Treatment"
  )
