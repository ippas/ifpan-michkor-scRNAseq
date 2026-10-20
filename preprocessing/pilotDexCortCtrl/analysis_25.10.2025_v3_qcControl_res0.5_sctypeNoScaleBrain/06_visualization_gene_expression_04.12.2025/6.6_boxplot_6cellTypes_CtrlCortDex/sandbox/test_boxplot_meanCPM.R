# ##############################################################################
# ---- prepare data ----
# ##############################################################################
combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$Astrocytes %>% 
  filter(gene_symbol %in% res_heatmap$genes_per_cluster$cluster_1) %>% head
  select(gene_symbol, meanCPM_CTRL, meanCPM_CORT, meanCPM_DEX)

res_heatmap$genes_per_cluster$cluster_1

 
library(ggplot2)
library(dplyr)

df_plot <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain[c(
  "Astrocytes",
  "Radial glial cells",
  "OPCs",
  "Oligodendrocytes",
  "Immature neurons",               
  "Mature neurons" 
)] %>%
  purrr::map(function(x) {
    x %>% 
      filter(gene_symbol %in% res_heatmap$genes_per_cluster$cluster_1) %>% 
      select(gene_symbol, meanCPM_CTRL, meanCPM_CORT, meanCPM_DEX) %>% 
      tidyr::pivot_longer(
        cols      = starts_with("meanCPM_"),
        names_to  = "treatment",
        values_to = "meanCPM"
      ) %>% 
      mutate(
        treatment = sub("^meanCPM_", "", treatment)
      )
  }) %>%
  bind_rows(.id = "celltype") %>%
  mutate(
    treatment = factor(treatment, levels = c("CTRL", "CORT", "DEX")),
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

# pomocnicza zmienna numeryczna tylko do linii/smooth
df_plot <- df_plot %>%
  mutate(
    treatment_num = as.numeric(factor(
      treatment,
      levels = c("CTRL", "CORT", "DEX")
    ))
  )

ggplot(df_plot) +
  ## boxplot (kategoryczna oś X)
  geom_boxplot(
    aes(
      x = treatment,
      y = log2(meanCPM),
      color = treatment
    ),
    outlier.shape = NA,
    width = 0.6,
    fill = NA,
    linewidth = 0.6
  ) +
  ## punkty
  geom_jitter(
    aes(
      x = treatment,
      y = log2(meanCPM),
      color = treatment
    ),
    width = 0.15,
    alpha = 0.5,
    size = 1
  ) +
  ## linie łączące ten sam gen
  geom_line(
    aes(
      x = treatment_num,
      y = log2(meanCPM),
      group = gene_symbol
    ),
    color = "grey75",
    linewidth = 0.5,
    alpha = 0.8
  ) +
  ## smooth (trend)
  geom_smooth(
    aes(
      x = treatment_num,
      y = log2(meanCPM),
      group = 1
    ),
    method = "loess",
    se = FALSE,
    color = "black",
    linewidth = 0.8
  ) +
  facet_wrap(~ celltype, scales = "free_y") +
  scale_x_discrete(limits = c("CTRL", "CORT", "DEX")) +
  scale_color_manual(
    values = c(
      "CTRL" = "grey40",
      "CORT" = "firebrick2",
      "DEX"  = "steelblue3"
    )
  ) +
  theme_classic() +
  labs(
    x = "Treatment",
    y = "log2(Mean CPM)"
  )





# ##############################################################################

df_plot <- df_plot %>%
  mutate(
    log2meanCPM = log2(meanCPM),
    treatment   = factor(treatment, levels = c("CTRL", "CORT", "DEX")),
    celltype    = factor(celltype)
  )

ggplot(df_plot) +
  
  ## BOXPLOTY – celltype obok siebie w ramach treatment
  geom_boxplot(
    aes(
      x = treatment,
      y = log2meanCPM,
      color = celltype,
      group = interaction(treatment, celltype)
    ),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.65,
    fill = NA,
    linewidth = 0.6
  ) +
  
  ## PUNKTY
  geom_jitter(
    aes(
      x = treatment,
      y = log2meanCPM,
      color = celltype
    ),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.8
    ),
    alpha = 0.4,
    size = 0.8
  ) +
  
  ## LINIE – łączą ten sam gen w obrębie celltype
  geom_line(
    aes(
      x = as.numeric(treatment),
      y = log2meanCPM,
      group = interaction(gene_symbol, celltype)
    ),
    color = "grey80",
    linewidth = 0.4,
    alpha = 0.6
  ) +
  
  ## TREND globalny (opcjonalnie – raczej delikatny)
  geom_smooth(
    aes(
      x = as.numeric(treatment),
      y = log2meanCPM,
      group = 1
    ),
    method = "loess",
    se = FALSE,
    color = "black",
    linewidth = 0.8
  ) +
  
  scale_x_discrete(limits = c("CTRL", "CORT", "DEX")) +
  
  theme_classic() +
  labs(
    x = "Treatment",
    y = "log2(Mean CPM)",
    color = "Cell type"
  )
