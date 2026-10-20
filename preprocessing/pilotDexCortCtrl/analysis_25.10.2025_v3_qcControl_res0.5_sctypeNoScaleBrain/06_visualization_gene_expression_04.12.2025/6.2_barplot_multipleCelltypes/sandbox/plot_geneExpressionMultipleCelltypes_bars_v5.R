###############################################################
### A) PACKAGES
###############################################################
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)
library(tibble)

###############################################################
### B) PARAMETR — ZMIENIASZ TYLKO TO
###############################################################
gene_symbol <- "ZBTB16"   # można zmienić na dowolny gen

###############################################################
### C) WYCIĄGNIĘCIE GENU DLA ASTROCYTÓW
###############################################################
df_gene <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$Astrocytes %>%
  filter(gene_symbol == !!gene_symbol)

if(nrow(df_gene) == 0){
  stop(paste("Brak genu:", gene_symbol))
}

###############################################################
### D) SAMPLE MAP
###############################################################
sample_map <- integrated_pilotDexCortCtrl_filtered@meta.data %>%
  select(bc1_well, treatment) %>%
  distinct()

###############################################################
### E) LONG FORMAT (WSZYSTKIE A1–A12 → CPM)
###############################################################
df_long_ast <- df_gene %>%
  select(gene_symbol, starts_with("A")) %>%
  pivot_longer(
    cols = starts_with("A"),
    names_to = "bc1_well",
    values_to = "CPM"
  ) %>%
  left_join(sample_map, by = "bc1_well") %>%
  mutate(
    celltype  = "Astrocytes",
    treatment = factor(treatment, levels = c("CTRL", "CORT", "DEX"))
  )

###############################################################
### F) SUMMARY (mean, sd, se, p-values, log2FC, labels)
###############################################################
df_sum_ast <- df_long_ast %>%
  group_by(celltype, treatment) %>%
  nest() %>%
  mutate(
    mean_CPM = map_dbl(data, ~ mean(.x$CPM, na.rm = TRUE)),
    sd_CPM   = map_dbl(data, ~ sd(.x$CPM, na.rm = TRUE)),
    n        = map_dbl(data, ~ sum(!is.na(.x$CPM))),
    se_CPM   = sd_CPM / sqrt(n),
    
    # p-values z edgeR
    pval = case_when(
      treatment == "CORT" ~ df_gene$PValue_edgerCORT,
      treatment == "DEX"  ~ df_gene$PValue_edgerDEX,
      TRUE ~ NA_real_
    ),
    
    # log2FC z edgeR
    logfc = case_when(
      treatment == "CORT" ~ df_gene$log2FC_edgerCORT,
      treatment == "DEX"  ~ df_gene$log2FC_edgerDEX,
      TRUE ~ NA_real_
    ),
    
    # format p-value
    p_label = case_when(
      is.na(pval) ~ "",
      pval < 1e-4 ~ paste0("p = ", format(pval, scientific = TRUE, digits = 2)),
      TRUE        ~ paste0("p = ", signif(pval, 2))
    ),
    
    # final label p + log2FC
    label_fc = case_when(
      is.na(logfc) ~ p_label,
      TRUE         ~ paste0(p_label, " ; log2FC = ", round(logfc, 2))
    )
  ) %>%
  select(-data)

###############################################################
### G) USTAWIENIE POZIOMEJ WYSOKOŚCI BRACKETÓW
###############################################################
y_max <- max(df_sum_ast$mean_CPM + df_sum_ast$se_CPM)

df_brackets <- tibble(
  xstart = c(1, 1),
  xend   = c(2, 3),
  y      = c(y_max * 1.15, y_max * 1.28),
  label  = c(
    df_sum_ast$label_fc[df_sum_ast$treatment == "CORT"],
    df_sum_ast$label_fc[df_sum_ast$treatment == "DEX"]
  )
)

###############################################################
### H) WYKRES FINALNY
###############################################################
ggplot() +
  
  # punkty
  geom_jitter(
    data = df_long_ast,
    aes(x = treatment, y = CPM, color = treatment),
    width = 0.15,
    size = 3,
    alpha = 0.85
  ) +
  
  # bary
  geom_col(
    data = df_sum_ast,
    aes(x = treatment, y = mean_CPM, fill = treatment),
    width = 0.6,
    alpha = 0.75
  ) +
  
  # errorbary
  geom_errorbar(
    data = df_sum_ast,
    aes(
      x = treatment,
      ymin = mean_CPM - se_CPM,
      ymax = mean_CPM + se_CPM
    ),
    width = 0.15,
    linewidth = 0.7,
    color = "black"
  ) +
  
  ###############################################################
# BRACKETY — LINIJKI
###############################################################
geom_segment(
  data = df_brackets,
  aes(x = xstart, xend = xend, y = y, yend = y),
  linewidth = 0.8
) +
  geom_segment(
    data = df_brackets,
    aes(x = xstart, xend = xstart, y = y - 0.02*y_max, yend = y),
    linewidth = 0.8
  ) +
  geom_segment(
    data = df_brackets,
    aes(x = xend, xend = xend, y = y - 0.02*y_max, yend = y),
    linewidth = 0.8
  ) +
  
  # podpisy
  geom_text(
    data = df_brackets,
    aes(x = (xstart + xend)/2, y = y + 0.03*y_max, label = label),
    size = 5,
    fontface = "bold"
  ) +
  
  ###############################################################
### STYLE
###############################################################
scale_fill_manual(values = c("CTRL" = "grey40", "CORT" = "indianred2", "DEX" = "steelblue3")) +
  scale_color_manual(values = c("CTRL" = "grey40", "CORT" = "indianred2", "DEX" = "steelblue3")) +
  theme_classic(base_size = 15) +
  labs(
    title = paste0(gene_symbol, " — Astrocytes"),
    x = "Treatment",
    y = "CPM"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
