###############################################################
### 0) PACKAGES
###############################################################
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)

###############################################################
### 1) PARAMETRY
###############################################################
gene_symbol <- "FKBP5"

###############################################################
### 2) WYCIĄGNIĘCIE GENU DLA ASTROCYTÓW
###############################################################
df_gene <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$Astrocytes %>%
  filter(gene_symbol == !!gene_symbol)

if(nrow(df_gene) == 0){
  stop(paste("Brak genu:", gene_symbol))
}

###############################################################
### 3) SAMPLE MAP (CTRL / CORT / DEX)
###############################################################
sample_map <- integrated_pilotDexCortCtrl_filtered@meta.data %>%
  select(bc1_well, treatment) %>%
  distinct()

###############################################################
### 4) LONG FORMAT
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
    celltype = "Astrocytes",
    treatment = factor(treatment, levels = c("CTRL", "CORT", "DEX"))
  )

###############################################################
### 5) SUMMARY (nest → mean, SD, SE, p, log2FC)
###############################################################
df_sum_ast <- df_long_ast %>%
  group_by(celltype, treatment) %>%
  nest() %>%
  mutate(
    mean_CPM = map_dbl(data, ~ mean(.x$CPM, na.rm = TRUE)),
    sd_CPM   = map_dbl(data, ~ sd(.x$CPM, na.rm = TRUE)),
    n        = map_dbl(data, ~ sum(!is.na(.x$CPM))),
    se_CPM   = sd_CPM / sqrt(n),
    
    ### p-value edgeR
    pval = case_when(
      treatment == "CORT" ~ df_gene$PValue_edgerCORT,
      treatment == "DEX"  ~ df_gene$PValue_edgerDEX,
      TRUE ~ NA_real_
    ),
    
    ### log2FC edgeR
    logfc_edger = case_when(
      treatment == "CORT" ~ df_gene$log2FC_edgerCORT,
      treatment == "DEX"  ~ df_gene$log2FC_edgerDEX,
      TRUE ~ NA_real_
    ),
    
    ### label p-value
    p_label = case_when(
      is.na(pval) ~ "",
      pval < 1e-4 ~ paste0("p = ", format(pval, scientific = TRUE, digits = 2)),
      TRUE ~ paste0("p = ", signif(pval, 2))
    ),
    
    ### FINAL LABEL: p + log2FC
    label_fc = ifelse(
      is.na(logfc_edger),
      p_label,
      paste0(p_label, " ; log2FC = ", round(logfc_edger, 2))
    )
  ) %>%
  select(-data)

###############################################################
### 6) GLOBAL MAX (do ustawienia wysokości tekstów)
###############################################################
global_max <- max(df_sum_ast$mean_CPM + df_sum_ast$se_CPM, na.rm = TRUE)

###############################################################
### 7) WYKRES
###############################################################
ggplot() +
  
  ### RAW POINTS
  geom_jitter(
    data = df_long_ast,
    aes(x = treatment, y = CPM, color = treatment),
    width = 0.15,
    size = 3,
    alpha = 0.85
  ) +
  
  ### MEAN BARS
  geom_col(
    data = df_sum_ast,
    aes(x = treatment, y = mean_CPM, fill = treatment),
    width = 0.6,
    alpha = 0.75
  ) +
  
  ### SE ERROR BARS
  geom_errorbar(
    data = df_sum_ast,
    aes(
      x = treatment,
      ymin = mean_CPM - se_CPM,
      ymax = mean_CPM + se_CPM
    ),
    width = 0.15,
    color = "black",
    linewidth = 0.7
  ) +
  
  ### LABEL: p + log2FC
  geom_text(
    data = df_sum_ast %>% filter(treatment != "CTRL"),
    aes(
      x = treatment,
      y = mean_CPM + se_CPM + 0.07 * global_max,
      label = label_fc
    ),
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  
  ### STYL
  scale_fill_manual(values = c("CTRL" = "grey40", "CORT" = "firebrick2", "DEX" = "steelblue3")) +
  scale_color_manual(values = c("CTRL" = "grey40", "CORT" = "firebrick2", "DEX" = "steelblue3")) +
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
