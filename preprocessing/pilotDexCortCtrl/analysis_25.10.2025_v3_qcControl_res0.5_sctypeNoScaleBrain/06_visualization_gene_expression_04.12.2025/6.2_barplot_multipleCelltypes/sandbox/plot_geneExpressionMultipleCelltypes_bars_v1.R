library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)

# ============================================================
# 1. DF CARMIL1 — z map(~ ...)
# ============================================================
df_carmil1 <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain %>% 
  purrr::map(~ dplyr::filter(.x, gene_symbol %in% res_heatmap$cluster_genes$Cluster_1)) %>% 
  .[celltypes] %>% 
  dplyr::bind_rows(.id = "celltype") %>% 
  dplyr::filter(gene_symbol == "CARMIL1")

# ============================================================
# 2. SAMPLE MAP
# ============================================================
sample_map <- integrated_pilotDexCortCtrl_filtered@meta.data %>% 
  dplyr::select(bc1_well, treatment) %>% 
  distinct()

# ============================================================
# 3. LONG FORMAT
# ============================================================
df_long <- df_carmil1 %>%
  dplyr::select(celltype, gene_symbol, starts_with("A")) %>%
  tidyr::pivot_longer(
    cols = starts_with("A"),
    names_to = "bc1_well",
    values_to = "CPM"
  ) %>%
  dplyr::left_join(sample_map, by = "bc1_well")

# ============================================================
# 4. NAPRAWA NAZW CELLTYPE
# ============================================================
df_long$celltype <- df_long$celltype %>%
  str_trim() %>% 
  str_replace("Oligodendrocyte precursor cells", "OPC")

df_carmil1$celltype <- df_carmil1$celltype %>%
  str_trim() %>% 
  str_replace("Oligodendrocyte precursor cells", "OPC")

# ============================================================
# 5. LEVELS DLA CELL TYPES — TAK JAK CHCIAŁEŚ
# ============================================================
celltype_levels <- c(
  "Radial glial cells",
  "Astrocytes",
  "OPC",
  "Oligodendrocytes",
  "Immature neurons",
  "Mature neurons"
)

df_long$celltype <- factor(df_long$celltype, levels = celltype_levels)

# ============================================================
# 6. ORDER TREATMENT — CTRL → CORT → DEX
# ============================================================
df_long$treatment <- factor(df_long$treatment, levels = c("CTRL", "CORT", "DEX"))

# ============================================================
# 7. SUMMARY (MEAN + SE) — z nest + map
# ============================================================
df_summary <- df_long %>%
  dplyr::group_by(celltype, treatment) %>%
  tidyr::nest() %>%
  mutate(
    mean_CPM = map_dbl(data, ~ mean(.x$CPM, na.rm = TRUE)),
    sd_CPM   = map_dbl(data, ~ sd(.x$CPM, na.rm = TRUE)),
    n        = map_dbl(data, ~ sum(!is.na(.x$CPM))),
    se_CPM   = sd_CPM / sqrt(n)
  ) %>%
  dplyr::select(-data)

# ============================================================
# 8. WYKRES (JITTER + MEAN + SE)
# ============================================================
ggplot() +
  
  # RAW POINTS
  geom_jitter(
    data = df_long,
    aes(x = celltype, y = CPM, color = treatment),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.75
    ),
    size = 2,
    alpha = 0.85
  ) +
  
  # MEAN BARS
  geom_col(
    data = df_summary,
    aes(x = celltype, y = mean_CPM, fill = treatment),
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 0.75
  ) +
  
  # SE ERROR BARS
  geom_errorbar(
    data = df_summary,
    aes(
      x = celltype,
      ymin = mean_CPM - se_CPM,
      ymax = mean_CPM + se_CPM,
      group = treatment
    ),
    position = position_dodge(width = 0.75),
    width = 0.2,
    color = "black",
    linewidth = 0.7
  ) +
  
  # WYGLĄD
  theme_classic(base_size = 15) +
  labs(
    title = "CARMIL1 — CPM per treatment across cell types",
    x = "Cell type",
    y = "CPM"
  ) +
  scale_fill_manual(values = c(
    "CTRL" = "grey40",
    "CORT" = "firebrick2",
    "DEX"  = "steelblue3"
  )) +
  scale_color_manual(values = c(
    "CTRL" = "grey40",
    "CORT" = "firebrick2",
    "DEX"  = "steelblue3"
  )) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )
