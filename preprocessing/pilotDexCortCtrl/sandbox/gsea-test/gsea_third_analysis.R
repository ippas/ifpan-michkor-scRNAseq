
pmid33510160_geneMarkers <- read_xlsx_allSheets(file = "data/celltype_references/PMID33510160/41467_2021_20892_MOESM6_ESM.xlsx",
                    header = TRUE,
                    skip_dict = list("All differentially expressed" = 1))


pmid33510160_geneMarkers$OPC2$gene



# Step 1: Prepare average expression matrix (protein-coding v114, per cluster)
avg_mat <- prepare_avg_expr_by_sample_cluster(
  seurat_obj  = integrated_pilotDexCortCtrl,
  cluster_col = "seurat_clusters_res0.13",
  assay       = "RNA",
  # slot        = "counts",       # or "data" depending on what you want
  biomart_version = "112",
  select_samples =c(
    "BG60BLUE_A1", "BG60BLUE_A2", "BG60BLUE_A3", "BG60BLUE_A4",
    "BG60GREEN_A1", "BG60GREEN_A2", "BG60GREEN_A3", "BG60GREEN_A4"
  )
)


fkbp5_all
fkbp5_ctrl


# Załóżmy, że obie zmienne są wektorami numeric
fkbp5_all  <- as.numeric(fkbp5_all)
fkbp5_ctrl <- as.numeric(fkbp5_ctrl)

# ramka do ggplot
df <- data.frame(
  cluster = paste0("cluster_", 0:9),
  all     = fkbp5_all,
  ctrl    = fkbp5_ctrl
)

# korelacja
cor(df$all, df$ctrl, method = "pearson")

# scatterplot
library(ggplot2)
ggplot(df, aes(x = all, y = ctrl, label = cluster)) +
  geom_point(color = "blue", size = 3) +
  geom_text(vjust = -0.8, size = 3) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
  labs(x = "FKBP5 (all samples)", y = "FKBP5 (CTRL only)",
       title = "Comparison of FKBP5 expression across clusters") +
  theme_minimal()


# ramka z wartościami
df <- data.frame(
  cluster = paste0("cluster_", 0:9),
  all     = as.numeric(fkbp5_all),
  ctrl    = as.numeric(fkbp5_ctrl)
)

# --- normalizacja min-max w obrębie każdego warunku ---
df <- df %>%
  mutate(
    all_norm  = (all  - min(all))  / (max(all)  - min(all)),
    ctrl_norm = (ctrl - min(ctrl)) / (max(ctrl) - min(ctrl)),
    diff_norm = all_norm - ctrl_norm   # różnica znormalizowana
  )

# --- klastry z największą różnicą ---
df %>% arrange(desc(abs(diff_norm)))

# Check dimensions
dim(avg_mat)
# rows = protein-coding genes, cols = clusters

# Step 2: Compute specificity scores
spec <- compute_gene_specificity_score(
  expr_mat   = avg_mat,
  tau_variant= "mean",          # or "yanai"
  log1p      = F,           # set TRUE if you want log-transform
  return_zero_rows = T      # remove all-zero genes
)


fgsea_out <- run_fgsea_on_score_matrix(
  # mat        =  tmp %>% na.omit(),
  mat = spec$combined,
  gene_lists = panglaoDB_list_allMarkers[c("Pericytes",   "Astrocytes", "Bergmann_glia", "Microglia",
                                             "Oligodendrocyte_progenitor_cells", "Radial_glia_cells",
                                             "Satellite_glial_cells", "Tanycytes", "Pericytes",
                                             "Neurons", "Immature_neurons", "Neural_stem/precursor_cells")],   # or your custom marker sets
  nperm      = 1000,
  do_fisher = T
)


spec$combined
fgsea_out$first_per_cluster %>% arrange(cluster)


fgsea_out$results$cluster_4
