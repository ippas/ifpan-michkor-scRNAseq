# ============================================================
# 📘 Correlation analysis of logFC values between clusters
# ============================================================

# 🔹 Required packages
library(purrr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# 🧩 Helper function: compute correlation & draw heatmap
# ============================================================

plot_cluster_correlation <- function(results_list, contrast_name, title_suffix) {
  # 1️⃣ Extract logFC tables for the given contrast
  logfc_list <- map(
    results_list,
    ~ .x[[contrast_name]] %>% select(name_id, logFC)
  )
  names(logfc_list) <- names(results_list)
  
  # 2️⃣ Merge all data frames by 'name_id'
  logfc_df <- purrr::reduce(logfc_list, full_join, by = "name_id")
  
  # 3️⃣ Rename columns: first = name_id, others = cluster names
  colnames(logfc_df) <- c("name_id", names(logfc_list))
  
  # 4️⃣ Compute correlation matrix (excluding 'name_id')
  cor_matrix <- logfc_df %>%
    select(-name_id) %>%
    cor(use = "pairwise.complete.obs", method = "pearson")
  
  # 5️⃣ Define color scale (white → red → gray)
  col_fun <- circlize::colorRamp2(
    breaks = c(0, 0.9, 1),
    colors = c("white", "red", "grey50")
  )
  
  # 6️⃣ Draw the heatmap with text labels
  ComplexHeatmap::Heatmap(
    cor_matrix,
    name = "Pearson r",
    col = col_fun,
    column_title = paste0("Cluster correlation (", title_suffix, ")"),
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    rect_gp = gpar(col = "grey60"),
    heatmap_legend_param = list(
      title = "Pearson r",
      at = c(0, 0.5, 0.9, 1),
      labels = c("0", "0.5", "0.9", "1")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(
        sprintf("%.2f", cor_matrix[i, j]),  # correlation value
        x, y,
        gp = gpar(fontsize = 9, col = "black")
      )
    }
  )
}

# ============================================================
# 🧠 Input data
# ============================================================

# Assumes the object:
# edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results
# contains lists with elements like: $CORT_vs_CTRL and $DEX_vs_CTRL

# ============================================================
# 🔹 Correlation heatmap for CORT_vs_CTRL
# ============================================================

heatmap_cort <- plot_cluster_correlation(
  results_list = edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results,
  contrast_name = "CORT_vs_CTRL",
  title_suffix = "CORT_vs_CTRL"
)

# ============================================================
# 🔹 Correlation heatmap for DEX_vs_CTRL
# ============================================================

heatmap_dex <- plot_cluster_correlation(
  results_list = edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results,
  contrast_name = "DEX_vs_CTRL",
  title_suffix = "DEX_vs_CTRL"
)

# ============================================================
# 🎨 Display both heatmaps
# ============================================================

ComplexHeatmap::draw(heatmap_cort + heatmap_dex)
