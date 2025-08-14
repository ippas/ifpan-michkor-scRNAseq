write_cluster_listAndTwoPlots_to_xlsx <- function(
    cluster_list,                   # lista data.frame'ów per klaster
    plot_list_a,                    # lista 1: wykresy per klaster
    plot_list_b,                    # lista 2: wykresy per klaster
    file_path,                      # ścieżka do XLSX
    img_width_a = 18,               # szerokość 1. obrazka (cale)
    img_height_a = 15,              # wysokość 1. obrazka (cale)
    img_width_b = 18,               # szerokość 2. obrazka (cale)
    img_height_b = 15,              # wysokość 2. obrazka (cale)
    img_dpi = 300,                   # DPI dla obu
    img_offset_cols = 3,            # ile kolumn odstępu od tabeli
    start_row_first = 2,            # wiersz startowy dla pierwszego obrazka
    gap_rows = 25,                   # odstęp (wiersze) między pierwszym a drugim
    auto_colwidth = TRUE,            # auto szerokość kolumn
    round_digits_fc = 5,             # liczba miejsc po przecinku dla avg_log2FC
    pval_digits = 5                  # liczba miejsc po przecinku dla p_val i p_val_adj (scientific)
) {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  
  sorted_names <- names(cluster_list)[order(as.integer(gsub(".*_(\\d+)$", "\\1", names(cluster_list))))]
  wb <- openxlsx::createWorkbook()
  pb <- txtProgressBar(min = 0, max = length(sorted_names), style = 3)
  
  for (i in seq_along(sorted_names)) {
    nm <- sorted_names[i]
    df <- cluster_list[[nm]]
    
    # --- Zaokrąglanie i formatowanie kolumn ---
    if ("avg_log2FC" %in% names(df)) {
      df$avg_log2FC <- round(df$avg_log2FC, round_digits_fc)
    }
    if ("p_val" %in% names(df)) {
      df$p_val <- formatC(df$p_val, format = "e", digits = pval_digits)
    }
    if ("p_val_adj" %in% names(df)) {
      df$p_val_adj <- formatC(df$p_val_adj, format = "e", digits = pval_digits)
    }
    
    safe_sheet <- substr(gsub("[:\\\\/\\?\\*\\[\\]]", "_", nm), 1, 31)
    openxlsx::addWorksheet(wb, safe_sheet)
    openxlsx::writeData(wb, sheet = safe_sheet, x = df)
    openxlsx::freezePane(wb, sheet = safe_sheet, firstRow = TRUE)
    
    if (isTRUE(auto_colwidth)) {
      openxlsx::setColWidths(wb, sheet = safe_sheet, cols = 1:ncol(df), widths = "auto")
    }
    
    start_col <- ncol(df) + img_offset_cols
    
    # --- Pierwszy obrazek ---
    if (!is.null(plot_list_a[[nm]])) {
      tmp_png_a <- tempfile(fileext = ".png")
      ggplot2::ggsave(
        filename = tmp_png_a,
        plot = plot_list_a[[nm]],
        width = img_width_a, height = img_height_a, dpi = img_dpi, units = "in"
      )
      openxlsx::insertImage(
        wb, sheet = safe_sheet, file = tmp_png_a,
        startRow = start_row_first, startCol = start_col,
        width = img_width_a, height = img_height_a
      )
    }
    
    # --- Drugi obrazek ---
    if (!is.null(plot_list_b[[nm]])) {
      tmp_png_b <- tempfile(fileext = ".png")
      ggplot2::ggsave(
        filename = tmp_png_b,
        plot = plot_list_b[[nm]],
        width = img_width_b, height = img_height_b, dpi = img_dpi, units = "in"
      )
      start_row_second <- start_row_first + gap_rows
      openxlsx::insertImage(
        wb, sheet = safe_sheet, file = tmp_png_b,
        startRow = start_row_second, startCol = start_col,
        width = img_width_b, height = img_height_b
      )
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  openxlsx::saveWorkbook(wb, file = file_path, overwrite = TRUE)
  invisible(TRUE)
}
