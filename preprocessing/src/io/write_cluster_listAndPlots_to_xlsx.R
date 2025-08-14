write_cluster_listAndPlots_to_xlsx <- function(
    cluster_list,               # lista data.frame'ów per klaster
    plot_list,                  # lista wykresów (np. top9_markersUp_plots_list)
    file_path,                  # ścieżka do XLSX
    img_width = 18,              # szerokość obrazka (cale)
    img_height = 15,             # wysokość obrazka (cale)
    img_dpi = 300,              # DPI obrazka
    img_offset_cols = 3,        # ile kolumn odstępu od tabeli
    start_row = 2,              # wiersz startowy dla obrazka
    auto_colwidth = TRUE        # automatyczna szerokość kolumn
) {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  
  sorted_names <- names(cluster_list)[order(as.integer(gsub(".*_(\\d+)$", "\\1", names(cluster_list))))]
  
  wb <- openxlsx::createWorkbook()
  
  # utwórz progress bar
  pb <- txtProgressBar(min = 0, max = length(sorted_names), style = 3)
  
  for (i in seq_along(sorted_names)) {
    nm <- sorted_names[i]
    df <- cluster_list[[nm]]
    
    safe_sheet <- substr(gsub("[:\\\\/\\?\\*\\[\\]]", "_", nm), 1, 31)
    
    openxlsx::addWorksheet(wb, safe_sheet)
    openxlsx::writeData(wb, sheet = safe_sheet, x = df)
    openxlsx::freezePane(wb, sheet = safe_sheet, firstRow = TRUE)
    
    if (isTRUE(auto_colwidth)) {
      openxlsx::setColWidths(wb, sheet = safe_sheet, cols = 1:ncol(df), widths = "auto")
    }
    
    if (!is.null(plot_list) && nm %in% names(plot_list) && !is.null(plot_list[[nm]])) {
      tmp_png <- tempfile(fileext = ".png")
      ggplot2::ggsave(
        filename = tmp_png, plot = plot_list[[nm]],
        width = img_width, height = img_height, dpi = img_dpi, units = "in"
      )
      
      start_col <- ncol(df) + img_offset_cols
      
      openxlsx::insertImage(
        wb, sheet = safe_sheet, file = tmp_png,
        startRow = start_row, startCol = start_col,
        width = img_width, height = img_height
      )
    }
    
    setTxtProgressBar(pb, i)  # aktualizacja progress baru
  }
  
  close(pb)  # zamknięcie progress baru
  
  openxlsx::saveWorkbook(wb, file = file_path, overwrite = TRUE)
  invisible(TRUE)
}
