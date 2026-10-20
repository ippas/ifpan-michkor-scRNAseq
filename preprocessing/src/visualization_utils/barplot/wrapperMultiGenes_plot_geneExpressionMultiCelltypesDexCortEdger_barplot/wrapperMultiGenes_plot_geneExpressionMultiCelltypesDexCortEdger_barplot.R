wrapperMultiGenes_plot_geneExpressionMultiCelltypesDexCortEdger_barplot <- function(
    genes,
    output_dir,
    combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
    integrated_pilotDexCortCtrl_filtered,
    force_create_output_dir = FALSE,
    svg_width  = 12,
    svg_height = 6,
    ...   # dodatkowe parametry przekazywane do funkcji bazowej
) {
  # -------------------------------------------------------------------
  # A) Pakiety
  # -------------------------------------------------------------------
  library(ggplot2)
  library(fs)
  
  # -------------------------------------------------------------------
  # B) Sprawdzenie katalogu wyjściowego
  # -------------------------------------------------------------------
  if (!dir_exists(output_dir)) {
    if (force_create_output_dir) {
      dir_create(output_dir, recurse = TRUE)
      message("📁 Utworzono katalog: ", output_dir)
    } else {
      stop("Katalog nie istnieje: ", output_dir,
           "\nUżyj force_create_output_dir = TRUE aby go utworzyć.")
    }
  }
  
  # -------------------------------------------------------------------
  # C) Inicjalizacja listy wyjściowej
  # -------------------------------------------------------------------
  output_files <- list()
  
  # -------------------------------------------------------------------
  # D) Pętla po genach
  # -------------------------------------------------------------------
  for (gene in genes) {
    message("▶️ Generuję wykres dla genu: ", gene)
    
    # --- Wywołanie Twojej funkcji bazowej ---
    p <- plot_geneExpressionMultiCelltypesDexCortEdger_barplot(
      gene_symbol = gene,
      combined_summaryEdgerChi2_res0.3RNAnoScaledBrain = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
      integrated_pilotDexCortCtrl_filtered = integrated_pilotDexCortCtrl_filtered,
      ...
    )
    
    # --- Nazwa pliku ---
    file_path <- file.path(output_dir, paste0(gene, ".svg"))
    
    # --- Zapis do SVG ---
    svg(filename = file_path, width = svg_width, height = svg_height)
    print(p)
    dev.off()
    
    output_files[[gene]] <- file_path
    message("💾 Zapisano: ", file_path)
  }
  
  message("\n✅ Wszystkie wykresy wygenerowane.")
  
  return(output_files)
}



# genes_to_plot <- c("FKBP5", "HTRA1", "ZBTB16")
# 
# wrapperMultiGenes_plot_geneExpressionMultiCelltypesDexCortEdger_barplot(
#   genes = genes_to_plot,
#   output_dir = "results/figures/barplots_07.12.2025/test/",
#   combined_summaryEdgerChi2_res0.3RNAnoScaledBrain = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
#   integrated_pilotDexCortCtrl_filtered = integrated_pilotDexCortCtrl_filtered,
#   force_create_output_dir = TRUE,
#   svg_width  = 14,
#   svg_height = 7,
#   # poniżej możesz przekazać dowolne argumenty do funkcji głównej:
#     celltypes_order = c(
#       "Radial glial cells",
#       "Astrocytes",
#       "OPCs",
#       "Oligodendrocytes",
#       "Immature neurons",
#       "Mature neurons"
#     ),
#     bar_alpha         = 0.5,
#     jitter_alpha      = 0.85,
#     jitter_size       = 2.5,
#     text_size         = 10,          # base_size ggplot
#     bracket_text_size = 10 / 2.83,   # ~px
#     axis_text_size    = 12,
#     axis_title_size   = 12,
#     strip_text_size   = 12,
#     title_size        = 16
# )
