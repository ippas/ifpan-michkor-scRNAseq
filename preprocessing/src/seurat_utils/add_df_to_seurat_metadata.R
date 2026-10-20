add_df_to_seurat_metadata <- function(
    seurat_obj,
    df,
    join_by_meta,           # column name in Seurat metadata to join by
    join_by_df = "cluster_id",  # column name in df to join by (default: "cluster_id")
    exclude_columns = NULL       # optional: vector of columns to exclude from df
){
  # --- Input checks ---
  stopifnot(join_by_df %in% colnames(df))  # ensure join_by_df exists in df
  
  # --- Ensure join column is character to avoid join issues ---
  df <- df %>% mutate("{join_by_df}" := as.character(.data[[join_by_df]]))
  
  # --- Optionally drop specified columns from df before joining ---
  if (!is.null(exclude_columns)) {
    df <- df %>% select(-any_of(exclude_columns))
  }
  
  # --- Prepare metadata to add ---
  # 1. Extract Seurat metadata and keep cell IDs as a column
  # 2. Ensure join_by_meta is character to match df
  # 3. Left join df by matching join_by_df (from df) to join_by_meta (from metadata)
  # 4. Reorder rows to match original Seurat object cell order
  # 5. Keep only columns from df (excluding join key)
  meta_to_add <- seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate("{join_by_meta}" := as.character(.data[[join_by_meta]])) %>%
    left_join(df, by = setNames(join_by_df, join_by_meta)) %>%
    arrange(match(cell_id, rownames(seurat_obj@meta.data))) %>%
    select(all_of(setdiff(colnames(df), join_by_df))) %>%
    as.data.frame()
  
  # --- Add new metadata columns back to the Seurat object ---
  seurat_obj <- AddMetaData(seurat_obj, metadata = meta_to_add)
  
  # --- Return updated Seurat object ---
  return(seurat_obj)
}
