read_xlsx_allSheets <- function(
    file,
    header = TRUE,
    skip_dict = NULL,
    include_sheets = NULL,
    exclude_sheets = NULL
) {
  # file: path to .xlsx file
  # header: whether the first row should be used as column names (TRUE/FALSE)
  # skip_dict: named list; key = sheet name or index, value = number of rows to skip
  # include_sheets: vector of sheet names or indices to include (default: all)
  # exclude_sheets: vector of sheet names or indices to exclude (default: none)
  
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required. Install it with: install.packages('readxl')")
  }
  
  # Get all sheet names
  sheet_names <- readxl::excel_sheets(file)
  
  # If include_sheets is specified
  if (!is.null(include_sheets)) {
    sheet_names <- sheet_names[include_sheets %in% sheet_names | 
                                 suppressWarnings(as.numeric(include_sheets)) %in% seq_along(sheet_names)]
  }
  
  # If exclude_sheets is specified
  if (!is.null(exclude_sheets)) {
    sheet_names <- setdiff(
      sheet_names,
      c(
        sheet_names[exclude_sheets %in% sheet_names],
        sheet_names[suppressWarnings(as.numeric(exclude_sheets))]
      )
    )
  }
  
  # Read sheets
  sheet_list <- lapply(seq_along(sheet_names), function(i) {
    sheet_name <- sheet_names[i]
    
    # Default skip = 0
    skip_val <- 0
    
    # Check skip_dict (by name or index)
    if (!is.null(skip_dict)) {
      if (!is.null(skip_dict[[sheet_name]])) {
        skip_val <- skip_dict[[sheet_name]]
      } else if (!is.null(skip_dict[[as.character(i)]])) {
        skip_val <- skip_dict[[as.character(i)]]
      }
    }
    
    message(sprintf("Reading sheet '%s' (skip = %d)...", sheet_name, skip_val))
    
    df <- readxl::read_excel(
      path = file,
      sheet = sheet_name,
      col_names = header,
      skip = skip_val
    )
    
    return(df)
  })
  
  names(sheet_list) <- sheet_names
  return(sheet_list)
}
