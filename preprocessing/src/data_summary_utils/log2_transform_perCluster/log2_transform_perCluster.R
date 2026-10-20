log2_transform_perCluster <- function(
    count_list,
    pseudocount = 1,
    progress = TRUE
) {
  if (!is.list(count_list)) {
    stop("Input must be a list of matrices (e.g., from normalize_CPM_perCluster())")
  }
  
  out_list <- vector("list", length(count_list))
  names(out_list) <- names(count_list)
  
  if (progress) pb <- utils::txtProgressBar(min = 0, max = length(count_list), style = 3)
  
  for (i in seq_along(count_list)) {
    m <- count_list[[i]]
    if (!is.matrix(m) && !inherits(m, "Matrix")) {
      m <- as.matrix(m)
    }
    out_list[[i]] <- log2(m + pseudocount)
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  
  if (progress) close(pb)
  
  return(out_list)
}