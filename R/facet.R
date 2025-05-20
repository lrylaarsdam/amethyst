######################################################################
#' @title loadWindows
#' @description Aggregate c and t observations over smoothed genomic windows of a user-defined size.
#' Optionally calculate the % methylation for each window.
#'
#' @param obj Amethyst object containing the path(s) to calculate.
#' @param type Type of methylation - e.g. "CG" or "CH" - to calculate. Typically will be "CG."
#' @param threads Optional number of threads for parallelization. Basically required for this step.
#' @param index Name of chr index in the index slot.
#' This reduces memory constraints by processing one chromosome at a time.
#' @param step Width of the genomic window to calculate. Default is 500 bp.
#' @param smooth Number of windows to include surrounding the target region; i.e. produces a sliding window matrix.
#' Default parameter is 3, resulting in a 1500 x 500 bp sliding window matrix.
#' @param genome Genome build of the organism(s) being analyzed. Options are currently "hg19", "hg38", "mm10", or "mm39".
#' @param futureType If using parallelization, should R multithread using "multicore" or "multisession".
#' @param groupBy Parameter contained in the metadata over which to aggregate observations. Default is "cluster_id".
#' @param returnSumMatrix Whether or not the function should return the matrix of summed c and t observations. Required for testDMR input.
#' @param returnPctMatrix Whether or not the function should calculate % methylation over each genomic window. Required for heatMap input.
#' @return if returnSumMatrix = TRUE, returns a data.table with the genomic location as chr, start, and end columns,
#' plus aggregated c and t observations for each groupBy value. If returnPctMatrix = TRUE, returns a data.table with the
#' genomic location as chr, start, and end columns, plus one column of the mean % methylation for each window per groupBy value.
#' If both are true, returns a list containing each result.
#' @importFrom dplyr filter pull mutate select left_join
#' @importFrom stringr str_count
#' @importFrom rhdf5 h5read h5ls
#' @importFrom future plan multicore sequential
#' @importFrom furrr future_map future_pmap
#' @importFrom plyr round_any .
#' @importFrom tidyr separate
#' @importFrom data.table setnames setDT tstrsplit := setorder data.table frollapply frollsum frollmean shift copy
#' @importFrom purrr reduce
#' @importFrom utils read.table
#' @export
#'
#' @examples output <- calcSmoothedWindows(obj, returnSumMatrix = TRUE, returnPctMatrix = FALSE) # example to calculate testDMR input
#' @examples output <- calcSmoothedWindows(obj, returnSumMatrix = FALSE, returnPctMatrix = TRUE) # example to calculate heatMap input
loadSmoothedWindows <- function(
    obj,
    type = "CG",
    threads = 1,
    name = "500",
    smoothTrim = 500,
    futureType = "multicore",
    groupBy = "cluster_id",
    nmin = 10,
    chunkSize = 50,
    returnSumMatrix = TRUE,
    returnPctMatrix = TRUE,
    save = FALSE) {
  
  # get barcodes and paths from amethyst object
  if (is.null(obj@h5paths)) {
    stop("Please generate the path list for each barcode and store in the obj@h5paths slot.")
  } else {
    h5paths <- obj@h5paths
  }
  
  # set up multithreading
  if (threads > 1) {
    if (futureType == "multicore") {
      future::plan(future::multicore, workers = threads)
    }
    if (futureType == "multisession") {
      future::plan(future::multisession, workers = threads)
    } else if (!(futureType %in% c("multicore", "multisession"))) {
      stop("Options for parallelizing include multicore or multisession.")
    }
  }
  
  # get number of clusters
  membership <- obj@metadata |> dplyr::select(groupBy)
  colnames(membership) <- "membership"
  membership <- membership |> dplyr::filter(!is.na(membership))
  
  valid_groups <- membership |> 
    dplyr::count(membership) |> 
    dplyr::filter(n > nmin) |> 
    dplyr::pull(membership)
  
  groups <- as.list(unique(membership$membership[membership$membership %in% valid_groups]))
  group_results <- list()
  
  for (gr in groups) {
    members <- rownames(membership |> dplyr::filter(membership == gr))
    member_paths <- h5paths[match(members, rownames(h5paths)), "paths"]
    
    # Split members and paths into chunks
    n_chunks <- ceiling(length(members) / chunkSize)
    member_chunks <- split(seq_along(members), ceiling(seq_along(members) / chunkSize))
    
    group_result_chunks <- list()
    
    for (i in seq_along(member_chunks)) {
      idx <- member_chunks[[i]]
      chunk_members <- members[idx]
      chunk_paths <- member_paths[idx]
      
      chunk_results <- furrr::future_pmap(
        .l = list(chunk_paths, chunk_members),
        .f = function(path, barcode) {
          tryCatch({
            barcode_name <- sub("\\..*$", "", barcode)
            data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode, "/", name)))
          }, error = function(e) {
            cat("Error processing data for barcode", barcode, ":", conditionMessage(e), "\n")
            return(NULL)
          })
        },
        .progress = TRUE
      )
      
      # Filter out NULLs
      chunk_results <- purrr::compact(chunk_results)
      chunk_results <- data.table::rbindlist(chunk_results)
      
      if (nrow(chunk_results)) {
        chunk_results <- chunk_results[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = .(chr, start, end)]
        group_result_chunks[[i]] <- chunk_results
      }
    }
    
    # Combine all chunks for this group
    if (length(group_result_chunks)) {
      member_results <- data.table::rbindlist(group_result_chunks)
      member_results <- member_results[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = .(chr, start, end)]
      data.table::setnames(member_results, c("chr", "start", "end", paste0(gr, "_c"), paste0(gr, "_t")))
      group_results[[gr]] <- member_results
    } else {
      cat("No valid data found for group", gr, "\n")
    }
  }
  
  count_matrix <- Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = FALSE), group_results)
  data.table::setorder(count_matrix, chr, start)
  
  count_matrix <- count_matrix[, `:=`(start = (start + smoothTrim), end = (end - smoothTrim))]
  
  if (threads > 1) {
    future::plan(future::sequential)
    gc()
  }
  
  # Calculate cluster track from count_matrix
  if (returnPctMatrix) {
    pctm_matrix <- data.table::copy(count_matrix)
    for (gr in groups) {
      group_c_cols <- paste0(gr, "_c")
      group_t_cols <- paste0(gr, "_t")
      pctm_matrix[, m := round(.SD[[group_c_cols]] * 100 / (.SD[[group_c_cols]] + .SD[[group_t_cols]]), 2), .SDcols = c(group_c_cols, group_t_cols)]
      data.table::setnames(pctm_matrix, "m", as.character(gr))
    }
    
    pctm_matrix <- pctm_matrix[, c(paste0(groups, "_c"), paste0(groups, "_t")) := NULL]
  }
  
  if (returnSumMatrix & returnPctMatrix) {
    return(list(sum_matrix = count_matrix, pct_matrix = pctm_matrix))
  } else if (!returnSumMatrix & returnPctMatrix) {
    return(pctm_matrix)
  } else if (returnSumMatrix & ! returnPctMatrix) {
    return(count_matrix)
  }
}


######################################################################
#' @title loadSmoothedWindows
#' @description Aggregate c and t observations over smoothed genomic windows of a user-defined size.
#' Optionally calculate the % methylation for each window.
#'
#' @param obj Amethyst object containing the path(s) to calculate.
#' @param type Type of methylation - e.g. "CG" or "CH" - to calculate. Typically will be "CG."
#' @param threads Optional number of threads for parallelization. Basically required for this step.
#' @param index Name of chr index in the index slot.
#' This reduces memory constraints by processing one chromosome at a time.
#' @param step Width of the genomic window to calculate. Default is 500 bp.
#' @param smooth Number of windows to include surrounding the target region; i.e. produces a sliding window matrix.
#' Default parameter is 3, resulting in a 1500 x 500 bp sliding window matrix.
#' @param genome Genome build of the organism(s) being analyzed. Options are currently "hg19", "hg38", "mm10", or "mm39".
#' @param futureType If using parallelization, should R multithread using "multicore" or "multisession".
#' @param groupBy Parameter contained in the metadata over which to aggregate observations. Default is "cluster_id".
#' @param returnSumMatrix Whether or not the function should return the matrix of summed c and t observations. Required for testDMR input.
#' @param returnPctMatrix Whether or not the function should calculate % methylation over each genomic window. Required for heatMap input.
#' @return if returnSumMatrix = TRUE, returns a data.table with the genomic location as chr, start, and end columns,
#' plus aggregated c and t observations for each groupBy value. If returnPctMatrix = TRUE, returns a data.table with the
#' genomic location as chr, start, and end columns, plus one column of the mean % methylation for each window per groupBy value.
#' If both are true, returns a list containing each result.
#' @importFrom dplyr filter pull mutate select left_join
#' @importFrom stringr str_count
#' @importFrom rhdf5 h5read h5ls
#' @importFrom future plan multicore sequential
#' @importFrom furrr future_map future_pmap
#' @importFrom plyr round_any .
#' @importFrom tidyr separate
#' @importFrom data.table setnames setDT tstrsplit := setorder data.table frollapply frollsum frollmean shift copy
#' @importFrom purrr reduce
#' @importFrom utils read.table
#' @export
#'
#' @examples output <- calcSmoothedWindows(obj, returnSumMatrix = TRUE, returnPctMatrix = FALSE) # example to calculate testDMR input
#' @examples output <- calcSmoothedWindows(obj, returnSumMatrix = FALSE, returnPctMatrix = TRUE) # example to calculate heatMap input
loadSmoothedWindows <- function(
    obj,
    type = "CG",
    threads = 1,
    name = "500",
    smoothTrim = 500,
    futureType = "multicore",
    groupBy = "cluster_id",
    nmin = 10,
    chunkSize = 50,
    returnSumMatrix = TRUE,
    returnPctMatrix = TRUE,
    save = FALSE) {
  
  # get barcodes and paths from amethyst object
  if (is.null(obj@h5paths)) {
    stop("Please generate the path list for each barcode and store in the obj@h5paths slot.")
  } else {
    h5paths <- obj@h5paths
  }
  
  # set up multithreading
  if (threads > 1) {
    if (futureType == "multicore") {
      future::plan(future::multicore, workers = threads)
    }
    if (futureType == "multisession") {
      future::plan(future::multisession, workers = threads)
    } else if (!(futureType %in% c("multicore", "multisession"))) {
      stop("Options for parallelizing include multicore or multisession.")
    }
  }
  
  # get number of clusters
  membership <- obj@metadata |> dplyr::select(groupBy)
  colnames(membership) <- "membership"
  membership <- membership |> dplyr::filter(!is.na(membership))
  
  valid_groups <- membership |> 
    dplyr::count(membership) |> 
    dplyr::filter(n > nmin) |> 
    dplyr::pull(membership)
  
  groups <- as.list(unique(membership$membership[membership$membership %in% valid_groups]))
  group_results <- list()
  
  for (gr in groups) {
    members <- rownames(membership |> dplyr::filter(membership == gr))
    member_paths <- h5paths[match(members, rownames(h5paths)), "paths"]
    
    # Split members and paths into chunks
    n_chunks <- ceiling(length(members) / chunkSize)
    member_chunks <- split(seq_along(members), ceiling(seq_along(members) / chunkSize))
    
    group_result_chunks <- list()
    
    for (i in seq_along(member_chunks)) {
      idx <- member_chunks[[i]]
      chunk_members <- members[idx]
      chunk_paths <- member_paths[idx]
      
      chunk_results <- furrr::future_pmap(
        .l = list(chunk_paths, chunk_members),
        .f = function(path, barcode) {
          tryCatch({
            barcode_name <- sub("\\..*$", "", barcode)
            data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode, "/", name)))
          }, error = function(e) {
            cat("Error processing data for barcode", barcode, ":", conditionMessage(e), "\n")
            return(NULL)
          })
        },
        .progress = TRUE
      )
      
      # Filter out NULLs
      chunk_results <- purrr::compact(chunk_results)
      chunk_results <- data.table::rbindlist(chunk_results)
      
      if (nrow(chunk_results)) {
        chunk_results <- chunk_results[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = .(chr, start, end)]
        group_result_chunks[[i]] <- chunk_results
      }
    }
    
    # Combine all chunks for this group
    if (length(group_result_chunks)) {
      member_results <- data.table::rbindlist(group_result_chunks)
      member_results <- member_results[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = .(chr, start, end)]
      data.table::setnames(member_results, c("chr", "start", "end", paste0(gr, "_c"), paste0(gr, "_t")))
      group_results[[gr]] <- member_results
    } else {
      cat("No valid data found for group", gr, "\n")
    }
  }
  
  count_matrix <- Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = FALSE), group_results)
  data.table::setorder(count_matrix, chr, start)
  
  count_matrix <- count_matrix[, `:=`(start = (start + smoothTrim), end = (end - smoothTrim))]
  
  if (threads > 1) {
    future::plan(future::sequential)
    gc()
  }
  
  # Calculate cluster track from count_matrix
  if (returnPctMatrix) {
    pctm_matrix <- data.table::copy(count_matrix)
    for (gr in groups) {
      group_c_cols <- paste0(gr, "_c")
      group_t_cols <- paste0(gr, "_t")
      pctm_matrix[, m := round(.SD[[group_c_cols]] * 100 / (.SD[[group_c_cols]] + .SD[[group_t_cols]]), 2), .SDcols = c(group_c_cols, group_t_cols)]
      data.table::setnames(pctm_matrix, "m", as.character(gr))
    }
    
    pctm_matrix <- pctm_matrix[, c(paste0(groups, "_c"), paste0(groups, "_t")) := NULL]
  }
  
  if (returnSumMatrix & returnPctMatrix) {
    return(list(sum_matrix = count_matrix, pct_matrix = pctm_matrix))
  } else if (!returnSumMatrix & returnPctMatrix) {
    return(pctm_matrix)
  } else if (returnSumMatrix & ! returnPctMatrix) {
    return(count_matrix)
  }
}


