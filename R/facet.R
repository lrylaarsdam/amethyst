# Authors: Lauren Rylaarsdam, PhD; Ben Skubi (Facet)
# 2024-2025

######################################################################
#' @title loadWindows
#' @description Load pre-computed c and t observations aggregated with Facet
#'
#' @param obj Amethyst object containing the h5paths to the pre-aggregated windows
#' @param type Type of methylation - e.g. "CG" or "CH" - to calculate. Typically will be "CG"
#' @param threads Optional number of threads for parallelization
#' @param metric Calculate either methylation "percent" (out of 100), "score" (normalized metric from -1 to 1 with
#' 0 as the mean), or "ratio" (value / global methylation)
#' @param name Name of aggregated dataset in h5 file
#' @param nmin Minimum number of observations to include window value
#' @param chunkSize How many cells to process at a time (for memory management purposes)
#' @param futureType If using parallelization, should R multithread using "multicore" or "multisession".
#' @return Returns a data frame with columns as cells, rows as genomic windows, and values as aggregated methylation levels
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
#' @examples
#' \dontrun{
#'   obj@genomeMatrices[["cg_100k_score_facet"]] <- loadWindows(obj, name = "100000", nmin = 2, metric = "score")
#' }
loadWindows <- function(
    obj,
    type = "CG",
    metric = "percent",
    threads = 1,
    name,
    futureType = "multicore",
    nmin = 10,
    chunkSize = 50) {

  # get barcodes and paths from amethyst object

  if (is.null(obj@h5paths$path) | is.null(obj@h5paths$barcode)) {
    stop("\nPlease make sure the obj@h5paths slot is constructed correctly. Path and/or barcode columns not found.")
  } else {
    barcodes <- as.list(obj@h5paths$barcode)
    paths <- as.list(obj@h5paths$path)
    if (!is.null(obj@h5paths$prefix)) {
      prefixes <- obj@h5paths$prefix
    } else {
      prefixes <- as.list(rep("", length(barcodes)))
    }
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

  # Split into chunks
  n_chunks <- ceiling(length(barcodes) / chunkSize)
  barcode_chunks <- split(seq_along(barcodes), ceiling(seq_along(barcodes) / chunkSize))

  results <- list()

  for (i in seq_along(barcode_chunks)) {
    idx <- barcode_chunks[[i]]
    chunk_members <- barcodes[idx]
    chunk_paths <- paths[idx]
    chunk_prefixes <- prefixes[idx]

    chunk_results <- furrr::future_pmap(
      .l = list(chunk_paths, chunk_members, chunk_prefixes),
      .f = function(path, barcode, prefix) {
        tryCatch({

          # read in data
          barcode_name <- paste0(prefix, sub("\\..*$", "", barcode))
          h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode, "/", name)))
          h5 <- unique(h5) # added 250922 to filter duplicate rows

          if (metric != "percent") {
            meth_cell <- obj@metadata[barcode_name, paste0("m", tolower(type), "_pct")]/100 # pull global methylation level from metadata
          }

          meth_window <- h5[, .(value = round(sum(c_nz) / (sum(c_nz) + sum(t_nz)), 3), n = sum(c_nz + t_nz, na.rm = TRUE)), by = c("chr", "start", "end")]
          meth_window <- meth_window[n >= nmin, .(chr, start, end, value)]

          if (metric == "percent") { summary <- meth_window[, value := (value * 100)] }
          if (metric == "score") { summary <- meth_window[, value := round((ifelse((value - meth_cell > 0),
                                                                                   (value - meth_cell)/(1 - meth_cell),
                                                                                   (value - meth_cell)/meth_cell)), 3)] }
          if (metric == "ratio") { summary <- meth_window[, value := round(value / meth_cell, 3)] }
          data.table::setnames(summary, c("chr", "start", "end", barcode_name))

        }, error = function(e) {
          cat("Error processing data for barcode", barcode, ":", conditionMessage(e), "\n")
          return(NULL)
        })
      },
      .progress = TRUE
    )

    # 251015 - filter step added to save join when no data is present for a cell 251015
    chunk_results <- Filter(function(x) !is.null(x) && nrow(x) > 0, chunk_results)

    results[[i]] <- Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = FALSE), chunk_results)
    data.table::setorder(results[[i]], chr, start)

  }

 results <- Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = FALSE), results)
 data.table::setorder(results, chr, start)

 results <- results[, window := paste0(chr, "_", start, "_", end)]
 results <- results[, c("chr", "start", "end") := NULL]
 results <- results |> tibble::column_to_rownames(var = "window")

}


######################################################################
#' @title loadSmoothedWindows
#' @description Aggregate c and t observations over smoothed genomic windows of a user-defined size.
#' Optionally calculate the % methylation for each window.
#'
#' @param obj Amethyst object containing the h5paths to the pre-aggregated windows
#' @param type Type of methylation - e.g. "CG" or "CH" - to calculate. Typically will be "CG"
#' @param threads Optional number of threads for parallelization (multithreading strongly suggested)
#' @param futureType If using parallelization, should R multithread using "multicore" or "multisession"
#' @param groupBy Parameter contained in the metadata over which to aggregate observations. Default is "cluster_id"
#' @param returnSumMatrix Whether or not the function should return the matrix of summed c and t observations. Required for testDMR input
#' @param name Name of aggregated dataset in h5 file
#' @param smoothTrim Width of values to trim for output in the "start" and "end" columns. For example, we typically
#' calculate a 1500x500bp smoothed sliding window matrix, so 500 is trimmed off either end to produce smoothed unique values
#' @param nmin Minimum number of observations to include window value
#' @param chunkSize How many cells to process at a time (for memory management purposes)
#' @param returnPctMatrix Whether or not the function should calculate % methylation over each genomic window. Required for heatMap input
#'
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
#' @examples
#' \dontrun{
#'   cluster1kbwindows <- loadSmoothedWindows(obj, type = "CG", threads = 10, name = "500", smoothTrim = 500)
#' }
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
    returnPctMatrix = TRUE) {

  # get barcodes and paths from amethyst object
  if (is.null(obj@h5paths$path) | is.null(obj@h5paths$barcode)) {
    stop("\nPlease make sure the obj@h5paths slot is constructed correctly. Path and/or barcode columns not found.")
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
    dplyr::filter(n > 1) |>
    dplyr::pull(membership)

  groups <- as.list(unique(membership$membership[membership$membership %in% valid_groups]))
  group_results <- list()

  for (gr in groups) {
    members <- rownames(membership |> dplyr::filter(membership == gr))

    if (is.null(obj@h5paths$prefix)) {
      member_paths <- h5paths[match(members, h5paths$barcode), "path"]
    } else {
      h5paths$cell_ids <- paste0(obj@h5paths$prefix, obj@h5paths$barcode)
      member_paths <- h5paths[match(members, h5paths$cell_ids), "path"]
      members <- h5paths[match(members, h5paths$cell_ids), "barcode"]
    }

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
        chunk_results <- chunk_results[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE), n = sum(c_nz + t_nz, na.rm = TRUE)), by = .(chr, start, end)]
        chunk_results <- chunk_results[n >= nmin, .(chr, start, end, c, t)]
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


