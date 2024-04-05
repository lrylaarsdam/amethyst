############################################################################################################################
#' @title geneGeneM
#' @description Helper function to retrieve the relevant portion of the hdf5 file for a given gene
#' @param obj The object containing the path to the hdf5 file
#' @param gene Which gene to fetch hdf5 indexes
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
#' @param cells Optional parameter to retrieve indexes for a subset of cells
#' @return Returns a data frame listing the start and count for each barcode with reads covering the requested gene
#' @export
#' @examples getGeneM(obj = obj, gene = "RBFOX3", type = "CH", cells = c("ATTGAGGATACCAATTCCATCACGAGCG", "ATTGAGGATAACGCTTCGTCCGCCGATC"))
getGeneM <- function(obj,
                     gene,
                     type,
                     cells = NULL) {
  # Check if the gene has been indexed. Return error if not.
  if (!(gene %in% names(obj@index[[type]]))) {
    print(paste0("Error: ", gene, " has not been indexed."))
  }

  # Make data frame of all indexes for the gene
  indexes <- as.data.frame(obj@index[[type]][[gene]])
  if (type == "promoters") {
    type <- "CG"
  }

  # Retrieve barcodes of cells to pull from hdf5. Only pull barcodes containing values covering the gene of interest.
  if (!is.null(cells)) {
    ids <- intersect(indexes$cell_id, cells)
  } else {
    ids <- indexes$cell_id
  }

  # Read in the portion of the hdf5 file with reads covering the given gene for all requested cells
  gene_m <- purrr::map(
    .x = ids,
    .f = function(x) {
      dataset_name <- paste0(type, "/", x)
      start <- indexes |> dplyr::filter(cell_id == x) |> dplyr::pull(start)
      count <- indexes |> dplyr::filter(cell_id == x) |> dplyr::pull(count)
      rhdf5::h5read(obj@h5path, name = dataset_name, start = start, count = count, stride = 1)
    })

  # The result will be a list. Name the list with the corresponding cell barcodes.
  names(gene_m) <- ids

  # Row bind the list into one big list
  gene_m <- do.call(rbind, gene_m)

  # Make a new column containing the barcodes and remove unnecessary characters added automatically by R
  gene_m <- gene_m |>
    tibble::rownames_to_column(var = "cell_id") |>
    dplyr::mutate(cell_id = sub("\\..*$", "", cell_id))
  output <- gene_m
}


############################################################################################################################
#' @title makeWindows
#' @description Calculate methylation levels across a fixed genomic window or gene body.
#'
#' @param obj Object for which to calculate methylation windows
#' @param groupBy Optional metadata parameter to group cells by. If included, the function will calculate mean methylation of a given window over the group
#' @param genes Optional genes to include. Do not specify both genes and stepsize
#' @param stepsize If specified, the methylation levels will be calculated over genomic windows of the indicated size
#' @param metric Calculate either percent, score, or ratio
#' @param threads Enables multithreading
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
#'
#' @return Returns a data frame of the calculated methylation windows with the genomic region or gene name as rows and cell barcode or group as columns
#' @export
#' @examples makeWindows(obj = obj, groupBy = "cluster_id", genes = c("SATB2", "GAD1"), type = "CH", metric = "percent", threads = 10)
#' @examples makeWindows(obj = obj, stepsize = 100000, type = "CG", metric = "score")
makeWindows <- function(
    obj,
    type,
    metric = "percent",
    groupBy = NULL,
    genes = NULL,
    stepsize = NULL,
    threads = 1) {

  # check appropriate metric was specified
  if (!metric %in% c("percent", "score", "ratio")) {
    stop("Metric calculated must be one of either percent, score, or ratio.")
  }

  # set up multithreading
  if (threads > 1) {
    future::plan(future::multicore, workers = threads)
  }

  # get all cell IDs in h5 file
  if (is.null(obj@metadata)) {
    barcodes <- rhdf5::h5ls(obj@h5path)
    barcodes <- as.list(unique(barcodes |> dplyr::filter(otype == "H5I_DATASET") |> dplyr::pull(name)))
  } else {
    barcodes <- as.list(rownames(obj@metadata))
  }

  # If stepsize is specified, make fixed size genomic windows
  if (!is.null(stepsize)) {
    windows <- future_map(barcodes, function(bar) {
      h5 <- data.table::data.table(rhdf5::h5read(obj@h5path, name = paste0({{type}}, "/", bar)))
      h5 <- h5[, window := paste0(chr, "_", plyr::round_any(pos, stepsize, floor), "_", plyr::round_any(pos, stepsize, ceiling))]
      meth_cell <- h5[, round(sum(c != 0) / (sum(c != 0) + sum(t != 0)), 4)] # determine global methylation level
      meth_window <- h5[, .(value = round(sum(c != 0) / (sum(c != 0) + sum(t != 0)), 3)), by = window]

      if (metric == "percent") {
        summary <- meth_window[, value := (value * 100)]
      }
      if (metric == "score") {
        summary <- meth_window[, value := round((ifelse(value - meth_cell > 0,
                                                        (value - meth_cell)/(1 - meth_cell),
                                                        (value - meth_cell)/meth_cell)), 3)]
      }
      if (metric == "ratio") {
        summary <- meth_window[, value := round(value / meth_cell, 3)]
      }
      setnames(summary, "value", bar)
    })
    windows_merged <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE), windows)

    # order window rows by genomic location
    windows_merged <- windows_merged[!grepl("([^_]*_){3,}", window)] # removes alternative loci
    windows_merged[, c("chr", "start", "end") := {
      split_parts <- tstrsplit(window, "_", fixed = TRUE)
      list(split_parts[[1]], as.numeric(split_parts[[2]]), as.numeric(split_parts[[3]]))
    }]
    setorder(windows_merged, chr, start)
    windows_merged[, c("chr", "start", "end") := NULL]
    windows_merged <- windows_merged |> tibble::column_to_rownames(var = "window")
  }

  if (!is.null(genes)) {
    windows <- future_map(as.list(genes), function(x) {
      genem <- data.table::data.table(getGeneM(obj = {{obj}}, gene = x, type = {{type}}))
      if (metric == "percent") {
        summary <- genem[, .(value = round(sum(c != 0) * 100 / (sum(c != 0) + sum(t != 0)), 2)), by = cell_id]
      }
      if (metric == "score" | metric == "ratio") {
        if (type == "CG" | type == "promoters") {
          glob_m <- obj@metadata |> dplyr::mutate(pctm = mcg_pct/100) |> dplyr::select(pctm) |> tibble::rownames_to_column(var = "cell_id")
        } else if (type == "CH") {
          glob_m <- obj@metadata |> dplyr::mutate(pctm = mch_pct/100) |> dplyr::select(pctm) |> tibble::rownames_to_column(var = "cell_id")
        }
        summary <- genem[, .(value = round(sum(c != 0) / (sum(c != 0) + sum(t != 0)), 2)), by = cell_id]
        summary <- merge(summary, glob_m, by = "cell_id")
        if (metric == "score") {
          summary <- summary[, value := round((ifelse(value - pctm > 0,
                                                      (value - pctm)/(1 - pctm),
                                                      (value - pctm)/pctm)), 3), by = cell_id][, pctm := NULL]
        }
        if (metric == "ratio") {
          summary <- summary[, value := round(value / pctm, 3), by = cell_id][, pctm := NULL]
        }
      }
      setnames(summary, "value", x)
    })
    windows_merged <- Reduce(function(x, y) merge(x, y, by = "cell_id", all = TRUE), windows) |> tibble::column_to_rownames(var = "cell_id")
    windows_merged <- as.data.frame(t(windows_merged))
  }

  if (!is.null(groupBy)) {
    membership <- obj@metadata |> dplyr::select(groupBy)
    colnames(membership) <- "membership"
    groups <- as.list(unique(membership$membership))
    aggregated <- list()
    for (i in groups) {
      members <- rownames(membership |> dplyr::filter(membership == i))
      aggregated[[i]] <- as.data.frame(rowMeans(windows_merged[, members], na.rm = TRUE))
      colnames(aggregated[[i]]) <- paste0(i)
      aggregated[[i]] <- aggregated[[i]] |> tibble::rownames_to_column(var = "window")
    }
    windows_merged <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE), aggregated) |> tibble::column_to_rownames(var = "window")
  }

  if (threads > 1) {
    future::plan(NULL)
    gc()
  }
  windows_merged

}


############################################################################################################################
#' @title addMetadata
#' @description Add new cell metadata to the metadata slot of a pre-existing object. Make sure row names are cell IDs
#'
#' @param obj Object to add metadata to
#' @param metadata Data frame, with row names as cell IDs, containing the new metadata to be added
#' @return Updates obj to contain the new metadata concatenated to the old metadata
#' @export
#' @examples obj <- addMetadata(obj = obj, metadata = annotations)
addMetadata <- function(obj,
                        metadata) {
  if (length(intersect(rownames(obj@metadata), rownames(metadata))) == 0) {
    print("No intersection detected; check metadata structure. Make sure row names are cell IDs.")
    output <- obj
  } else {
    obj@metadata <- dplyr::left_join(
      obj@metadata |> tibble::rownames_to_column(var = "cell_id"),
      metadata |> tibble::rownames_to_column(var = "cell_id"),
      by = "cell_id"
    ) |>
      tibble::column_to_rownames(var = "cell_id")
    output <- obj
  }
}

############################################################################################################################
#' @title addCellInfo
#' @description Add the information contained in the cellInfo file outputs to the obj@metadata slot
#'
#' @param obj Amethyst object to add cell info metadata
#' @param file Path of the cellInfo.txt file
#' @return Returns the same amethyst object with information contained in the cellInfo file added to the obj@metadata slot
#' @export
#' @examples obj <- addCellInfo(obj = obj, file = "~/Downloads/cellInfo.txt")
addCellInfo <- function(obj,
                        file) {

  cellInfo <- read.table(file, sep = "\t", header = F)
  if (ncol(cellInfo) == 6) {
    colnames(cellInfo) <- c("cell_id", "cov", "cg_cov", "mcg_pct", "ch_cov", "mch_pct")
  } else if (ncol(cellInfo) == 9) {
    colnames(cellInfo) <- c("cell_id", "cov", "cg_cov", "mcg_pct", "ch_cov", "mch_pct", "n_frag", "n_pair", "n_single")
  } else {
    stop("The cell file is not in the expected format. Please check your input.")
  }
  obj <- addMetadata(obj, metadata = cellInfo |> tibble::column_to_rownames(var = "cell_id"))
  output <- obj
}

############################################################################################################################
#' @title impute
#' @description Impute values with the Rmagic package https://github.com/KrishnaswamyLab/MAGIC
#'
#' @param obj Amethyst object to perform imputation on
#' @param matrix Name of the matrix contained in the genomeMatrices slot to perform imputation on
#' @param npca Number of principle components to use
#' @return Adds a new data frame to the genomeMatrices slot with imputed values
#' @export
#' @examples obj <- impute(obj, matrix = "gene_ch")
impute <- function(obj,
                   matrix,
                   npca = 10) {
  name <- matrix # store name of matrix for output
  matrix <- obj@genomeMatrices[[matrix]]
  num_na <- sum(is.na(matrix))
  pct_na <- round(sum(is.na(matrix)) * 100 / (ncol(matrix) * nrow(matrix)), 3)
  if (num_na > 0) {
    print(paste0("Warning: Replacing ", num_na, " NAs with zeros (", pct_na, "% of matrix) to perform imputation."))
  }
  matrix[is.na(matrix)] <- 0
  matrix_imputed <- Rmagic::magic(t(matrix), npca = npca)[["result"]]
  obj@genomeMatrices[[paste0(name, "_imputed")]] <- as.data.frame(t(matrix_imputed))
  output <- obj
}

############################################################################################################################
#' @title aggregateMatrix
#' @description Aggregate genomic window matrices by a categorical variable contained in the metadata
#'
#' @param obj Amethyst object containing the matrix to be aggregated
#' @param matrix Name of matrix contained in the genomeMatrices slot to aggregate
#' @param groupBy Parameter in the metadata to aggregate over
#' @return Returns a matrix containing the mean value per group
#' @export
#' @examples obj <- aggregateMatrix(obj, matrix = "ch_100k_pct", groupBy = cluster_id, name = "cluster_ch_100k_pct")
aggregateMatrix <- function(obj,
                            matrix,
                            groupBy) {

  matrix <- obj@genomeMatrices[[matrix]]
  membership <- obj@metadata |> dplyr::select(groupBy)
  colnames(membership) <- "membership"
  groups <- as.list(unique(membership$membership))
  aggregated <- list()

  for (i in groups) {
    members <- rownames(membership |> dplyr::filter(membership == i))
    aggregated[[i]] <- as.data.frame(rowMeans(matrix[, members], na.rm = TRUE))
    colnames(aggregated[[i]]) <- paste0(i)
    aggregated[[i]] <- aggregated[[i]] |> tibble::rownames_to_column(var = "window")
  }
  aggregated <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE), aggregated)

  # reorder if genomic windows
  if (sum(grepl("chr", aggregated$window)) > 0) {
    # order
    setDT(aggregated)
    aggregated <- aggregated[!grepl("([^_]*_){3,}", window)] # removes alternative loci
    aggregated[, c("chr", "start", "end") := {
      split_parts <- tstrsplit(window, "_", fixed = TRUE)
      list(split_parts[[1]], as.numeric(split_parts[[2]]), as.numeric(split_parts[[3]]))
    }]
    setorder(aggregated, chr, start)
    aggregated[, c("chr", "start", "end") := NULL]
  }
  aggregated <- aggregated |> tibble::column_to_rownames(var = "window")
}

############################################################################################################################
#' @title makeSlidingWindows
#' @description Calculate mean percent methylation over 1500 base pair genomic windows shifted every 500 bp and aggregate by cluster.
#' Intended to generate input for the tileGeneM function.
#'
#' @param obj Amethyst object containing the path to the h5 file over which to calculate 500 bp windows
#' @param type Type of methylation to calculate; e.g. "CG" or "CH"
#' @param threads Enables multithreading
#' @return Returns a data frame with % methylation over 1500 base pair genomic windows shifted every 500 bp aggregated by cluster
#' @export
#' @examples cluster_cg_500bp_windows <- makeSlidingWindows(obj = obj, type = "CG", threads = 20)
makeSlidingWindows <- function(
    obj,
    type,
    threads = 1) {

  if (threads > 1) {
    plan(multicore, workers = threads)
  }

  # get all cell IDs in h5 file
  if (is.null(obj@metadata)) {
    barcodes <- h5ls(obj@h5path)
    barcodes <- as.list(unique(barcodes %>% dplyr::filter(otype == "H5I_DATASET") %>% dplyr::pull(name)))
  } else {
    barcodes <- as.list(rownames(obj@metadata))
  }

  # get all possible 500 bp positions so you don't accidentally average across missing values
  bed <- read.table("/home/groups/ravnica/projects/amethyst/bedfiles/hg38_500bp_intervals.bed", sep = "\t")
  colnames(bed) <- c("chr", "start", "end")
  bed <- bed %>% arrange(chr, start, end) %>% dplyr::mutate(window = paste0(chr, "_", start, "_", end)) %>% dplyr::select(window)

  # calculate % methylation over 500 bp windows for all cells
  windows <- future_map(.x = barcodes, .f = function(x) {
    barcode_name <- sub("\\..*$", "", x)
    data <- h5read(obj@h5path, name = paste0(type, "/", x))
    setDT(data)
    data[, window := paste0(chr, "_", plyr::round_any(pos, 500, floor), "_", plyr::round_any(pos, 500, ceiling))]
    result <- data[, .(value = round(sum(c != 0) * 100 / (sum(c != 0) + sum(t != 0)), 1)), by = window]
    setnames(result, "value", barcode_name)
    result <- left_join(bed, result, by = c("window"))
  }, .progress = TRUE)

  if (all(sapply(windows[-1], function(df) identical(windows[[1]]$window, df$window)))) { # check to make sure all the gene columns are in the same order before cbinding
    result <- do.call(cbind, c(windows[[1]][, 1:2], lapply(windows[-1], function(df) df[, 2]))) # cbind second column of all dataframes with first
  } else {
    result <- purrr::reduce(windows, full_join, by = c("window"))
  }

  # Calculate the rolling mean for each column
  result <- separate(result, col = "window", into = c("chr", "start", "end"), sep = "_")
  smoothed <- result
  setDT(smoothed)

  smooth_start <- smoothed[, lapply(.SD, function(x) frollapply(x, n = 3, FUN = function(y) y[1], align = "center", fill = NA)), .SDcols ="start"]
  smooth_end <- smoothed[, lapply(.SD, function(x) frollapply(x, n = 3, FUN = function(y) y[3], align = "center", fill = NA)), .SDcols ="end"]
  averages <- smoothed[, lapply(.SD, function(x) round(frollmean(x, n = 3, align = "center", fill = NA, na.rm = TRUE), 2)), .SDcols = -c("window", "chr", "start", "end")]

  # Add the averages as new columns to your original data.table
  smoothed <- smoothed[, .SD, .SDcols = 1:4]
  smoothed[, paste0("smooth_start") := smooth_start]
  smoothed[, paste0("smooth_end") := smooth_end]
  smoothed$smooth_window <- paste0(smoothed$chr, "_", smoothed$smooth_start, "_", smoothed$smooth_end)
  smoothed[, flag := ifelse(shift(chr, 1) != shift(chr, -2), "FLAG", "OK")] # remove rows over the chromosome junctions
  smoothed[, names(averages) := averages]
  smoothed <- smoothed[flag == "OK"]
  smoothed <- smoothed[, !("flag") , with = FALSE]

  # If aggregating by group
  groups <- as.list(unique(obj@metadata[["cluster_id"]]))
  aggregated <- list()
  for (i in groups) {
    aggregated[[i]] <- rowMeans(as.data.frame(smoothed)[rownames(obj@metadata[obj@metadata$cluster_id == i , ])], na.rm = TRUE)
  }
  aggregated <- do.call(data.frame, aggregated)
  colnames(aggregated) <- groups
  aggregated <- cbind(as.data.frame(smoothed)[1:7], aggregated)
  output <- aggregated

}

############################################################################################################################
#' @title makeFuzzyGeneWindows
#' @description Rapidly calculate % methylation over all genes by averaging 500 bp window values.
#' Slightly less accurate than the alternative indexing approach for short genes.
#'
#' @param obj Amethyst object for which to calculate % methylation over all genes
#' @param type Type of methylation to calculate; e.g. "CH".
#' "CG" is an option but may be less biologically relevant
#' @param threads Enables multithreading
#' @return Returns a matrix just like the output of makeWindows with rows as genes, columns as cells, and values as % methylation
#' @export
#'
#' @examples gene_ch <- makeFuzzyGeneWindows(obj = obj, type = "CH", threads = 20)
makeFuzzyGeneWindows <- function(
    obj,
    type,
    threads = 1) {

  if (threads > 1) {
    plan(multicore, workers = threads)
  }

  # download bed file
  bed <- read.table("/home/groups/ravnica/projects/amethyst/bedfiles/hg38_mainChr_500bp_geneannot.bed", header = T, sep = '\t')
  setDT(bed)
  bed <- bed[, .(window, gene)]

  # get all cell IDs in h5 file
  if (is.null(obj@metadata)) {
    barcodes <- rhdf5::h5ls(obj@h5path)
    barcodes <- as.list(unique(barcodes %>% dplyr::filter(otype == "H5I_DATASET") %>% dplyr::pull(name)))
  } else {
    barcodes <- as.list(rownames(obj@metadata))
  }

  # calculate % methylation over 500 bp windows for all cells
  windows <- future_map(.x = barcodes, .f = function(x) {
    barcode_name <- sub("\\..*$", "", x)
    data <- data.table::data.table(rhdf5::h5read(obj@h5path, name = paste0(type, "/", x)))
    data[, window := paste0(chr, "_", plyr::round_any(pos, 500, floor), "_", plyr::round_any(pos, 500, ceiling))]
    result <- data[, .(value = round(sum(c != 0) * 100 / (sum(c != 0) + sum(t != 0)), 1)), by = window]
    setnames(result, "value", barcode_name)
    result <- left_join(bed, result, by = c("window"))
    result <- result[, lapply(.SD, function(x) round(mean(x, na.rm = TRUE), 2)), by = gene, .SDcols = barcode_name]
    setorder(result, gene)
  }, .progress = TRUE)

  if (all(sapply(windows[-1], function(df) identical(windows[[1]]$gene, df$gene)))) { # check to make sure all the gene columns are in the same order before cbinding
    result <- do.call(cbind, c(windows[[1]][, 1:2], lapply(windows[-1], function(df) df[, 2]))) # cbind second column of all dataframes with first
  } else {
    result <- purrr::reduce(windows, full_join, by = c("gene"))
  }
  result <- result |> column_to_rownames(var = "gene")
}


