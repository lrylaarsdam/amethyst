############################################################################################################################
#' @title makeWindows
#' @description Calculate methylation levels across a fixed genomic window, bed file coordinates, or gene body.
#'
#' @param obj Object for which to calculate methylation windows
#' @param type What type of methylation to retrieve; i.e. "CH" or "CG"
#' @param genes Optional genes to calculate % methylation over body
#' @param stepsize If specified, the methylation levels will be calculated over genomic windows of the indicated size
#' @param bed Optional bed file to input. If specified, methylation levels will be calculated over input regions
#' @param metric Calculate either methylation percent, score, or ratio
#' @param groupBy Optional metadata parameter to group cells by. If included, the function will calculate mean methylation of a given window over the group
#' @param threads Enables multithreading
#' @param index If calculating genomic windows, specify the name of the chr index in the index slot.
#' This index should contain the coordinates in the hdf5 file corresponding to each chromosome. Reduces memory constraints.
#' @param futureType Method of parallelization, i.e. "multicore" (recommended) or "multisession". Multisession is more
#' memory-efficient, but make sure your workspace is as small as possible before using.
#' @param nmin Minimum number of observations for the window to be included.
#' @param save Boolean indicating whether to save the intermediate output by chromosome. Good for very large datasets.
#' @return Returns a data frame with columns as cells and rows as either genes or genomic windows
#' @export
#' @examples makeWindows(obj, type = "CH", genes = c("SATB2", "TBR1", "PACS1"), metric = "percent")
#' @importFrom data.table data.table tstrsplit := setorder setDT rbindlist
#' @importFrom dplyr pull filter select mutate
#' @importFrom furrr future_map
#' @importFrom future plan multicore sequential
#' @importFrom plyr round_any
#' @importFrom rhdf5 h5ls h5read
#' @importFrom tibble column_to_rownames rownames_to_column
makeWindows <- function(
    obj,
    type,
    genes = NULL,
    promoter = FALSE,
    stepsize = NULL,
    bed = NULL,
    metric = "percent",
    index = paste0("chr_", tolower(type)),
    groupBy = NULL,
    threads = 1,
    futureType = "multicore",
    nmin = 2,
    save = FALSE) {

  # check appropriate metric was specified
  if (!metric %in% c("percent", "score", "ratio")) {
    stop("Metric calculated must be one of either percent, score, or ratio.")
  }

  # check parameters
  if (sum(!is.null(bed), !is.null(stepsize), !is.null(genes)) > 1) {
    stop("Please only specify a fixed step size, gene list, or input bed file.")
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

  # check paths exist
  if (is.null(obj@h5paths)) {
    stop("Please generate the path list for each barcode and store in the obj@h5paths slot.")
  }

  # check index
  if (is.null(obj@index[[index]])) {
    stop("Please index which rows in the h5 files correspond to each chromosome using indexChr.")
  }


  #define chr groups
  chr_groups <- as.list(names(obj@index[[index]]))

  by_chr <- list()

  if (!is.null(bed)) {

    if (sum(class(bed) %in% "character") > 0) {
      bed <- data.table::data.table(read.table(bed))
    } else if (sum(class(bed) %in% c("data.frame", "data.table")) > 0) {
      data.table::setDT(bed)
    }
    if (ncol(bed) > 3) {
      stop("\nPlease make sure the bed file consists of three columns - chr, start, end - with no header.\n")
    }
    setnames(bed, c("chr", "start", "end"))
    chr_groups <- as.list(unique(as.character(bed$chr)))

    bed <- split(bed, bed$chr)

    for (chr in chr_groups) {

      # get index positions
      sites <- obj@index[[index]][[chr]] # get chr index for h5 file

      # get paths
      barcodes <- as.list(rownames(obj@h5paths[rownames(obj@h5paths) %in% sites$cell_id, , drop = FALSE]))
      paths <- as.list(obj@h5paths[rownames(obj@h5paths) %in% sites$cell_id, , drop = FALSE]$paths)

      # add up sum c and sum t in member cells
      windows <- furrr::future_pmap(.l = list(paths, barcodes), .f = function(path, barcode) {
        tryCatch({
          bed_tmp <- copy(bed[[chr]])
          barcode_name <- sub("\\..*$", "", barcode)

          meth_cell <- obj@metadata[barcode, paste0("m", tolower(type), "_pct")]/100 # pull global methylation level from metadata
          h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode),
                                                     start = sites$start[sites$cell_id == barcode],
                                                     count = sites$count[sites$cell_id == barcode])) # read in 1 chr at a time

          meth_window <- h5[bed_tmp, on = .(chr = chr, pos >= start, pos <= end), nomatch = 0L, .(chr, start, end, pos = x.pos, c, t)]
          meth_window <- meth_window[, .(value = round((sum(c != 0)/ (sum(c != 0) + sum(t != 0))), 4),
                                         n = sum(c + t, na.rm = TRUE)), by = .(chr, start, end)]
          meth_window <- meth_window[n >= nmin, ][, n := NULL]

          if (metric == "percent") { summary <- meth_window[, value := (value * 100)] }
          if (metric == "score") { summary <- meth_window[, value := round((ifelse((value - meth_cell) > 0,
                                                                                   (value - meth_cell)/(1 - meth_cell),
                                                                                   (value - meth_cell)/meth_cell)), 3)] }
          if (metric == "ratio") { summary <- meth_window[, value := round(value / meth_cell, 3)] }
          data.table::setnames(summary, "value", barcode_name)
        }, error = function(e) {
          cat("Error processing data for barcode", barcode, ":", conditionMessage(e), "\n")
          return(NULL)  # Return NULL or any other value indicating failure
        })
      }, .progress = TRUE)

      # Reduce merge the data.tables per cell in chunks (way faster for large datasets)
      windows <- split(windows, ceiling(seq_along(windows)/1000))
      windows_merged <- Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = FALSE),
                               furrr::future_map(windows, ~ Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = FALSE), .x), .progress = TRUE))
      by_chr[[chr]] <- windows_merged
      cat("\nCompleted ", chr,"\n")
    }

    windows_merged <- data.table::rbindlist(by_chr, fill = TRUE, use.names = TRUE)
    data.table::setorder(windows_merged, chr, start)
    windows_merged <- windows_merged[, window := paste0(chr, "_", start, "_", end)]
    windows_merged[, c("chr", "start", "end") := NULL]
    windows_merged <- windows_merged |> tibble::column_to_rownames(var = "window")

  }

  # If stepsize is specified, make fixed size genomic windows
  if (!is.null(stepsize)) {

    for (chr in chr_groups) {

      # get index positions
      sites <- obj@index[[index]][[chr]] # get chr index for h5 file

      # get paths
      barcodes <- as.list(rownames(obj@h5paths[rownames(obj@h5paths) %in% sites$cell_id, , drop = FALSE]))
      paths <- as.list(obj@h5paths[rownames(obj@h5paths) %in% sites$cell_id, , drop = FALSE]$paths)

      windows <- furrr::future_pmap(.l = list(paths, barcodes), .f = function(path, barcode) {
        if (futureType == "multisession") {
          options(scipen = 999)
        }
        tryCatch({
          h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode),
                                                     start = sites$start[sites$cell_id == barcode],
                                                     count = sites$count[sites$cell_id == barcode])) # read in 1 chr at a time
          h5 <- h5[pos %% stepsize == 0, pos := pos + 1] # otherwise sites exactly divisible by stepsize will be their own window
          h5 <- h5[, window := paste0(chr, "_", plyr::round_any(pos, stepsize, floor), "_", plyr::round_any(pos, stepsize, ceiling))]
          meth_cell <- obj@metadata[barcode, paste0("m", tolower(type), "_pct")]/100 # pull global methylation level from metadata
          meth_window <- h5[, .(value = round(sum(c != 0) / (sum(c != 0) + sum(t != 0)), 3), n = sum(c + t, na.rm = TRUE)), by = window]
          meth_window <- meth_window[n >= nmin, .(window, value)]

          if (metric == "percent") { summary <- meth_window[, value := (value * 100)] }
          if (metric == "score") { summary <- meth_window[, value := round((ifelse((value - meth_cell > 0),
                                                                                   (value - meth_cell)/(1 - meth_cell),
                                                                                   (value - meth_cell)/meth_cell)), 3)] }
          if (metric == "ratio") { summary <- meth_window[, value := round(value / meth_cell, 3)] }
          data.table::setnames(summary, "value", barcode)
        }, error = function(e) {
          # Handle the error here, for example:
          cat("Error processing data for barcode", barcode, ":", conditionMessage(e), "\n")
          return(NULL)  # Return NULL or any other value indicating failure
        })
      }, .progress = TRUE)

      # Reduce merge the data.tables per cell in chunks (way faster for large datasets)
      windows <- split(windows, ceiling(seq_along(windows)/1000))
      windows_merged <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE),
                               furrr::future_map(windows, ~ Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE), .x), .progress = TRUE))
      if (save) {
        saveRDS(windows_merged, paste0("tmp_", type, "_", metric, "_", chr, "_", stepsize, "kb_windows_nmin", nmin, ".RData"))
      }
      by_chr[[chr]] <- windows_merged
      cat("\nCompleted ", chr,"\n")
    }

    windows_merged <- data.table::rbindlist(by_chr, fill = TRUE, use.names = TRUE)

    # sort by chr and start just in case
    windows_merged <- windows_merged[, c("chr", "start", "end") := {
      split_parts <- data.table::tstrsplit(window, "_", fixed = TRUE)
      list(split_parts[[1]], as.numeric(split_parts[[2]]), as.numeric(split_parts[[3]]))
    }]
    windows_merged <- windows_merged[(end - start) == stepsize]
    data.table::setorder(windows_merged, chr, start)
    windows_merged[, c("chr", "start", "end") := NULL]

    # filter alternative loci and sex chromosomes
    windows_merged <- windows_merged |> tibble::column_to_rownames(var = "window")
    windows_merged <- windows_merged[!sapply(rownames(windows_merged), function(name) length(strsplit(name, "_")[[1]]) > 3 || grepl("chrEBV|chrM|KI", name)), ]

  }

  if (!is.null(genes)) {
    if (is.null(obj@ref)) {
      stop("Please add a genome annotation file to the obj@ref slot using the makeRef function. \nMake sure the build is the same as what was used for alignment.")
    }

    if (!promoter) {
      bed <- obj@ref |> dplyr::filter(gene_name %in% genes & type == "gene") |>
        dplyr::select(seqid, start, end, gene_name) |>
        dplyr::distinct(gene_name, .keep_all = TRUE) |>
        dplyr::rename("chr" = "seqid", "gene" = "gene_name") |>
        data.table::data.table()
    } else if (promoter) {
      bed <- obj@ref |> dplyr::filter(gene_name %in% genes & type == "gene") |>
        dplyr::select(seqid, start, end, gene_name, strand) |>
        dplyr::distinct(gene_name, .keep_all = TRUE) |>
        dplyr::rename("chr" = "seqid", "gene" = "gene_name") |>
        data.table::data.table()
      bed[, `:=`(p_start = ifelse(strand == "+", start - 1500, end - 1500),
                 p_end = ifelse(strand == "+", start + 1500, end + 1500))] # make bed file for promoter cg
      bed[, `:=`(start = NULL, end = NULL, strand = NULL)]
      setnames(bed, c("chr", "gene", "start", "end"))
    }

    chr_groups <- as.list(unique(as.character(bed$chr)))
    bed <- split(bed, bed$chr)

    for (chr in chr_groups) {

      # get index positions
      sites <- obj@index[[index]][[chr]] # get chr index for h5 file

      # get paths
      barcodes <- as.list(rownames(obj@h5paths[rownames(obj@h5paths) %in% sites$cell_id, , drop = FALSE]))
      paths <- as.list(obj@h5paths[rownames(obj@h5paths) %in% sites$cell_id, , drop = FALSE]$paths)

      # add up sum c and sum t in member cells
      windows <- furrr::future_pmap(.l = list(paths, barcodes), .f = function(path, barcode) {
        tryCatch({
          bed_tmp <- copy(bed[[chr]])
          barcode_name <- sub("\\..*$", "", barcode)
          data <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode),
                                                       start = sites$start[sites$cell_id == barcode],
                                                       count = sites$count[sites$cell_id == barcode])) # read in 1 chr at a time
          result <- data[bed_tmp, on = .(chr = chr, pos >= start, pos <= end), nomatch = 0L, .(chr, start, end, pos = x.pos, c, t, gene)]
          summary <- result[, .(cell_id =
                                  round(((sum(c != 0) * 100 )/ (sum(c != 0) + sum(t != 0))), 2),
                                n = sum(c + t, na.rm = TRUE)),
                            by = .(gene)]
          summary <- summary[n >= nmin, .(gene, cell_id)]
          setnames(summary, "cell_id", barcode_name)
        }, error = function(e) {
          cat("Error processing data for barcode", barcode, ":", conditionMessage(e), "\n")
          return(NULL)  # Return NULL or any other value indicating failure
        })
      }, .progress = TRUE)

      # Reduce merge the data.tables per cell in chunks (way faster for large datasets)
      windows <- split(windows, ceiling(seq_along(windows)/1000))
      windows_merged <- Reduce(function(x, y) merge(x, y, by = c("gene"), all = TRUE, sort = FALSE),
                               furrr::future_map(windows, ~ Reduce(function(x, y) merge(x, y, by = c("gene"), all = TRUE, sort = FALSE), .x), .progress = TRUE))
      by_chr[[chr]] <- windows_merged
      cat(paste0("\nCompleted ", bed[[chr]][["gene"]], " in ", chr))
    }

    windows_merged <- data.table::rbindlist(by_chr, fill = TRUE, use.names = TRUE)
    windows_merged <- windows_merged |> tibble::column_to_rownames(var = "gene")

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
    windows_merged <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE), aggregated) |> tibble::column_to_rownames(var = "window")
  }

  if (threads > 1) {
    future::plan(future::sequential)
    gc()
  }

  return(windows_merged)

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
#' @importFrom dplyr left_join
#' @importFrom tibble column_to_rownames rownames_to_column
addMetadata <- function(obj,
                        metadata) {
  if (is.null(obj@metadata)) {
    obj@metadata <- metadata
  } else {
    if (length(intersect(rownames(obj@metadata), rownames(metadata))) == 0) {
      stop("No intersection detected; check metadata structure. Make sure row names are cell IDs.")
    } else {
      obj@metadata <- dplyr::left_join(
        obj@metadata |> tibble::rownames_to_column(var = "cell_id"),
        metadata |> tibble::rownames_to_column(var = "cell_id"),
        by = "cell_id"
      ) |>
        tibble::column_to_rownames(var = "cell_id")
    }
  }
  return(obj)
}

############################################################################################################################
#' @title addAnnot
#' @description Add information containined in the .annot file to the Amethyst object
#'
#' @param file Path to .annot file
#' @param obj Amethyst object to add metadata to
#' @param name Column name to add in the metadata
#'
#' @return Updates obj to contain the new metadata concatenated to the old metadata
#' @export
#' @examples obj <- addMetadata(obj = obj, metadata = annotations)
#' @importFrom dplyr left_join
#' @importFrom tibble column_to_rownames rownames_to_column
addAnnot <- function(obj,
                     file,
                     name) {
  annot <- utils::read.table(file, sep = "\t", header = F)
  if (ncol(annot) == 2) {
    colnames(annot) <- c("cell_id", name)
  } else if (ncol(annot) > 2) {
    cat("\nAdditional annot columns will need to be renamed.\n")
    colnames(annot) <- c("cell_id", name)
  }
  if (is.null(obj@metadata)) {
    obj@metadata <- annot |> tibble::column_to_rownames(var = "cell_id")
  } else {
    obj <- addMetadata(obj, metadata = annot |> tibble::column_to_rownames(var = "cell_id"))
  }
  output <- obj
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
#' @importFrom tibble column_to_rownames
#' @importFrom utils read.table
addCellInfo <- function(obj,
                        file) {
  cellInfo <- utils::read.table(file, sep = "\t", header = F)
  if (ncol(cellInfo) == 6) {
    colnames(cellInfo) <- c("cell_id", "cov", "cg_cov", "mcg_pct", "ch_cov", "mch_pct")
  } else if (ncol(cellInfo) == 10) {
    colnames(cellInfo) <- c("cell_id", "cov", "cg_cov", "mcg_pct", "ch_cov", "mch_pct", "n_frag", "n_pair", "n_single", "xpct_m")
  } else if (ncol(cellInfo) != 6 && ncol(cellInfo) != 10) {
    cat("\nThe cell file is not in the expected format of 6 or 10 columns. Columns after 10 will not be named.\n")
    names <- c("cell_id", "cov", "cg_cov", "mcg_pct", "ch_cov", "mch_pct", "n_frag", "n_pair", "n_single", "xpct_m")
    colnames(cellInfo) <- names[1:ncol(cellInfo)]
  }
  if (is.null(obj@metadata)) {
    obj@metadata <- cellInfo |> tibble::column_to_rownames(var = "cell_id")
  } else {
    obj <- addMetadata(obj, metadata = cellInfo |> tibble::column_to_rownames(var = "cell_id"))
  }
  output <- obj
}

############################################################################################################################
#' @title impute
#' @description Impute values with the Rmagic package https://github.com/KrishnaswamyLab/MAGIC
#'
#' @param obj Amethyst object to perform imputation on
#' @param matrix Name of the matrix contained in the genomeMatrices slot to perform imputation on
#' @param replaceNA Rmagic can't accept NA values. Replace NA values with 0, 1, "mch_pct", or "mcg_pct".
#' @param npca Number of principle components to use
#'
#' @return Returns a data frame of imputed values.
#' @export
#' @examples obj <- impute(obj, matrix = "gene_ch")
#' @importFrom Rmagic magic
impute <- function(obj,
                   matrix,
                   npca = 10,
                   replaceNA = 0) {
  name <- matrix # store name of matrix for output
  matrix <- obj@genomeMatrices[[matrix]]
  num_na <- sum(is.na(matrix))
  pct_na <- round(sum(is.na(matrix)) * 100 / (ncol(matrix) * nrow(matrix)), 3)
  if (num_na > 0) {
    print(paste0("Warning: Replacing ", num_na, " NAs (", pct_na, "% of matrix) to perform imputation."))
  }
  if (replaceNA == 0) {
    matrix[is.na(matrix)] <- 0
  } else if (replaceNA == 1) {
    matrix[is.na(matrix)] <- 1
  } else if (replaceNA %in% c("mch_pct", "mcg_pct")) {
    matrix <- matrix[, colnames(matrix) %in% rownames(obj@metadata)]
    if (replaceNA == "mch_pct") {
      replacement_vector <- obj@metadata$mch_pct[match(colnames(matrix), rownames(obj@metadata))]
    }
    if (replaceNA == "mcg_pct") {
      replacement_vector <- obj@metadata$mcg_pct[match(colnames(matrix), rownames(obj@metadata))]
    }
    # Replace NA values in the matrix with corresponding mch_pct values
    na_indices <- which(is.na(matrix), arr.ind = TRUE)
    matrix[na_indices] <- replacement_vector[na_indices[, 2]]
  } else {
    stop("Default NA replacement must be 0, 1, 'mch_pct', or 'mcg_pct'.")
  }
  matrix_imputed <- Rmagic::magic(t(matrix), npca = npca)[["result"]]
  matrix_imputed <- as.data.frame(t(matrix_imputed))

  return(matrix_imputed)
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
#' @examples cluster_ch_100k_pct <- aggregateMatrix(obj, matrix = "ch_100k_pct", groupBy = "cluster_id")
#' @importFrom data.table setDT tstrsplit
#' @importFrom dplyr select filter
#' @importFrom tibble rownames_to_column column_to_rownames
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
    data.table::setDT(aggregated)
    aggregated <- aggregated[!grepl("([^_]*_){3,}", window)] # removes alternative loci
    aggregated[, c("chr", "start", "end") := {
      split_parts <- data.table::tstrsplit(window, "_", fixed = TRUE)
      list(split_parts[[1]], as.numeric(split_parts[[2]]), as.numeric(split_parts[[3]]))
    }]
    data.table::setorder(aggregated, chr, start)
    aggregated[, c("chr", "start", "end") := NULL]
  }
  aggregated <- aggregated |> tibble::column_to_rownames(var = "window")
}

############################################################################################################################
#' @title subsetObject
#' @description Subset all slots of an amethyst object to only include information from specified barcodes
#'
#' @param obj Amethyst object to subset
#' @param cells Character vector of barcodes to include in the subset
#' @return Returns a new amethyst object with all slots subsetted
#' @export
#' @examples subset <- subsetObject(obj = obj, cells = c("ATTGAGGATAACGCTTCGTCCGCCGATC", "ATTGAGGATAACGCTTCGTCTAAGGTCA"))
subsetObject <- function(obj,
                         cells) {
  subset <- obj
  subset@h5paths <- subset@h5paths[rownames(subset@h5paths) %in% cells, , drop = FALSE]

  # get genomeMatrices to subset
  contains_cells <- sapply(subset@genomeMatrices, function(matrix) {
    any(colnames(matrix) %in% cells)
  })
  which_matrices <- names(subset@genomeMatrices)[contains_cells]
  if (length(which_matrices) != length(subset@genomeMatrices)) {
    cat("\nWarning: any aggregated matrices will still contain information from non-subsetted cells.\n")
  }

  for (i in which_matrices) {
    subset@genomeMatrices[[i]] <- subset@genomeMatrices[[i]][, colnames(subset@genomeMatrices[[i]]) %in% cells]
  }

  # subset reductions
  for (i in names(subset@reductions)) {
    subset@reductions[[i]] <- subset@reductions[[i]][rownames(subset@reductions[[i]]) %in% cells , ]
  }

  # subset indexes
  for (i in names(subset@index)) {
    subset@index[[i]] <- lapply(subset@index[[i]], function(x) {
      x <- x[x$cell_id %in% cells, ]
    })
  }

  # subset metadata
  subset@metadata <- subset@metadata[rownames(subset@metadata) %in% cells, , drop = FALSE]

  return(subset)

}

############################################################################################################################
#' @title combineObject
#' @description Merge a list of amethyst objects into one
#'
#' @param objList List of amethyst objects to merge. Make sure barcode names are unique.
#' @param matrices List of genomeMatrices contained within each amethyst object to merge.
#' @return Returns a new object with merged slots. Dimensionality reductions should be re-run post merging.
#' @export
#' @importFrom tibble column_to_rownames rownames_to_column
#' @examples new <- combineObject(objList = list(obj1, obj2, obj3), genomeMatrices = c("ch_100k_pct", "gene_ch"))
combineObject <- function(objList,
                          genomeMatrices) {
  cat("\nMake sure barcodes are unique when combining objects.")

  combined <- createObject()

  # rbind h5paths
  combined@h5paths <- do.call(rbind, lapply(objList, function(x) {x@h5paths}))

  # merge genomeMatrices
  for (i in genomeMatrices) {
    tomerge <- lapply(objList, function(x) {x@genomeMatrices[[i]] |> tibble::rownames_to_column(var = "window")})
    combined@genomeMatrices[[i]] <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE), tomerge) |>
      tibble::column_to_rownames(var = "window")
  }

  # print reductions warning
  cat("\nDimensionality reductions (i.e., irlba, umap, tsne) should not be merged and will need to be re-run.")

  # merge indexes
  cat("\nOnly indexes and subindexes found in all objects will be merged.")
  common_indexes <- Reduce(intersect, lapply(objList, function(x) names(x@index)))
  for (i in common_indexes) {
    common_subindexes <- Reduce(intersect, lapply(objList, function(x) names(x@index[[i]])))
    for (j in common_subindexes) {
      combined@index[[i]][[j]] <- do.call(rbind, lapply(objList, function(x) {x@index[[i]][[j]]}))
    }
  }

  # merge metadata
  combined@metadata <- do.call(rbind, lapply(objList, function(x) {x@metadata}))

  # pull ref
  if (!is.null(objList[[1]]@ref)) {
    combined@ref <- objList[[1]]@ref
  } else {
    cat("\nThe first object did not contain a reference. The combined ref slot will remain empty.")
  }

  return(combined)

}


############################################################################################################################
#' @title getGeneM
#' @description Helper function to retrieve the relevant portion of the hdf5 file for a given gene
#' @param obj The object containing the path to the hdf5 file
#' @param gene Which gene to fetch hdf5 indexes
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
#' @param cells Optional parameter to retrieve indexes for a subset of cells
#' @return Returns a data frame listing the start and count for each barcode with reads covering the requested gene
#' @export
#' @examples getGeneM(obj = obj, gene = "RBFOX3", type = "CH", cells = c("ATTGAGGATACCAATTCCATCACGAGCG", "ATTGAGGATAACGCTTCGTCCGCCGATC"))
#' @importFrom dplyr pull filter mutate
#' @importFrom purrr map
#' @importFrom rhdf5 h5read
#' @importFrom tibble rownames_to_column
getGeneM <- function(obj,
                     gene,
                     type,
                     threads = 1,
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

  # get barcodes and paths from amethyst object
  if (is.null(obj@h5paths)) {
    stop("Please generate the path list for each barcode and store in the obj@h5paths slot.")
  } else {
    h5paths <- obj@h5paths[rownames(obj@h5paths) %in% ids, , drop = FALSE]
    barcodes <- as.list(rownames(h5paths))
    paths <- as.list(h5paths$paths)
  }


  # set up multithreading
  if (threads > 1) {
    future::plan(future::multicore, workers = threads)
  }

  # Read in the portion of the hdf5 file with reads covering the given gene for all requested cells
  gene_m <- furrr::future_pmap(
    .l = list(paths, barcodes), .f = function(path, bar) {
      tryCatch({
        dataset_name <- paste0(type, "/", bar)
        start <- indexes |> dplyr::filter(cell_id == bar) |> dplyr::pull(start)
        count <- indexes |> dplyr::filter(cell_id == bar) |> dplyr::pull(count)
        rhdf5::h5read(path, name = dataset_name, start = start, count = count, stride = 1)
      }, error = function(e) {
        cat("Error processing data for barcode", bar, ":", conditionMessage(e), "\n")
        return(NULL)  # Return NULL or any other value indicating failure
      })
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

