# Authors: Lauren Rylaarsdam, PhD
# 2024-2025

############################################################################################################################
#' @title makeWindows
#' @description Calculate methylation levels across fixed genomic windows, bed file coordinates, gene bodies, or promoter regions
#'
#' @param obj Object for which to calculate methylation windows
#' @param type What type of methylation to retrieve; e.g. "CH" or "CG"
#' @param genes Optional genes to calculate % methylation over body
#' @param stepsize If specified, the methylation levels will be calculated over genomic windows of the indicated size
#' @param bed Optional bed file to input. If specified, methylation levels will be calculated over input regions
#' @param metric Calculate either methylation "percent" (out of 100), "score" (normalized metric from -1 to 1 with
#' 0 as the mean), or "ratio" (value / global methylation)
#' @param groupBy Optional metadata parameter to group cells by. If included, the function will calculate mean methylation of a given window over the group
#' @param threads Enables multithreading
#' @param index If calculating genomic windows, specify the name of the chr index in the index slot.
#' This index should contain the coordinates in the hdf5 file corresponding to each chromosome. Reduces memory constraints.
#' @param futureType Method of parallelization, i.e. "multicore" (recommended) or "multisession". Multisession is more
#' memory-efficient, but make sure your workspace is as small as possible before using.
#' @param nmin Minimum number of observations for the window to be included
#' @param promoter Boolean; option to calculate values over the predicted promoter region of each gene.
#' Promoter is defined as strand-aware start site +/- 1500 bp. This approach will not be appropriate for many genes
#' @param chrList Optional; chromosome whitelist (only use if necessary)
#' @param save Boolean indicating whether to save the intermediate output by chromosome. Good for very large datasets
#' @return Returns a data frame with columns as cells, rows as genomic windows, and values as aggregated methylation levels
#' @export
#' @importFrom data.table data.table tstrsplit := setorder setDT rbindlist copy
#' @importFrom dplyr pull filter select mutate
#' @importFrom furrr future_map
#' @importFrom future plan multicore sequential
#' @importFrom plyr round_any
#' @importFrom rhdf5 h5ls h5read
#' @importFrom tibble column_to_rownames rownames_to_column
#' @examples
#' \dontrun{
#'   obj@genomeMatrices[["cg_100k_score"]] <- makeWindows(obj, stepsize = 100000, type = "CG", metric = "score")
#'   obj@genomeMatrices[["pbmc_vmrs"]] <- makeWindows(obj, bed = "~/Downloads/pbmc_vmr.bed", type = "CG", metric = "percent", threads = 10, index = "chr_cg", nmin = 2)
#'   obj@genomeMatrices[["cg_promoters"]] <- makeWindows(obj, genes = protein_coding, promoter = TRUE, type = "CG", metric = "percent", index = "chr_cg", nmin = 5)
#' }

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
    save = FALSE,
    chrList = NULL) {

  # check appropriate metric was specified
  if (!metric %in% c("percent", "score", "ratio")) {
    stop("Metric calculated must be one of either percent, score, or ratio.")
  }

  # check parameters
  if (sum(!is.null(bed), !is.null(stepsize), !is.null(genes)) > 1) {
    stop("Please only specify a fixed step size, gene list, or input bed file.")
  }

  # check type compatibility
  if (type != "CG" && type != "CH" && metric != "percent") {
    cat("To calculate methylation score or ratio in a non-standard context,
        \nplease check that a corresponding context_pct column for each cell is in the metadata.
        \nFor example: cac_pct for mCAC.")
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

  by_chr <- list()

  if (!is.null(bed)) {
    file_name <- sub(".*/(.*)\\.bed$", "\\1", as.character(bed))

    if (sum(class(bed) %in% "character") > 0) {
      bed <- data.table::data.table(read.table(bed))
    } else if (sum(class(bed) %in% c("data.frame", "data.table")) > 0) {
      data.table::setDT(bed)
    }
    if (ncol(bed) > 3) {
      stop("\nPlease make sure the bed file consists of three columns - chr, start, end - with no header.\n")
    }
    setnames(bed, c("chr", "start", "end"))

    if (is.null(chrList)) {
      chr_groups <- as.list(unique(as.character(bed$chr)))
    } else {
      chr_groups <- chrList
    }

    bed <- split(bed, bed$chr)

    for (chr in chr_groups) {

      # get index positions
      sites <- obj@index[[index]][[chr]] # get chr index for h5 file

      # get paths
      barcodes <- as.list(obj@h5paths[obj@h5paths$barcode %in% sites$cell_id, , drop = FALSE]$barcode)
      paths <- as.list(obj@h5paths[obj@h5paths$barcode %in% sites$cell_id, , drop = FALSE]$path)
      if (!is.null(obj@h5paths$prefix)) {
        prefixes <- as.list(obj@h5paths[obj@h5paths$barcode %in% sites$cell_id, , drop = FALSE]$prefix)
      } else {
        prefixes <- as.list(rep("", length(barcodes)))
      }

      # add up sum c and sum t in member cells
      windows <- furrr::future_pmap(.l = list(paths, barcodes, prefixes), .f = function(path, barcode, prefix) {
        tryCatch({
          bed_tmp <- data.table::copy(bed[[chr]])
          barcode_name <- paste0(prefix, sub("\\..*$", "", barcode))

          if (metric != "percent") {
            meth_cell <- obj@metadata[barcode, paste0("m", tolower(type), "_pct")]/100 # pull global methylation level from metadata
          }

          h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode, "/1"),
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
      if (save) {
        saveRDS(windows_merged, paste0("tmp_", type, "_", metric, "_", chr, "_", file_name, "_nmin", nmin, ".RData"))
      }
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

    #define chr groups
    if (is.null(chrList)) {
      chr_groups <- as.list(names(obj@index[[index]]))
    } else {
      chr_groups <- chrList
    }

    for (chr in chr_groups) {

      # get index positions
      sites <- obj@index[[index]][[chr]] # get chr index for h5 file

      # get paths
      barcodes <- as.list(obj@h5paths[obj@h5paths$barcode %in% sites$cell_id, , drop = FALSE]$barcode)
      paths <- as.list(obj@h5paths[obj@h5paths$barcode %in% sites$cell_id, , drop = FALSE]$path)
      if (!is.null(obj@h5paths$prefix)) {
        prefixes <- as.list(obj@h5paths[obj@h5paths$barcode %in% sites$cell_id, , drop = FALSE]$prefix)
      } else {
        prefixes <- as.list(rep("", length(barcodes)))
      }

      windows <- furrr::future_pmap(.l = list(paths, barcodes, prefixes), .f = function(path, barcode, prefix) {
        if (futureType == "multisession") {
          options(scipen = 999)
        }
        tryCatch({
          barcode_name <- paste0(prefix, sub("\\..*$", "", barcode))

          h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode, "/1"),
                                                     start = sites$start[sites$cell_id == barcode_name],
                                                     count = sites$count[sites$cell_id == barcode_name])) # read in 1 chr at a time
          h5 <- h5[pos %% stepsize == 0, pos := pos + 1] # otherwise sites exactly divisible by stepsize will be their own window
          h5 <- h5[, window := paste0(chr, "_", plyr::round_any(pos, stepsize, floor), "_", plyr::round_any(pos, stepsize, ceiling))]
          if (metric != "percent") {
            meth_cell <- obj@metadata[barcode, paste0("m", tolower(type), "_pct")]/100 # pull global methylation level from metadata
          }
          meth_window <- h5[, .(value = round(sum(c != 0) / (sum(c != 0) + sum(t != 0)), 3), n = sum(c + t, na.rm = TRUE)), by = window]
          meth_window <- meth_window[n >= nmin, .(window, value)]

          if (metric == "percent") { summary <- meth_window[, value := (value * 100)] }
          if (metric == "score") { summary <- meth_window[, value := round((ifelse((value - meth_cell > 0),
                                                                                   (value - meth_cell)/(1 - meth_cell),
                                                                                   (value - meth_cell)/meth_cell)), 3)] }
          if (metric == "ratio") { summary <- meth_window[, value := round(value / meth_cell, 3)] }
          data.table::setnames(summary, "value", barcode_name)
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

    if (is.null(chrList)) {
      chr_groups <- as.list(unique(as.character(bed$chr)))
    } else {
      chr_groups <- chrList
    }

    bed <- split(bed, bed$chr)

    for (chr in chr_groups) {

      # get index positions
      sites <- obj@index[[index]][[chr]] # get chr index for h5 file

      # get paths
      barcodes <- as.list(obj@h5paths[obj@h5paths$barcode %in% sites$cell_id, , drop = FALSE]$barcode)
      paths <- as.list(obj@h5paths[obj@h5paths$barcode %in% sites$cell_id, , drop = FALSE]$path)
      if (!is.null(obj@h5paths$prefix)) {
        prefixes <- as.list(obj@h5paths[obj@h5paths$barcode %in% sites$cell_id, , drop = FALSE]$prefix)
      } else {
        prefixes <- as.list(rep("", length(barcodes)))
      }

      # add up sum c and sum t in member cells
      windows <- furrr::future_pmap(.l = list(paths, barcodes, prefixes), .f = function(path, barcode, prefix) {
        tryCatch({
          bed_tmp <- data.table::copy(bed[[chr]])
          barcode_name <- paste0(prefix, sub("\\..*$", "", barcode))
          if (metric != "percent") {
            meth_cell <- obj@metadata[barcode, paste0("m", tolower(type), "_pct")]/100 # pull global methylation level from metadata
          }
          h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode, "/1"),
                                                     start = sites$start[sites$cell_id == barcode],
                                                     count = sites$count[sites$cell_id == barcode])) # read in 1 chr at a time

          meth_gene <- h5[bed_tmp, on = .(chr = chr, pos >= start, pos <= end), nomatch = 0L, .(chr, start, end, pos = x.pos, c, t, gene)]
          meth_gene <- meth_gene[, .(value = round((sum(c != 0)/ (sum(c != 0) + sum(t != 0))), 4),
                                     n = sum(c + t, na.rm = TRUE)), by = gene]
          meth_gene <- meth_gene[n >= nmin, ][, n := NULL]

          if (metric == "percent") { summary <- meth_gene[, value := (value * 100)] }
          if (metric == "score") { summary <- meth_gene[, value := round((ifelse((value - meth_cell) > 0,
                                                                                 (value - meth_cell)/(1 - meth_cell),
                                                                                 (value - meth_cell)/meth_cell)), 3)] }
          if (metric == "ratio") { summary <- meth_gene[, value := round(value / meth_cell, 3)] }
          data.table::setnames(summary, "value", barcode_name)
        }, error = function(e) {
          cat("Error processing data for barcode", barcode, ":", conditionMessage(e), "\n")
          return(NULL)  # Return NULL or any other value indicating failure
        })
      }, .progress = TRUE)

      # Reduce merge the data.tables per cell in chunks (way faster for large datasets)
      windows <- split(windows, ceiling(seq_along(windows)/1000))
      windows_merged <- Reduce(function(x, y) merge(x, y, by = c("gene"), all = TRUE, sort = FALSE),
                               furrr::future_map(windows, ~ Reduce(function(x, y) merge(x, y, by = c("gene"), all = TRUE, sort = FALSE), .x), .progress = TRUE))
      if (save) {
        saveRDS(windows_merged, paste0("tmp_", type, "_", metric, "_", chr, "_", length(genes), "genes_nmin", nmin, ".RData"))
      }
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
#' @importFrom dplyr left_join
#' @importFrom tibble column_to_rownames rownames_to_column
#' @examples
#' \dontrun{
#'   obj <- addMetadata(obj = obj, metadata = annotations)
#' }
addMetadata <- function(obj,
                        metadata) {

  if (!is.null(obj@h5paths$prefix)) {
    cat("\nWarning: Make sure your prefixes are already appended to the metadata you are adding.")
  }

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
#' @param names Column name(s) to add in the metadata
#' @return Updates obj to contain the new metadata concatenated to the old metadata
#' @export
#' @importFrom dplyr left_join
#' @importFrom tibble column_to_rownames rownames_to_column
#' @examples
#' \dontrun{
#'   obj <- addMetadata(obj = obj, metadata = annotations, names = c("plate"))
#' }
addAnnot <- function(obj,
                     file,
                     names) {

  if (!is.null(obj@h5paths$prefix)) {
    cat("\nWarning: Make sure your prefixes are already appended to the metadata you are adding.")
  }

  annot <- utils::read.table(file, sep = "\t", header = F)
  colnames(annot) <- c("cell_id", names)
  if (is.null(obj@metadata)) {
    obj@metadata <- annot |> tibble::column_to_rownames(var = "cell_id")
  } else {
    obj <- addMetadata(obj, metadata = annot |> tibble::column_to_rownames(var = "cell_id"))
  }
  output <- obj
}

############################################################################################################################
#' @title addCellInfo
#' @description Helper function to the information contained in the Premethyst cellInfo file intermediate output to the obj@metadata slot
#'
#' @param obj Amethyst object to add cell info metadata
#' @param file Path of the cellInfo.txt file
#' @return Returns the same amethyst object with information contained in the cellInfo file added to the obj@metadata slot
#' @export
#' @importFrom tibble column_to_rownames
#' @importFrom utils read.table
#' @examples
#' \dontrun{
#'   obj <- addCellInfo(obj = obj, file = "~/Downloads/cellInfo.txt")
#' }

addCellInfo <- function(obj,
                        file) {
  cat("This is a helper function intended to process the Premethyst pipeline cellInfo.txt file intermediate output. If you are using the ScaleMethyl pipeline and used the createScaleObject helper function, metadata should be added automatically. If neither of these pipelines apply, please add metadata to slot manually where rownames are cell barcodes.\n")

  if (!is.null(obj@h5paths$prefix)) {
    cat("Warning: Make sure your prefixes are already appended to the metadata you are adding.\n")
  }

  cellInfo <- utils::read.table(file, sep = "\t", header = F)
  if (ncol(cellInfo) == 6) {
    colnames(cellInfo) <- c("cell_id", "cov", "cg_cov", "mcg_pct", "ch_cov", "mch_pct")
  } else if (ncol(cellInfo) == 10) {
    colnames(cellInfo) <- c("cell_id", "cov", "cg_cov", "mcg_pct", "ch_cov", "mch_pct", "n_frag", "n_pair", "n_single", "xpct_m")
  } else if (ncol(cellInfo) != 6 && ncol(cellInfo) != 10) {
    cat("The cell file is not in the expected format of 6 or 10 columns. Columns after 10 will not be named.\n")
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
#' Warning: we recommend using with caution. Ideally, any results would also be visible without imputation, and this
#' would be used to help smooth values.
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
    print(paste0("Warning: Replacing ", num_na, " NAs (", pct_na, "% of matrix) to perform imputation.\n"))
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
    stop("Default NA replacement must be 0, 1, 'mch_pct', or 'mcg_pct'.\n")
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
#' @param nmin Minimum number of cells in group to include in the result
#' @param groupBy Parameter in the metadata to aggregate over
#' @return Returns a matrix containing the mean value per group
#' @export
#' @importFrom data.table setDT tstrsplit
#' @importFrom dplyr select filter
#' @importFrom tibble rownames_to_column column_to_rownames
#' @examples
#' \dontrun{
#'   obj@genomeMatrices[["ch_100k_pct_cluster"]] <- aggregateMatrix(obj, matrix = "ch_100k_pct", groupBy = "cluster_id", nmin = 10)
#' }

aggregateMatrix <- function(obj,
                            matrix,
                            groupBy,
                            nmin = 5) {

  matrix <- obj@genomeMatrices[[matrix]]
  membership <- obj@metadata |> dplyr::select(groupBy)
  colnames(membership) <- "membership"
  groups <- as.list(unique(membership$membership))
  aggregated <- list()
  num_members <- membership |> dplyr::group_by(membership) |> dplyr::summarise(n = n())

  for (i in groups) {
    # Check if n is greater than or equal to nmin
    if (num_members[num_members$membership == i, "n"] >= nmin) {
      if (num_members[num_members$membership == i, "n"] > 1) {
        members <- rownames(membership |> dplyr::filter(membership == i))
        aggregated[[i]] <- as.data.frame(rowMeans(matrix[, members], na.rm = TRUE))
        colnames(aggregated[[i]]) <- paste0(i)
        aggregated[[i]] <- aggregated[[i]] |> tibble::rownames_to_column(var = "window")
      } else if (num_members[num_members$membership == i, "n"] == 1) {
        members <- rownames(membership |> dplyr::filter(membership == i))
        aggregated[[i]] <- as.data.frame(matrix[, members, drop = FALSE]) |> tibble::rownames_to_column(var = "window")
      }
    } else {
      # Optionally, you can add a message or log when skipping a group.
      message("Skipping group ", i, " because n < nmin.")
    }
  }
  aggregated <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE), aggregated)

  # reorder if genomic windows
  if (all(grepl("chr", aggregated$window) & grepl("_", aggregated$window))) {
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
#' @examples
#' \dontrun{
#'   subset <- subsetObject(obj = obj, cells = c("ATTGAGGATAACGCTTCGTCCGCCGATC", "ATTGAGGATAACGCTTCGTCTAAGGTCA"))
#' }
subsetObject <- function(obj,
                         cells) {
  subset <- obj

  subset@h5paths$cell_ids <- paste0(subset@h5paths$prefix, subset@h5paths$barcode)
  subset@h5paths <- subset@h5paths[subset@h5paths$cell_ids %in% cells, , drop = FALSE]

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
#' @param objList List of amethyst objects to merge. Make sure barcode names are unique
#' @param genomeMatrices List of genomeMatrices to merge from each object
#'
#' @return Returns a new object with merged slots. Dimensionality reductions should be re-run post merging.
#' @export
#' @importFrom tibble column_to_rownames rownames_to_column
#' @examples
#' \dontrun{
#'   new <- combineObject(objList = list(obj1, obj2, obj3), genomeMatrices = c("ch_100k_pct", "gene_ch"))
#' }
combineObject <- function(objList,
                          genomeMatrices) {
  cat("\nMake sure barcodes are unique when combining objects. If not, please see our vignette on how to merge projects with overlapping barcodes.")

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
#' @title convertObject
#' @description Shortcut for converting old amethyst object format (v0.0.0.9000) to new (v1.0.0)
#' Adds "tracks" and "results" slots, converts obj@ref to a data.table, moves tracks and umap/tsne
#' coordinates to correct locations if the names contain "tracks" and "umap" or "tsne", respectively
#'
#' @param obj Old amethyst object to convert
#' @return Returns an updated amethyst object with slots:
#' h5paths, genomeMatrices, tracks, reductions, index, metadata, ref, results
#' @export
#' @importFrom data.table data.table
#' @examples
#' \dontrun{
#'  new_obj <- convertObject(old_obj)
#'}

convertObject <- function(obj) {

  methods::setClass("amethyst", slots = c(
    h5paths = "ANY",
    genomeMatrices = "ANY",
    tracks = "ANY",
    reductions = "ANY",
    index = "ANY",
    metadata = "ANY",
    ref = "ANY",
    results = "ANY"
  ))

  paths <- obj@h5paths |> dplyr::rename("path" = "paths")
  paths$barcode <- rownames(obj@h5paths)
  matrices <- obj@genomeMatrices[-(grep(pattern = "tracks", x = names(obj@genomeMatrices)))]
  tracks <- obj@genomeMatrices[(grep(pattern = "tracks", x = names(obj@genomeMatrices)))]
  ref <- data.table::data.table(obj@ref)

  new_obj <- methods::new(Class = "amethyst",
               h5paths = obj@h5paths,
               genomeMatrices = matrices,
               tracks = tracks,
               reductions = obj@reductions,
               index =  obj@index,
               metadata = obj@metadata,
               ref = ref,
               results = NULL)

  if (any(grepl("umap", colnames(obj@metadata)))) {
    new_obj@reductions[["umap"]] <- obj@metadata[, c("umap_x", "umap_y")]
    colnames(new_obj@reductions[["umap"]]) <- c("dim_x", "dim_y")
  }

  if (any(grepl("tsne", colnames(obj@metadata)))) {
    new_obj@reductions[["tsne"]] <- obj@metadata[, c("tsne_x", "tsne_y")]
    colnames(new_obj@reductions[["tsne"]]) <- c("dim_x", "dim_y")
  }

  return(new_obj)

}

############################################################################################################################
#' @title getGeneCoords
#' @description Helper function to fetch the appended chromosome, start, and end locations for a given gene

#' @param ref Genome annotation file
#' @param gene Name of gene for which to fetch coordinates
#' @return Returns the appended chromosome, start, and end locations for a given gene; e.g. "chr2_199269505_199471266"
#' @export
#' @examples
#' \dontrun{
#'  getGeneCoords(ref = obj@ref, gene = "SATB2")
#'}
getGeneCoords <- function(
    ref,
    gene) {

  coords <- ref$location[ref$type == "gene" & ref$gene_name == gene]
  return(coords)

}

