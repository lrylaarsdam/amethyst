#' @title Calculate Window Methylation Statistic
#' @description For a given genomic region, calculates a specific statistic
#' (such as percent methylated).
#' @param h5path Path to the .h5 file containing the methylation data for each
#' cell
#' @param barcode The barcode for a specific cell within the .h5 file
#' @param stepszie The size of the windows
#' @param stat One of `'percent'`, `'ration'` or `'score'`; the statistic to
#' calculate for the window
#' @returns Data.table of the cell methylation data updated to include the stat
calcWindowMethylationStat <- function(barcode, h5path, stepsize, stat) {
  .h5 <- data.table::data.table(rhdf5::h5read(h5path, name = paste0(type, "/", barcode)))
  .h5 <- .h5[, window := paste0(
    chr,
    "_",
    roundAny(pos, stepsize, floor),
    "_",
    roundAny(pos, stepsize, ceiling)
  )]

  if(stat == 'percent') {
    .h5 <- .h5[, .(
      value = round(sum(c != 0) * 100 / (sum(c != 0) + sum(t != 0)), 1)
    ), by = window]
    return(.h5)
  }
  .h5 <- .h5[, glob_m := sum(c != 0) / (sum(c != 0) + sum(t != 0))]
  if(stat == 'ratio') {
    .h5 <- .h5[, .(
        value = round((sum(c != 0) / (sum(c != 0) + sum(t != 0))) / mean(glob_m), 3)
      ), by = window]
  }
  if(stat == 'score') {
    .h5 <- .h5[, .(
      meth = sum(c != 0) / (sum(c != 0) + sum(t != 0)), glob_m = mean(glob_m)
    ), by = window]
    .h5 <- .h5[, diff := Map(`-`, meth, glob_m)]
    .h5 <- .h5[, value := round(diff / glob_m, 3)][diff > 0, value := round(diff / (1 - glob_m), 3)][]
  }
  .h5
}


calcWindowMethylationScore <- function(barcode, h5path, stepsize) {
  calcWindowMethylationStat(barcode, h5path, stepsize, stat = 'score')
}


calcWindowMethylationRatio <- function(barcode, h5path, stepsize) {
  calcWindowMethylationStat(barcode, h5path, stepsize, stat = 'ratio')
}


calcWindowMethylationPercent <- function(barcode, h5path, stepsize) {
  calcWindowMethylationStat(barcode, h5path, stepsize, stat = 'percent')
}


#' @title Round any
#' @description
#' Local implementation of the plyr function `round_any`.
#' @param x numeric vector to be rounded
#' @param accuracy number to round to
#' @param f rounding function
#' @returns x with numbers rounded to `accuracy` number of significant figures
#' as defined by `f`
#' @examples
#' roundAny(111.1111, 1)
#' roundAny(111.1111, 10)
#' roundAny(111.1111, 100)
#' roundAny(111.1111, 0.1)
#' roundAny(111.1111, 1, floor)
#' roundAny(111.1111, 0.1, floor)
#' roundAny(111.1111, 1, ceiling)
#' roundAny(111.1111, 100, ceiling)
#' roundAny(111.1111, 0.1, ceiling)
roundAny <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}


#' @title Get gene M
#' @description Helper function to retrieve the relevant portion of the hdf5 file for a given gene
#'
#' @param obj The object containing the path to the hdf5 file
#' @param gene Which gene to fetch hdf5 indexes
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
#' @param cells Optional parameter to retrieve indexes for a subset of cells
#'
#' @return Returns a data frame listing the start and count for each barcode with reads covering the requested gene
#' @export
#'
#' @examples
#' obj <-
#'   getGeneM(obj, gene = "RBFOX3", type = "CH")
getGeneM <- function(obj,
                     gene,
                     type,
                     cells = NULL) {
  # Check if the gene has been indexed. Return error if not.
  if (!(gene %in% names(obj@index[[type]]))) {
    print(paste0("Error: ", gene, " has not been indexed."))
  }

  # Make data frame of all indexes for the gene
  indexes <- data.table::data.table(as.data.frame(obj@index[[type]][[gene]]))
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
      rel_indexes <- indexes[cell_id == x, ]
      start <- rel_indexes$start
      count <- rel_indexes$count
      # start <- indexes |>
      #   dplyr::filter(cell_id == x) |>
      #   dplyr::pull(start)
      #
      # count <- indexes |>
      #   dplyr::filter(cell_id == x) |>
      #   dplyr::pull(count)
      out <- rhdf5::h5read(obj@h5path,
        name = dataset_name,
        start = start,
        count = count,
        stride = 1
      )
    }
  )

  # The result will be a list. Name the list with the corresponding cell barcodes.
  names(gene_m) <- ids

  # Row bind the list into one big list
  gene_m <- do.call(rbind, gene_m)
  # gene_m <- data.table::rbindlist(gene_m)

  # Make a new column containing the barcodes and remove unnecessary characters added automatically by R
  gene_m <- gene_m |>
    tibble::rownames_to_column(var = "cell_id") |>
    dplyr::mutate(cell_id = sub("\\..*$", "", cell_id))
  output <- gene_m
  output
}


#' @title Make windows
#' @description Calculate percent methylation across a genomic window of given
#' size or a gene body. If desired, aggregate average percent methylation
#' according to a categorical value in the metadata.
#' @param obj Object for which to calculate percent methylation windows
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene
#' body mCG, or promoter mCG.
#' @param groupBy Optional metadata parameter to group cells by. If included,
#' the function will calculate mean percent methylation of a given window over
#' the group.
#' @param genes If specified, the % methylation will be calculated over the
#' indicated gene bodies.
#' @param stepsize If specified, the % methylation will be calculated over
#' genomic windows of the indicated size.
#' @param nmin Optional. Minimum number of cytosines recorded in the window/gene
#' to be included in the output.
#'
#' @return Returns a data frame with the genomic region or gene name as column
#' names, cell barcode or group as row names, and values as percent methylated.
#' @export
#'
#' @examples
#' \dontrun{
#' makeWindows(obj, groupBy = cluster_id, genes = c("TBR2", "TBR1"), nmin = 10, type = "CH")
#' makeWindows(obj, stepsize = 100000, type = "CH")
#' }
makeWindows <- function(
    obj,
    type,
    groupBy = NULL,
    genes = NULL,
    stepsize = FALSE,
    nmin = NULL,
    metric = "percent",
    threads = 1,
    bed = NULL) {
  if (threads > 1) {
    future::plan(future::multicore, workers = threads)
  }

  # If stepsize is specified, make fixed size genomic windows
  if (stepsize) {
    # Optional aggregation over a categorical variable in obj@metadata
    if (!is.null(groupBy)) {
      groups <- as.list(unique(obj@metadata[[groupBy]]))
      aggregated <- list()
      for (i in groups) {
        # Retrieve barcodes for the cells belonging to the appropriate metadata
        barcodes <- rownames(obj@metadata |> dplyr::filter(.data[[groupBy]] == i))
        if (metric == "percent") {
          # read in each barcode and calculate % methylation over 100kb genomic
          # windows for each cell
          windows <- furrr::future_map(
            .x = barcodes,
            .f = calcWindowMethylationPercent,
            h5path = obj@h5path,
            stepsize = stepsize,
            .progress = FALSE
          )
        } else if (metric == "ratio") {
          # read in each barcode and calculate % methylation over 100kb genomic
          # windows for each cell
          windows <- furrr::future_map(
            .x = barcodes,
            .f = calcWindowMethylationRatio,
            h5path = obj@h5path,
            stepsize = stepsize,
            .progress = FALSE
          )
        } else if (metric == "score") {
          windows <- furrr::future_map(
            .x = barcodes,
            .f = calcWindowMethylationScore,
            h5path = obj@h5path,
            stepsize = stepsize,
            .progress = FALSE
          )
        }
        windows <- dplyr::bind_rows(windows)

        # summarise average methylation across the group for each window
        windows <- windows |>
          dplyr::group_by(window) |>
          dplyr::summarise(value = mean(value))

        # change the format so genomic windows are columns
        windows <- tidyr::pivot_wider(windows, names_from = window, values_from = value)

        # set the row name to be the group
        rownames(windows) <- paste(as.character(i))

        # add windows aggregated across the group to a list
        aggregated[[i]] <- windows
      }

      # aggregate the data frame after windows for each group are calculated.
      # Columns will be genomic windows; rows will be the group name.
      windows <- dplyr::bind_rows(aggregated)
    } else {
      barcodes <- rhdf5::h5ls(obj@h5path)
      # extract unique barcodes from the h5 file
      barcodes <- as.list(
        unique(
          barcodes |>
            dplyr::filter(otype == "H5I_DATASET") |>
            dplyr::pull(name)
        )
      )
      if (metric == "percent") {
        # read in each barcode and calculate % methylation over 100kb genomic
        # windows for each cell
        windows <- furrr::future_map(
          .x = barcodes,
          .f = calcWindowMethylationPercent,
          h5path = obj@h5path,
          stepsize = stepsize,
          .progress = FALSE
        )
      } else if (metric == "ratio") {
        # read in each barcode and calculate % methylation over 100kb genomic
        # windows for each cell
        windows <- furrr::future_map(
          .x = barcodes,
          .f = calcWindowMethylationRatio,
          h5path = obj@h5path,
          stepsize = stepsize,
          .progress = FALSE
        )
      } else if (metric == "score") {
        # read in each barcode and calculate % methylation over 100kb genomic
        # windows for each cell
        windows <- furrr::future_map(
          .x = barcodes,
          .f = calcWindowMethylationScore,
          h5path = obj@h5path,
          stepsize = stepsize,
          .progress = FALSE
        )
      }
      names(windows) <- barcodes
      windows <- dplyr::bind_rows(windows, .id = "cell_id")
      windows <- tidyr::pivot_wider(windows,
        id_cols = cell_id,
        names_from = window,
        values_from = value
      ) |>
        # change format so columns are genomic windows and rows are cells
        tibble::column_to_rownames(var = "cell_id")
    }
  }

  # If genes are specified, calculate percent methylation over the indicated gene bodies
  if (!is.null(genes)) { # based on protein coding genes from gtf. If genes with the same name are present, the longest is taken.
    # Get genomic locations (chr, start, and end) from the reference
    ranges <- obj@ref |>
      dplyr::filter(type == "gene") |>
      dplyr::mutate(length = end - start) |>
      dplyr::group_by(gene_name) |>
      dplyr::mutate(replicate = n()) |>
      dplyr::group_by(gene_name) |>
      dplyr::arrange(desc(length)) |>
      dplyr::filter(row_number() == 1) |>
      dplyr::select("seqid", "start", "end", "gene_name", "location") |>
      dplyr::arrange(seqid, start) |>
      dplyr::filter(gene_name %in% genes)
    if (!is.null(groupBy)) {
      groups <- as.list(unique(obj@metadata[[groupBy]]))
      aggregated <- list()
      for (i in groups) {
        barcodes <- rownames(
          obj@metadata |>
            dplyr::filter(.data[[groupBy]] == i)
        )
        windows <- furrr::future_map(
          .x = ranges$gene_name,
          .f = function(x) {
            getGeneM(obj = {{ obj }}, gene = x, type = {{ type }})
          }
        )
        names(windows) <- ranges$gene_name
        windows <- dplyr::bind_rows(windows, .id = "gene")
        windows <- windows |>
          dplyr::group_by(gene, cell_id) |>
          dplyr::summarise(n = n(), pct_m = sum(c != 0) * 100 / (sum(c != 0) + sum(t != 0)), .groups = "keep")
        if (!is.null(nmin)) {
          windows <- windows |> dplyr::filter(n >= nmin)
        }
        windows <- windows |>
          dplyr::group_by(gene) |>
          dplyr::summarise(pct_m = mean(pct_m))
        windows <- tidyr::pivot_wider(windows, names_from = gene, values_from = pct_m)
        rownames(windows) <- paste(as.character(i))
        aggregated[[i]] <- windows
      }
      windows <- dplyr::bind_rows(aggregated)
    } else {
      windows <- furrr::future_map(
        .x = ranges$gene_name,
        .f = function(x) {
          getGeneM(obj = {{ obj }}, gene = x, type = {{ type }}) |>
            dplyr::group_by(cell_id) |>
            dplyr::summarise(pct_m = sum(c != 0) * 100 / (sum(c != 0) + sum(t != 0)), n = n(), .groups = "keep")
        },
        .progress = TRUE
      )
      names(windows) <- ranges$gene_name
      windows <- dplyr::bind_rows(windows, .id = "gene")
      if (!is.null(nmin)) {
        windows <- windows |> dplyr::filter(n >= nmin)
      }
      windows <- tidyr::pivot_wider(windows,
        id_cols = cell_id,
        names_from = gene,
        values_from = pct_m
      ) |>
        tibble::column_to_rownames(var = "cell_id")
    }
  }

  # Calculate avg windows over bed file regions
  if (!is.null(bed)) {
    ranges <- bed
    colnames(ranges) <- c("chr", "start", "end")
    ranges <- ranges |>
      dplyr::mutate(window = paste0(chr, "_", start, "_", end))

    barcodes <- rhdf5::h5ls(obj@h5path)
    barcodes <- as.list(
      unique(
        barcodes |>
          dplyr::filter(otype == "H5I_DATASET") |>
          dplyr::pull(name)
      )
    )

    for (i in barcodes) {
      h5 <- rhdf5::h5read(obj@h5path, name = paste0(type, "/", i))
      glob_m <- sum(c != 0) * 100 / (sum(c != 0) + sum(t != 0))
      for (j in 1:length(ranges$window)) {
        ## continue editing
      }
    }

    # read in each barcode and calculate % methylation over 100kb genomic
    # windows for each cell
    windows <- furrr::future_map(
      .x = barcodes,
      .f = function(x) {
        rhdf5::h5read(obj@h5path, name = paste0(type, "/", x)) |>
          dplyr::mutate(
            window = paste0(
              chr,
              "_",
              plyr::round_any(pos, stepsize, floor),
              "_",
              plyr::round_any(pos, stepsize, ceiling)
            )
          ) |>
          dplyr::group_by(window) |>
          dplyr::summarise(value = round(sum(c != 0) * 100 / (sum(c != 0) + sum(t != 0)), 1), .groups = "keep")
      },
      .progress = TRUE
    )

    names(windows) <- barcodes
    windows <- dplyr::bind_rows(windows, .id = "cell_id")

    # change format so columns are genomic windows and rows are cells
    windows <- tidyr::pivot_wider(windows,
      id_cols = cell_id,
      names_from = window,
      values_from = value
    ) |> tibble::column_to_rownames(var = "cell_id")
  }

  if (threads > 1) {
    future::plan(NULL)
    gc()
  }

  output <- windows
  output
}

#' @title Add metadata
#' @description Add new cell metadata to the metadata slot of a pre-existing object. Make sure row names are cell IDs.
#' @param obj Object to add metadata to.
#' @param metadata Data frame, with row names as cell IDs, containing the new metadata to be added.
#'
#' @return Updates obj to contain the new metadata concatenated to the old metadata.
#' @export
#'
#' @examples obj <- addMetadata(obj, metadata = cellInfo)
addMetadata <- function(obj,
                        metadata) {
  if (length(intersect(rownames(obj@metadata), rownames(metadata))) == 0) {
    print("No intersection detected; check metadata structure. Make sure row names are cell IDs.")
    output <- obj
  } else {
    obj@metadata <- dplyr::left_join(
      obj@metadata |>
        tibble::rownames_to_column(var = "cell_id"),
      metadata |>
        tibble::rownames_to_column(var = "cell_id"),
      by = "cell_id"
    ) |>
      tibble::column_to_rownames(var = "cell_id")
    output <- obj
  }
  output
}

#' @title Add cell info
#' @description Add the information contained in the cellInfo file outputs to the obj@metadata
#' @param obj Amethyst object to add cell info metadata
#' @param file Path of the cellInfo.txt file
#'
#' @return Returns the same Amethyst object but with cellInfo added
#' @export
#'
#' @examples obj <- addCellInfo(obj, file = "~/Downloads/cellInfo.txt")
addCellInfo <- function(obj,
                        file) {
  cellInfo <- read.table(file, sep = "\t", header = F)
  colnames(cellInfo) <- c("cell_id", "cov", "cg_cov", "mcg_pct", "ch_cov", "mch_pct")
  obj <- addMetadata(obj, metadata = cellInfo |>
    tibble::column_to_rownames(var = "cell_id"))
  output <- obj
  output
}

#' @title Impute
#' @description Impute values with the Rmagic package https://github.com/KrishnaswamyLab/MAGIC
#' @param obj Amethyst object to perform imputation on
#' @param matrix Name of the matrix contained in the genomeMatrices slot to perform imputation on
#' @param npca Number of principle components to use
#'
#' @return Adds a new data frame to the genomeMatrices slot with imputed values
#' @export
#'
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
  matrix_imputed <- Rmagic::magic(matrix, npca = npca)[["result"]]
  obj@genomeMatrices[[paste0(name, "_imputed")]] <- matrix_imputed
  output <- obj
  output
}

#' @title Aggregate matrix
#'
#' @param obj Amethyst object containing the matrix to be aggregated
#' @param matrix Name of matrix contained in the genomeMatrices slot to aggregate
#' @param groupBy Parameter in the metadata to groupBy
#' @param name Name of output matrix to store
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- aggregateMatrix(obj, matrix = "ch_100k_pct", groupBy = cluster_id, name = "cluster_ch_100k_pct")
#' }
aggregateMatrix <- function(obj,
                            matrix,
                            groupBy,
                            name) {
  matrix <- obj@genomeMatrices[[matrix]]
  metadata <- obj@metadata |>
    dplyr::select({{ groupBy }})
  matrix <- merge(metadata, matrix, by = 0) |>
    tibble::column_to_rownames(var = "Row.names")
  matrix <- matrix |>
    dplyr::group_by({{ groupBy }}) |>
    dplyr::summarise(
      dplyr::across(
        where(is.numeric),
        mean,
        na.rm = TRUE
      ),
      .groups = "keep"
    ) |>
    tibble::column_to_rownames(var = (rlang::as_name(dplyr::ensym(groupBy))))
  obj@genomeMatrices[[name]] <- matrix
  output <- obj
  output
}
