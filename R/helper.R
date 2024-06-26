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
    .l = list(paths, barcodes),
    .f = function(path, bar) {
      dataset_name <- paste0(type, "/", bar)
      start <- indexes |> dplyr::filter(cell_id == bar) |> dplyr::pull(start)
      count <- indexes |> dplyr::filter(cell_id == bar) |> dplyr::pull(count)
      rhdf5::h5read(path, name = dataset_name, start = start, count = count, stride = 1)
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
#' @param index If calculating genomic windows, specify the name of the chr index in the index slot.
#' This index should contain the coordinates in the hdf5 file corresponding to each chromosome.
#' This reduces memory constraints by processing one chromosome at a time.
#' @param species "hg38" or "mm10". This determines which chromosome list is used.
#' @param futureType Method of parallelization, i.e. "multicore" (recommended) or "multisession". Multisession is more
#' memory-efficient, but make sure your workspace is as small as possible before using.
#' @param nmin Minimum number of observations for the window to be included.
#'
#' @return Returns a data frame of the calculated methylation windows with the genomic region or gene name as rows and cell barcode or group as columns
#' @export
#' @examples makeWindows(obj = obj, groupBy = "cluster_id", genes = c("SATB2", "GAD1"), type = "CH", metric = "percent", threads = 10)
#' @examples makeWindows(obj = obj, stepsize = 100000, type = "CG", metric = "score", species = "hg38", nmin = 10, index = "chr_cg")
#' @importFrom data.table data.table tstrsplit := setorder
#' @importFrom dplyr pull filter select mutate
#' @importFrom furrr future_map
#' @importFrom future plan multicore sequential
#' @importFrom plyr round_any
#' @importFrom rhdf5 h5ls h5read
#' @importFrom tibble column_to_rownames rownames_to_column
makeWindows <- function(
    obj,
    type,
    metric = "percent",
    index = NULL,
    species = NULL,
    groupBy = NULL,
    genes = NULL,
    stepsize = NULL,
    threads = 1,
    futureType = "multicore",
    nmin = 5,
    save = FALSE) {

  # check appropriate metric was specified
  if (!metric %in% c("percent", "score", "ratio")) {
    stop("Metric calculated must be one of either percent, score, or ratio.")
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

  # get barcodes and paths from amethyst object
  if (is.null(obj@h5paths)) {
    stop("Please generate the path list for each barcode and store in the obj@h5paths slot.")
  } else {
    barcodes <- as.list(rownames(obj@h5paths))
    paths <- as.list(obj@h5paths$paths)
  }

  # If stepsize is specified, make fixed size genomic windows
  if (!is.null(stepsize)) {

    #define chr groups
    if (species == "hg38") {
      chr_groups <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
    } else if (species == "mm10") {
      chr_groups <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX")
    } else if (!(species %in% c("hg38", "mm10"))) {
      stop("Species options are only hg38 or mm10 at this time.")
    }

    # check index
    if (is.null(obj@index[[index]])) {
      stop("Please index which rows in the h5 files correspond to each chromosome using indexChr.")
    }

    by_chr <- list()
    for (j in 1:length(chr_groups)) {

      # get index positions
      sites <- obj@index[[index]][[chr_groups[[j]]]] # get chr index for h5 file

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
          meth_cell <- h5[, round(sum(c != 0, na.rm = TRUE) / (sum(c != 0, na.rm = TRUE) + sum(t != 0, na.rm = TRUE)), 4)] # determine global methylation level
          meth_window <- h5[, .(value = round(sum(c != 0, na.rm = TRUE) / (sum(c != 0, na.rm = TRUE) + sum(t != 0, na.rm = TRUE)), 3), n = sum(c + t, na.rm = TRUE)), by = window]
          meth_window <- meth_window[n > nmin, .(window, value)]

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
        saveRDS(windows_merged, paste0("tmp_", type, "_", chr_groups[[j]], "_", stepsize, "kb_windows_nmin", nmin, ".RData"))
      }
      by_chr[[j]] <- windows_merged
      cat("\nCompleted ", chr_groups[[j]],"\n")
    }

    gc()
    windows_merged <- do.call(rbind, by_chr)

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
    windows_merged <- windows_merged[!sapply(rownames(windows_merged), function(name) length(strsplit(name, "_")[[1]]) > 3 || grepl("chrX|chrY|chrEBV|chrM|KI", name)), ]

  }

  if (!is.null(genes)) {
    windows <- furrr::future_map(as.list(genes), function(x) {
      genem <- data.table::data.table(getGeneM(obj = {{obj}}, gene = x, type = {{type}}))
      genem <- genem[, n := c + t, by = cell_id]
      genem <- genem[n > nmin]

      if (metric == "percent") {
        summary <- genem[, .(value = round(sum(c != 0, na.rm = TRUE) * 100 / (sum(c != 0, na.rm = TRUE) + sum(t != 0, na.rm = TRUE)), 2)), by = cell_id]
      }
      if (metric == "score" | metric == "ratio") {
        if (type == "CG" | type == "promoters") {
          glob_m <- obj@metadata |> dplyr::mutate(pctm = mcg_pct/100) |> dplyr::select(pctm) |> tibble::rownames_to_column(var = "cell_id")
        } else if (type == "CH") {
          glob_m <- obj@metadata |> dplyr::mutate(pctm = mch_pct/100) |> dplyr::select(pctm) |> tibble::rownames_to_column(var = "cell_id")
        }
        summary <- genem[, .(value = round(sum(c != 0, na.rm = TRUE) / (sum(c != 0, na.rm = TRUE) + sum(t != 0, na.rm = TRUE)), 2)), by = cell_id]
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
    }, .progress = TRUE)
    windows_merged <- Reduce(function(x, y) merge(x, y, by = "cell_id", all = TRUE, sort = FALSE), windows) |> tibble::column_to_rownames(var = "cell_id")
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
#' @param npca Number of principle components to use
#' @return Adds a new data frame to the genomeMatrices slot with imputed values
#' @export
#' @examples obj <- impute(obj, matrix = "gene_ch")
#' @importFrom Rmagic magic
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

######################################################################
#' @title makeClusterTracks
#' @description Calculate mean percent methylation over 1500 base pair genomic windows shifted every 500 bp and aggregate by cluster.
#' Intended to generate input for the tileGeneM function.
#'
#' @param obj Amethyst object containing the path(s) to calculate.
#' @param type Type of methylation - e.g. "CG" or "CH" - to calculate. Typically will be "CG."
#' @param threads Optional number of threads for parallelization. Basically required for this step.
#' @param bed Bed file containing genome in 500bp windows. Options are "hg38", "mm10", or path to bed file.
#' @param index Name of chr index in the index slot.
#' This reduces memory constraints by processing one chromosome at a time.
#'
#' @return Returns a data.table with the genomic location as chr, start, and end columns, plus one column of the mean % methylation per cluster
#' @importFrom dplyr filter pull mutate select left_join
#' @importFrom stringr str_count
#' @importFrom rhdf5 h5read h5ls
#' @importFrom future plan multicore sequential
#' @importFrom furrr future_map future_pmap
#' @importFrom plyr round_any
#' @importFrom tidyr separate
#' @importFrom data.table setnames setDT tstrsplit := setorder data.table
#' @importFrom purrr reduce
#' @importFrom utils read.table
#' @export
#'
#' @examples cluster_500bp_cg_tracks <- makeClusterTracks(obj = obj, type = "CG", threads = 20, bed = "hg38")
makeClusterTracks <- function(
    obj,
    type,
    threads = 1,
    bed,
    index) {

  # get barcodes and paths from amethyst object
  if (is.null(obj@h5paths)) {
    stop("Please generate the path list for each barcode and store in the obj@h5paths slot.")
  } else {
    h5paths <- obj@h5paths
  }

  # load and clean bed file
  if (bed == "hg38") {
    bed <- data.table::data.table(read.table("/home/groups/ravnica/projects/amethyst/bedfiles/hg38_500bp_intervals.bed", sep = "\t"))
  } else if (bed == "mm10") {
    bed <- data.table::data.table(read.table("/home/groups/ravnica/projects/amethyst/bedfiles/mm10_500bp_intervals.bed", sep = "\t"))
  } else if (!(bed %in% c("hg38", "mm10"))) {
    bed <- data.table::data.table(read.table(bed, sep = "\t"))
  }
  data.table::setnames(bed, c("chr", "start", "end"))
  data.table::setorder(bed, chr, start, end)
  bed <- bed[, window := paste0(chr, "_", start, "_", end)]
  bed <- bed[stringr::str_count(window, "_") == 2] # remove alternative loci
  chr_groups <- as.list(unique(bed$chr)) # store chr groups list to loop over
  chr_groups <-chr_groups[!grepl("chrY|chrM", chr_groups)]

  # set up multithreading
  if (threads > 1) {
    future::plan(future::multicore, workers = threads)
  }

  # get number of clusters
  cluster_groups <- as.list(unique(obj@metadata[["cluster_id"]])) # make list for each unique cluster
  aggregated <- furrr::future_map(cluster_groups, .f = function(i) {

    subset <- h5paths[rownames(h5paths) %in% rownames(obj@metadata |> dplyr::filter(cluster_id == i)), , drop = FALSE] # subset paths list to those in current cluster
    paths <- as.list(subset$paths)
    barcodes <- as.list(rownames(subset))
    by_chr <- list()

    for (j in 1:length(chr_groups)) {

      # get index positions
      sites <- obj@index[[index]][[chr_groups[[j]]]] # get chr index for h5 file

      # find 500bp window values for all cells in chr + cluster
      windows <- furrr::future_pmap(.l = list(paths, barcodes), .f = function(path, barcode) {
        barcode_name <- sub("\\..*$", "", barcode)
        data <- data.table::data.table(rhdf5::h5read(path, name = paste0("CG", "/", barcode),
                                                     start = sites$start[sites$cell_id == barcode],
                                                     count = sites$count[sites$cell_id == barcode])) # read in 1 chr at a time
        data <- data[pos %% 500 == 0, pos := pos + 1] # add 1 to anything exactly divisible by window size otherwise it will be its own window
        data <- data[, window := paste0(chr, "_", plyr::round_any(pos, 500, floor), "_", plyr::round_any(pos, 500, ceiling))][, c("chr", "pos", "pct") := NULL]
        data <- data[, .(value = as.integer(sum(c != 0, na.rm = TRUE) * 100 / (sum(c != 0, na.rm = TRUE) + sum(t != 0, na.rm = TRUE)))), by = window]
        data.table::setnames(data, "value", barcode_name)
      }, .progress = TRUE)

      windows <- split(windows, ceiling(seq_along(windows)/50))
      windows <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE),
                        furrr::future_map(windows, ~ Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE), .x), .progress = FALSE))

      # windows <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE), windows)
      data.table::setDT(windows)
      windows <- windows[, .(window, rowMean = round(rowMeans(.SD, na.rm = TRUE), 1)), .SDcols = -"window"] # calculate the summary statistic of each window for the cluster
      data.table::setnames(windows, "rowMean", paste0(i))

      by_chr[[j]] <- windows
      cat("\nCompleted chr", j, " cluster ", i, "\n")
    }
    aggregated <- do.call(rbind, by_chr)
  }, .progress = TRUE)

  if (threads > 1) {
    future::plan(future::sequential)
    gc()
  }

  aggregated <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE), aggregated)
  aggregated <- dplyr::left_join(bed, aggregated, by = "window")

  # Calculate the rolling mean for each column
  aggregated[, c("window") := NULL]
  data.table::setorder(aggregated, chr, start)

  averages <- aggregated[, lapply(.SD, function(x) round(frollmean(x, n = 3, align = "center", fill = NA, na.rm = TRUE), 2)), .SDcols = -c("chr", "start", "end")]

  # Add the averages as new columns to your original data.table
  aggregated[, colnames(aggregated)[colnames(aggregated) %in% cluster_groups] := NULL]
  aggregated[, smooth_start := lapply(.SD, function(x) frollapply(x, n = 3, FUN = function(y) y[1], align = "center", fill = NA)), .SDcols = "start"]
  aggregated[, smooth_end := lapply(.SD, function(x) frollapply(x, n = 3, FUN = function(y) y[3], align = "center", fill = NA)), .SDcols = "end"]
  aggregated[, names(averages) := averages]
  aggregated[shift(chr, 1) == shift(chr, -2)]
  aggregated[, c("smooth_start", "smooth_end") := NULL]
  return(aggregated)
}

############################################################################################################################
#' @title makeFuzzyGeneWindows
#' @description Rapidly calculate % methylation over all genes by averaging 500 bp window values.
#' Slightly less accurate than the alternative indexing approach for short genes.
#'
#' @param obj Amethyst object for which to calculate % methylation over all genes
#' @param type Type of methylation to calculate; e.g. "CH".
#' "CG" is an option but may be less biologically relevant
#' @param ref Species reference to obtain bed file in 500 bp genomic windows. "hg38", "mm10", or path to bed file.
#' @param threads Enables multithreading
#'
#' @return Returns a matrix just like the output of makeWindows with rows as genes, columns as cells, and values as % methylation
#' @export
#'
#' @examples gene_ch <- makeFuzzyGeneWindows(obj = obj, type = "CH", threads = 20)
#' @importFrom data.table setDT data.table
#' @importFrom dplyr filter pull
#' @importFrom furrr future_map
#' @importFrom plyr round_any
#' @importFrom purrr reduce
#' @importFrom rhdf5 h5ls h5read
#' @importFrom future plan sequential
#' @importFrom utils read.table
makeFuzzyGeneWindows <- function(
    obj,
    type,
    ref,
    threads = 1) {

  if (threads > 1) {
    plan(multicore, workers = threads)
  }

  # download bed file
  if (ref == "hg38") {
    bed <- utils::read.table("/home/groups/ravnica/projects/amethyst/bedfiles/hg38_mainChr_500bp_geneannot.bed", header = T, sep = '\t')
  } else if (ref == "mm10") {
    bed <- utils::read.table("/home/groups/ravnica/projects/amethyst/bedfiles/mm10_mainChr_500bp_geneannot.bed", header = T, sep = '\t')
  } else if (!(ref %in% c("hg38", "mm10"))) {
    bed <- data.table::data.table(read.table(ref, sep = "\t"))
  }
  data.table::setDT(bed)
  bed <- bed[, .(window, gene)]

  # get barcodes and paths from amethyst object
  if (is.null(obj@h5paths)) {
    stop("Please generate the path list for each barcode and store in the obj@h5paths slot.")
  } else {
    barcodes <- as.list(rownames(obj@h5paths))
    paths <- as.list(obj@h5paths$paths)
  }

  # calculate % methylation over 500 bp windows for all cells
  windows <- furrr::future_pmap(.l = list(paths, barcodes), .f = function(path, bar) {
    barcode_name <- sub("\\..*$", "", bar)
    data <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", bar)))
    data[, window := paste0(chr, "_", plyr::round_any(pos, 500, floor), "_", plyr::round_any(pos, 500, ceiling))]
    result <- data[, .(value = round(sum(c != 0, na.rm = TRUE) * 100 / (sum(c != 0, na.rm = TRUE) + sum(t != 0, na.rm = TRUE)), 1)), by = window]
    setnames(result, "value", barcode_name)
    result <- dplyr::left_join(bed, result, by = c("window"))
    result <- result[, lapply(.SD, function(x) round(mean(x, na.rm = TRUE), 2)), by = gene, .SDcols = barcode_name]
    data.table::setorder(result, gene)
  }, .progress = TRUE)

  if (all(sapply(windows[-1], function(df) identical(windows[[1]]$gene, df$gene)))) { # check to make sure all the gene columns are in the same order before cbinding
    result <- do.call(cbind, c(windows[[1]][, 1:2], lapply(windows[-1], function(df) df[, 2]))) # cbind second column of all dataframes with first
  } else {
    result <- purrr::reduce(windows, dplyr::full_join, by = c("gene"))
  }

  if (threads > 1) {
    future::plan(future::sequential)
    gc()
  }

  result <- result |> tibble::column_to_rownames(var = "gene")
}

