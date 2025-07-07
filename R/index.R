# Author: Lauren Rylaarsdam, PhD
# 2024-2025
############################################################################################################################
#' @title indexChr
#' @description If the whole hdf5 file had to be searched for relevant reads
#' every time gene-specific methylation information was needed, most functions
#' would take minutes to hours to run. The indexing process catalogs the
#' coordinates for each specified gene in every cell beforehand so the relevant
#' subset can quickly be retrieved for downstream processes. Typically, CG and CH
#' contexts would be indexed. This step can be skipped if using pre-computed
#' windows with Facet.
#'
#' @param type What type of methylation to retrieve; e.g. CG or CH.
#' @param obj Amethyst object containing h5 paths
#' @param chrList Optional; chromosome whitelist (only use if necessary)
#' @param threads Optional; number of threads to use if parallelizing (recommended)
#'
#' @return A list of data frames named by cell barcode. The data frames contain the
#' coordinates in the hdf5 file corresponding to each chromosome.
#' @export
#' @examples
#' \dontrun{
#'   obj@index[["chr_cg"]] <- indexChr(obj = obj, type = "CG", threads = 10)
#'   obj@index[["chr_ch"]] <- indexChr(obj = obj, type = "CG", chrList = c("chr1", "chr2"))
#' }
#' @importFrom data.table data.table rbindlist :=
#' @importFrom dplyr filter pull select
#' @importFrom furrr future_map
#' @importFrom future plan multicore sequential
#' @importFrom rhdf5 h5read h5ls
indexChr <- function(obj,
                     type,
                     chrList = NULL,
                     threads = 1) {

  # get barcodes and paths from amethyst object
  if (is.null(obj@h5paths)) {
    stop("Please generate the path list for each barcode and store in the obj@h5paths slot.")
  } else {
    barcodes <- obj@h5paths$barcode
    paths <- obj@h5paths$path
    if (!is.null(obj@h5paths$prefix)) {
      prefixes <- obj@h5paths$prefix
    } else {
      prefixes <- rep("", length(barcodes))
    }
  }

  if (threads > 1) {
    future::plan(future::multicore, workers = threads)
  }

  if (is.null(chrList)) {
    cat("Notice: If a chrList is not specified, chromosomes with underscores will be removed. This is to avoid indexing unwanted data like alternative contigs.\n")
  }

  output <- list()
  output <- furrr::future_pmap(.l = list(paths, barcodes, prefixes), .f = function(path, bar, prefix) {
    tryCatch({
      unique_id <- paste0(prefix, bar)
      h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", bar, "/1")))
      h5[, index := 1:.N]
      # for each gene, filter the reads based on corresponding chromosome number and make sure the position is in between the gene's start and end.
      # The minimum row number (index) is the start and the number of rows is the count. This portion of the hdf5 file can now be quickly read in for future gene-specific functions.

      #define chr groups
      if (is.null(chrList)) {
        chrs <- as.list(unique(h5$chr[!grepl("_|EBV|M", h5$chr)]))
      } else {
        chrs <- chrList
      }

      index <- furrr::future_map(chrs,
                                 .f = function(x) {
                                   ind <- h5[chr == x]
                                   ind <- ind[, .(cell_id = unique_id, chr = x, start = min(index), count = .N)]
                                   ind
                                 }, .progress = FALSE)
      index <- do.call(rbind, index) # bind to make one data frame
    }, error = function(e) {
      cat("Error processing data for barcode", bar, ":", conditionMessage(e), "\n")
      return(NULL)
    })
  }, .progress = TRUE) # show a progress bar

  if (threads > 1) {
    future::plan(future::sequential)
    gc()
  }
  output <- do.call(rbind, output)
  output <- split(output, output$chr)
  output
}



