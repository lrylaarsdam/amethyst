############################################################################################################################
#' @title indexChr
#' @description If the whole hdf5 file had to be searched for relevant reads
#' every time gene-specific methylation information was needed, most functions
#' would take minutes to hours to run. The indexing process catalogs the
#' coordinates for each specified gene in every cell beforehand so the relevant
#' subset can quickly be retrieved for downstream processes. mCG, mCH, and the
#' mCG status of the promoter must be indexed separately. Providing a gene
#' subset is recommended due to memory and time constraints.
#'
#' @param hdf5 Path to the hdf5 file containing base-level read information
#' organized by methylation type and barcode.
#' @param type What type of methylation to retrieve; i.e. mCH or mCG.
#' @param threads Optional threads parameter to use if parallelizing.
#' @return A list of data frames named by cell barcode. The data frames contain the
#' coordinates in the hdf5 file corresponding to each chromosome.
#' @export
#'
#' @examples
#' @importFrom data.table data.table rbindlist :=
#' @importFrom dplyr filter pull select
#' @importFrom furrr future_map
#' @importFrom future plan multicore sequential
#' @importFrom rhdf5 h5read h5ls
indexChr <- function(obj,
                     type,
                     threads = 1) {



  # get barcodes and paths from amethyst object
  if (is.null(obj@h5paths)) {
    stop("Please generate the path list for each barcode and store in the obj@h5paths slot.")
  } else {
    barcodes <- as.list(rownames(obj@h5paths))
    paths <- as.list(obj@h5paths$paths)
  }

  if (threads > 1) {
    future::plan(future::multicore, workers = threads)
  }

  output <- list()
  output <- furrr::future_pmap(.l = list(paths, barcodes), .f = function(path, bar) {
    tryCatch({
      h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", bar, "/1")))
      h5[, index := 1:.N]
      # for each gene, filter the reads based on corresponding chromosome number and make sure the position is in between the gene's start and end.
      # The minimum row number (index) is the start and the number of rows is the count. This portion of the hdf5 file can now be quickly read in for future gene-specific functions.
      chrs <- as.list(unique(h5$chr[!grepl("_|EBV|M", h5$chr)]))
      index <- furrr::future_map(chrs,
                                 .f = function(x) {
                                   ind <- h5[chr == x]
                                   ind <- ind[, .(cell_id = bar, chr = x, start = min(index), count = .N)]
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
