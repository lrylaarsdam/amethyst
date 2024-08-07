############################################################################################################################
#' @title findRanges
#'
#' @description Worker function to find the ranges used in the indexing function
#' @param gtf Genome annotation file with chromosome, start, and end position
#' information for genes of interest. See the "makeRef" function.
#' @param promoter Bool indicating if the requested type is for 'promoters';
#' used internally as a flag.
#' @returns A data.table of ranges
#' @importFrom data.table data.table :=
findRanges <- function(gtf, promoter, subset = NULL) {
  subset_passed <- !is.null(subset)
  gtf_dt <- data.table::data.table(gtf)
  if (subset_passed) {
    gtf_dt <- gtf_dt[gene_name %in% subset]
  }
  gtf_dt <- gtf_dt[(type == "gene" & gene_type == "protein_coding")][, length := end - start]
  gtf_dt <- gtf_dt[order(-length)]
  gtf_dt <- unique(gtf_dt, by = "gene_name")

  # determine chromosome number, start position, and stop position of genes to
  # index
  if (!promoter) {
    gtf_dt <- gtf_dt[, .(seqid, start, end, gene_name, location)]
    ranges <- data.table::setorderv(
      gtf_dt,
      c("seqid", "start"),
      c(1, 1)
    )
  }
  # if indexing promoter positions, determine first base pair position of the
  # gene (sense-aware) and get coordinates of the surrounding 3kb.
  if (promoter) {
    gtf_relevant <- gtf_dt[, tss := start][strand == "-", tss := end]
    gtf_relevant <- gtf_dt[, ":="(promoter_start = tss - 1500, promoter_end = tss + 1500)]
    gtf_relevant <- unique(gtf_dt, by = c("seqid", "start", "end", "gene_name"))
    gtf_relevant <- gtf_dt[seqid != "chrM", .(seqid, promoter_start, promoter_end, gene_name, location)]
    ranges <- data.table::setorderv(
      gtf_relevant,
      c("seqid", "promoter_start"),
      c(1, 1)
    )
  }
  ranges
}

############################################################################################################################
#' @title Index genes
#' @description If the whole hdf5 file had to be searched for relevant reads
#' every time gene-specific methylation information was needed, most functions
#' would take minutes to hours to run. The indexing process catalogs the
#' coordinates for each specified gene in every cell beforehand so the relevant
#' subset can quickly be retrieved for downstream processes. mCG, mCH, and the
#' mCG status of the promoter must be indexed separately. Providing a gene
#' subset is recommended due to memory and time constraints.
#'
#' @param obj Amethyst object containing paths to h5 files to be indexed
#' @param gtf Genome annotation file with chromosome, start, and end position
#' information for genes of interest. See the "makeRef" function.
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene
#' body mCG, or promoter mCG. Must be passed as one of "CG", "CH", or "promoter".
#' The latter will do CG over TSS +/- 1500 bp
#' @param threads Optional threads parameter to use if parallelizing. Highly
#' recommended as even small datasets will take a long time to index.
#' Defaults to 10.
#' @param subset List of genes to index. Optional but recommended.See the
#' `fetchMarkers` function.
#' @return A list of data frames named by gene. The data frames contain the
#' coordinates in the hdf5 file corresponding to the gene's genomic location for
#' each cell.
#' @export
#'
#' @examples index[["CH"]] <- indexGenes(
#'   hdf5 = "~/Downloads/test.h5",
#'   gtf = ref,
#'   type = "CH",
#'   threads = 10,
#'   subset = c("SATB2", "TBR1", "FOXG1")
#' )
#' @importFrom data.table data.table rbindlist := '%inrange%'
#' @importFrom dplyr filter pull
#' @importFrom furrr future_map
#' @importFrom future plan multicore
#' @importFrom rhdf5 h5read
indexGenes <- function(obj,
                       type,
                       threads = 1,
                       subset = NULL) {

  # Create a flag indicating if subsets passed and if the data type is 'promoters'
  subset_passed <- !is.null(subset)
  if (!(subset_passed)) {
    cat(paste("\nNo gene subset has been chosen. Consider only indexing relevant genes.
              \nAlternatively, consider using makeFuzzyGeneWindows to calculate genome-wide % methylation much more quickly.
              \n"))
  }
  if (is.null(obj@ref)) {
    stop("Please generate a reference using the makeRef function and add to the obj@ref slot.")
  }

  promoter_type <- type == "promoters"
  ranges <- findRanges(obj@ref, promoter = promoter_type, subset = subset)

  # File leading string used in the H5I_DATASET file names for the barcodes
  file_lead <- paste0(type, "/")
  if (promoter_type) {
    file_lead <- "CG/"
  }

  # get barcodes and paths from amethyst object
  if (is.null(obj@h5paths)) {
    stop("Please generate the path list for each barcode and store in the obj@h5paths slot.")
  } else {
    barcodes <- as.list(rownames(obj@h5paths))
    paths <- as.list(obj@h5paths$paths)
  }

  coords <- as.list(ranges$gene_name)

  if (threads > 1) {
    future::plan(future::multicore, workers = threads)
  }

  # read in each barcode
  output <- furrr::future_pmap(.l = list(paths, barcodes), .f = function(path, bar) {
    tryCatch({
      h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(file_lead, bar)))
      h5[, index := 1:.N]
      # for each gene, filter the reads based on corresponding chromosome number and make sure the position is in between the gene's start and end.
      # The minimum row number (index) is the start and the number of rows is the count. This portion of the hdf5 file can now be quickly read in for future gene-specific functions.
      coords <- as.list(ranges$gene_name)
      index <- furrr::future_map(.x = coords, .f = function(x) {
        rel_entries <- which(ranges$gene_name == x)
        rel_chr <- ranges$seqid[rel_entries]
        pos_lower_bound <- ranges[rel_entries, 2]
        pos_upper_bound <- ranges[rel_entries, 3]
        h5_loc <- h5[chr == rel_chr & pos %inrange% c(pos_lower_bound, pos_upper_bound)]
        ind <- h5_loc[, .(gene = x, chr = rel_chr, start = min(index), count = .N)]
        ind
      }, .progress = FALSE)
      names(index) <- as.list(ranges$gene_name) # rename the list according to gene names
      index <- do.call(rbind, index) # bind to make one data frame
      index$cell_id <- paste(bar) # append the corresponding cell barcode
      index <- index[index$count != 0, ]
    }, error = function(e) {
      cat("Error processing data for ", x, ":", conditionMessage(e), "\n")
      return(NULL)
    })
  }, .progress = TRUE) # show a progress bar

  if (threads > 1) {
    future::plan(NULL)
    gc()
  }
  output <- data.frame(data.table::rbindlist(output))
  output <- split(output, output$gene)
  output
}

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
      h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", bar)))
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
