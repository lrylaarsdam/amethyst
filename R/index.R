#' @title Find index
#' @description Worker function to find a single index used in the indexing function
#' @param curr_gene Gene name
#' @param ranges Ranges data.table
#' @param h5 A h5 file data.table
#' @returns An entry for the index list
findIndex <- function(curr_gene, ranges, h5) {
  rel_entries <- which(ranges$gene_name == curr_gene)
  rel_chr <- ranges$seqid[rel_entries]
  pos_lower_bound <- ranges[rel_entries, 2]
  pos_upper_bound <- ranges[rel_entries, 3]
  h5 <- h5[chr == rel_chr & pos %inrange% c(pos_lower_bound, pos_upper_bound)]
  h5 <- h5[, ":="(gene = curr_gene, start = min(index), count = .N), by = chr]
  ind <- h5[, .(gene, chr, start, count)]
  ind
}

#' @title Find ranges
#' @description Worker function to find the ranges used in the indexing function
#' @param gtf Genome annotation file with chromosome, start, and end position
#' information for genes of interest. See the "makeRef" function.
#' @param promoter Bool indicating if the requested type is for 'promoters';
#' used internally as a flag.
#' @returns A data.table of ranges
findRanges <- function(gtf, promoter, subset = NULL) {
  subset_passed <- !is.null(subset)
  gtf_dt <- data.table::data.table(gtf)
  if (subset_passed) {
    gtf_dt <- gtf_dt[gene_name %in% subset]
  }
  gtf_dt <- gtf_dt[type == "gene" & gene_type == "protein_coding", length := end - start]
  gtf_dt <- gtf_dt[order(-length)]
  gtf_dt <- unique(gtf_dt, by = "gene_name")

  # determine chromosome number, start position, and stop position of genes to
  # index
  if (!promoter) {
    gtf_dt <- gtf_dt[, .(seqid, start, end, gene_name, location)]
    ranges <- setorderv(
      gtf_dt,
      c("seqid", "start"),
      c(1, 1)
    )
  }
  # if indexing promoter positions, determine first base pair position of the
  # gene (sense-aware) and get coordinates of the surrounding 3kb.
  if (promoter) {
    gtf_relevant <- gtf_dt[, tss := start][strand == "+", tss := end]
    gtf_relevant <- gtf_dt[, ":="(promoter_start = tss - 1500, promoter_end = tss + 1500)]
    gtf_relevant <- unique(gtf_dt, by = c("seqid", "start", "end"))
    gtf_relevant <- gtf_dt[seqid != "chrM", .(seqid, promoter_start, promoter_end, gene_name, location)]
    ranges <- setorderv(
      gtf_dt,
      c("seqid", "promoter_start"),
      c(1, 1)
    )
  }
  ranges
}


#' @title Index genes
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
#'
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
indexGenes <- function(hdf5,
                       gtf,
                       type,
                       threads = 10,
                       subset = NULL,
                       method = "foreach",
                       verbose = FALSE) {
  if (verbose) {
    cat(paste("\nIndexing genes of type ", type))
  }
  # Check that an appropriate type has been passed
  wrong_type_passed <- !(type == "CG" | type == "CH" | type == "CAC" | type == "promoters")
  if (wrong_type_passed) {
    stop("`type` must be one of 'CG', 'CH', 'CAC' or 'promoters'.")
  }
  # Create a flag indicating if subsets passed and if the data type is 'promoters'
  subset_passed <- !is.null(subset)
  promoter_type <- type == "promoters"

  # Find the genomic ranges?
  if (verbose) {
    cat("\nFinding relevant genomic regions.")
  }
  ranges <- findRanges(gtf, promoter = promoter_type, subset = subset)

  # File leading string used in the H5I_DATASET file names for the barcodes
  file_lead <- paste0(type, "/")
  if (promoter_type) {
    file_lead <- "CG/"
  }

  # get a list of unique barcodes (cells) in the h5 file
  barcodes <- data.table::data.table(rhdf5::h5ls(hdf5))

  barcodes <- unlist(
    na.omit(
      unique(
        barcodes[
          otype == "H5I_DATASET",
          file_name := paste0(file_lead, name)
        ][, .(file_name)],
        by = "file_name"
      )
    )
  )

  # # optional set up multithreading
  # if (threads > 1) {
  #   plan(multicore, workers = threads)
  # }

  n_barcodes <- length(barcodes)
  barcodes_iter <- seq(1, n_barcodes)

  coords <- as.list(ranges$gene_name)
  n_coords <- length(coords)
  coords_iter <- seq(1, n_coords)

  if (verbose) {
    cat("\nCreate index.")
  }

  if (method == "foreach") {
    # Iterate over the file names reading in the file

    output <- foreach(bar = barcodes) %dofuture% {
      h5 <- data.table::data.table(rhdf5::h5read(hdf5, name = bar))
      h5[, index := 1:.N] # this is just a row number identifier. The hdf5 base calls are already in ascending chromosome and position order.
      index <- foreach(coord = coords) %dofuture% {
        .x <- findIndex(coord, ranges, h5)
      }
      names(index) <- as.list(ranges$gene_name) # rename the list according to gene names
      index <- data.table::rbindlist(index) # bind to make one data frame
      index$cell_id <- paste(bar) # append the corresponding cell barcode
      index
    }
  }

  if (method == "for") {
    # generate empty list to write to
    output <- vector("list", n_barcodes)

    # Iterate over the file names reading in the file
    for (ii in barcodes_iter) {
      .bar <- barcodes[[ii]]
      h5 <- data.table::data.table(rhdf5::h5read(hdf5, name = .bar))
      h5[, index := 1:.N] # this is just a row number identifier. The hdf5 base calls are already in ascending chromosome and position order.

      index <- vector("list", n_coords)
      for (jj in coords_iter) {
        index[[jj]] <- findIndex(coords[[jj]], ranges, h5)
      }
      names(index) <- as.list(ranges$gene_name) # rename the list according to gene names
      index <- data.table::rbindlist(index) # bind to make one data frame
      index$cell_id <- paste(.bar) # append the corresponding cell barcode
      output[[ii]] <- index
    }
    h5closeAll()
  }
  if (method == "lauren") {
    # read in each barcode
    output <- future_map(barcodes, function(x) {
      h5 <- h5read(hdf5, name = x)
      h5 <- h5 %>%
        dplyr::mutate(index = 1:nrow(.)) # this is just a row number identifier. The hdf5 base calls are already in ascending chromosome and position order.

      # for each gene, filter the reads based on corresponding chromosome number and make sure the position is in between the gene's start and end.
      # The minimum row number (index) is the start and the number of rows is the count. This portion of the hdf5 file can now be quickly read in for future gene-specific functions.
      coords <- as.list(ranges$gene_name)
      index <- future_map(coords,
        .f = function(x) {
          rel_entries <- which(ranges$gene_name == x)
          rel_chr <- ranges$seqid[rel_entries]
          pos_lower_bound <- ranges[rel_entries, 2]
          pos_upper_bound <- ranges[rel_entries, 3]
          h5_loc <- h5[chr == rel_chr & pos %inrange% c(pos_lower_bound, pos_upper_bound)]
          h5_loc <- h5_loc[, ":="(gene = x, start = min(index), count = .N), by = chr]
          ind <- h5_loc[, .(gene, chr, start, count)]
          ind
        },
        .progress = FALSE
      )

      names(index) <- as.list(ranges$gene_name) # rename the list according to gene names
      index <- do.call(rbind, index) # bind to make one data frame
      index$cell_id <- paste(x) # append the corresponding cell barcode
      return(index)
    }, .progress = FALSE) # show a progress bar
  }

  if (verbose) {
    cat("\nCompleted.")
  }

  # # return to one thread if multithreading
  # if (threads > 1) {
  #   future::plan(NULL)
  #   gc()
  # }

  output <- data.frame(data.table::rbindlist(output))

  # restructure the output to split by gene. This way the indexes for a given gene can easily be pulled for all cells for future functions.
  # output <- bind_rows(output)
  output <- split(output, output$gene)
  output
}

#' @title indexAll
#' @description Index gene body CG methylation, gene body CH methylation, and mCG status of promoters.
#' Promoters are defined by the first position of the gene body (strand-aware) +/- 1500 bp.
#' @param hdf5 Path to the hdf5 file containing base-level read information organized by methylation type and barcode.
#' @param subset List of genes to index. Optional but recommended.See the "fetchMarkers" function.
#' @param threads Optional threads parameter to use if parallelizing. Highly recommended as even small datasets will take a long time to index.
#' @param save Optional path to include for saving the intermediate index output as an RData file.
#'
#' @return Returns a list of coordinates in the hdf5 file corresponding to each gene. Organized by methylation type, gene, and cell.
#' @export
#'
#' @examples index <- indexAll(hdf5 = "~/Downloads/test.h5", subset = c("SATB2", "TBR2", "GAD1"), threads = 20, save = "~/Downloads/test_index.RData")
indexAll <- function(
    hdf5, # path to hdf5
    gtf,
    subset = markerGenes,
    threads = 10,
    save = NULL,
    method = "foreach",
    h5_approach = "read",
    choice = "all") {
  # Create an empty list to write output
  index <- list()

  do_ch <- choice %in% c("all", "CH")
  do_cg <- choice %in% c("all", "CG")
  do_promoters <- choice %in% c("all", "promoters")

  .ch <- .cg <- .p <- NULL

  # Doing one at a time since triple-nesting future_map resulted in memory issues
  if (do_ch) {
    cat("\nBeginning indexing gene body mCH status at ", format(Sys.time(), "%H:%M:%S"), ".")
    .ch <- indexGenes({{ hdf5 }}, gtf,
      type = "CH",
      subset = {{ subset }},
      threads = {{ threads }},
      method = method,
      h5_approach = h5_approach
    )
    cat(paste("\nGene body mCH is done indexing at ", format(Sys.time(), "%H:%M:%S"), ".\n"))
    print(str(.ch))
  }
  if (do_cg) {
    cat("\nBeginning indexing gene body mCG status at ", format(Sys.time(), "%H:%M:%S"), ".")
    .cg <- indexGenes({{ hdf5 }}, gtf,
      type = "CG",
      subset = {{ subset }},
      threads = {{ threads }},
      method = method,
      h5_approach = h5_approach
    )
    cat(paste("\nGene body mCG is done indexing at ", format(Sys.time(), "%H:%M:%S"), ".\n"))
    print(str(.cg))
  }
  if (do_promoters) {
    cat("\nBeginning indexing promoter mCG status at ", format(Sys.time(), "%H:%M:%S"), ".")
    .p <- indexGenes({{ hdf5 }}, gtf,
      type = "promoters",
      subset = {{ subset }},
      threads = {{ threads }},
      method = method,
      h5_approach = h5_approach
    )
    cat(paste("\nPromoter mCG status is done indexing at ", format(Sys.time(), "%H:%M:%S"), ".\n"))
    print(str(.p))
  }
  index[["CH"]] <- .ch
  index[["CG"]] <- .cg
  index[["promoters"]] <- .p

  # Write the resulting index to drive as an RData file if a path is provided
  if (!is.null(save)) {
    saveRDS(index, save)
  }
  output <- index
}
