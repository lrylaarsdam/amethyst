############################################################################################################################
#' indexGenes
#' If the whole hdf5 file had to be searched for relevant reads every time gene-specific methylation information was needed,
#' most functions would take minutes to hours to run. The indexing process catalogs the coordinates for each specified gene in
#' every cell beforehand so the relevant subset can quickly be retrieved for downstream processes. mCG, mCH, and the mCG status
#' of the promoter must be indexed separately. Providing a gene subset is recommended due to memory and time constraints.
#'
#' @param hdf5 Path to the hdf5 file containing base-level read information organized by methylation type and barcode.
#' @param gtf Genome annotation file with chromosome, start, and end position information for genes of interest. See the "makeRef" function.
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
#' @param threads Optional threads parameter to use if parallelizing. Highly recommended as even small datasets will take a long time to index.
#' @param subset List of genes to index. Optional but recommended.See the "fetchMarkers" function.
#'
#' @return Returns a list of dataframes named by gene. The dataframes contain the coordinates in the hdf5 file corresponding to the gene's genomic
#' location for each cell.
#' @export
#'
#' @examples index[["CH"]] <- indexGenes(hdf5 = "~/Downloads/test.h5", gtf = ref, type = "CH", threads = 10, subset = c("SATB2", "TBR1", "FOXG1"))
indexGenes <- function(hdf5, # path to hdf5
                       gtf = ref, # annotated genome file containing the chromosome number, start, and end positions for the genes to be indexed
                       type, # "CG", "CH", or "promoter" which will do CG over TSS +/- 1500 bp
                       threads = 10,
                       subset = markerGenes) {

  # determine chromosome number, start position, and stop position of genes to index
  if (type == "CG" | type == "CH" | type == "CAC") {
    ranges <- gtf %>% dplyr::filter(type == "gene" & gene_type == "protein_coding") %>% dplyr::mutate(length = end - start) %>% dplyr::group_by(gene_name) %>% dplyr::mutate(replicate = n()) %>%
      dplyr::group_by(gene_name) %>% arrange(desc(length)) %>% dplyr::filter(row_number()==1) %>% select("seqid", "start", "end", "gene_name", "location") %>% arrange(seqid, start)
    # if indexing promoter positions, determine first base pair position of the gene (sense-aware) and get coordinates of the surrounding 3kb.
  } else if (type == "promoters") {
    ranges <- gtf %>% filter(type == "gene" & gene_type == "protein_coding") %>% mutate(length = end - start) %>% group_by(gene_name) %>% dplyr::mutate(replicate = n()) %>%
      group_by(gene_name) %>% arrange(desc(length)) %>% dplyr::filter(row_number()==1) %>% mutate(tss = ifelse(strand == "+", start, end)) %>% mutate(promoter_start = tss - 1500, promoter_end = tss + 1500) %>%
      distinct(seqid, start, end, .keep_all = T) %>% filter(seqid != "chrM") %>% select("seqid", "promoter_start", "promoter_end", "gene_name", "location") %>%
      arrange(seqid, promoter_start)
  }
  # if subset is called, filter the coordinate list
  if(!is.null(subset)) {
    ranges <- ranges %>% dplyr::filter(gene_name %in% subset)
  }

  # get a list of unique barcodes (cells) in the h5 file
  barcodes <- h5ls(hdf5)
  barcodes <- as.list(unique(barcodes %>% dplyr::filter(otype == "H5I_DATASET") %>% dplyr::pull(name)))

  # generate empty list to write to
  output <- list()

  # optional set up multithreading
  if (threads > 1) {
    plan(multicore, workers = threads)
  }

  # read in each barcode
  output <- future_map(barcodes, function(x) {
    if (type == "CG" | type == "CH" | type == "CAC") {
      h5 <- h5read(hdf5, name = paste0(type, "/", x))
    } else if (type == "promoters") {
      h5 <- h5read(hdf5, name = paste0("CG", "/", x))
    }
    h5 <- h5 %>% dplyr::mutate(index = 1:nrow(.)) # this is just a row number identifier. The hdf5 base calls are already in ascending chromosome and position order.

    # for each gene, filter the reads based on corresponding chromosome number and make sure the position is in between the gene's start and end.
    # The minimum row number (index) is the start and the number of rows is the count. This portion of the hdf5 file can now be quickly read in for future gene-specific functions.
    coords <- as.list(ranges$gene_name)
    index <- future_map(coords,
                        .f = function(x) {
                          h5 %>% dplyr::filter(chr == paste(ranges$seqid[ranges$gene_name == x]) & pos >= as.numeric(ranges[ranges$gene_name == x, 2]) & pos <= as.numeric(ranges[ranges$gene_name == x, 3])) %>%
                            group_by(chr) %>% dplyr::summarise(gene = paste(x), start = min(index), count = nrow(.)) %>% dplyr::select(gene, chr, start, count)
                        }, .progress = TRUE)
    names(index) <- as.list(ranges$gene_name) # rename the list according to gene names
    index <- do.call(rbind, index) # bind to make one data frame
    index$cell_id <- paste(x) # append the corresponding cell barcode
    return(index)
  }, .progress = TRUE) # show a progress bar

  # return to one thread if multithreading
  if (threads > 1) {
    future::plan(NULL)
    gc()
  }
  # restructure the output to split by gene. This way the indexes for a given gene can easily be pulled for all cells for future functions.
  output <- bind_rows(output)
  output <- split(output, output$gene)
}

############################################################################################################################
#' indexAll
#' Index gene body CG methylation, gene body CH methylation, and mCG status of promoters.
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
  hdf5,  # path to hdf5
  subset = markerGenes,
  threads = 10,
  save = NULL) {

  # Create an empty list to write output
  index <- list()
  # Doing one at a time since triple-nesting future_map resulted in memory issues
  index[["CH"]] <- indexGenes({{hdf5}}, type = "CH", subset = {{subset}}, threads = {{threads}})
  print(paste("Gene body mCH is done indexing at ", format(Sys.time(), "%H:%M:%S"), ". Beginning indexing gene body mCG."))

  index[["CG"]] <- indexGenes({{hdf5}}, type = "CG", subset = {{subset}}, threads = {{threads}})
  print(paste("Gene body mCG is done indexing at ", format(Sys.time(), "%H:%M:%S"), ". Beginning indexing promoter mCG status."))

  index[["promoters"]] <- indexGenes({{hdf5}}, type = "promoters", subset = {{subset}}, threads = {{threads}})
  print(paste("Indexing complete at ", format(Sys.time(), "%H:%M:%S"), "! Finishing..."))

  # Write the resulting index to drive as an RData file if a path is provided
  if (!is.null(save)) {
    saveRDS(index, save)
  }

  output <- index
}

############################################################################################################################
# addIndex
#' Get indexes for a new gene and append it to the pre-existing index slot
#' @param obj Amethyst object to append an index
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG. Can be a list.
#' @param threads Optional threads parameter to use if parallelizing. If running locally, use with caution.
#' @param genes List of genes to index.
#'
#' @return
#' @export
#'
#' @examples
addIndex <- function(obj, # path to hdf5
                     type, # "CG", "CH", or "promoters". Can be a list as well: c("CG", "CH", "promoters")
                     threads = 1,
                     genes) {

  # make a list of the requested methylation types to index
  index <- as.list(type)

  # generate indexes for all requested genes and methylation types
  index <- future_map(index,
                      .f = function(x) {indexGenes({{hdf5}}, type = x, subset = {{genes}}, threads = {{threads}})})
  names(index) <- as.list(type)

  # for each of the three possible types of methylation, if a new gene(s) has been indexed, append it to the pre-existing index.
  if (exists(index[["CH"]])) {
    obj@index[["CH"]] <- c(obj@index[["CH"]], index[["CH"]])
  }
  if (exists(index[["CG"]])) {
    obj@index[["CG"]] <- c(obj@index[["CG"]], index[["CG"]])
  }
  if (exists(index[["promoters"]])) {
    obj@index[["promoters"]] <- c(obj@index[["promoters"]], index[["promoters"]])
  }

  # redefine the object as the output.
  output <- obj
}
