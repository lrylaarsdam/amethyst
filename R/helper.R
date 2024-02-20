############################################################################################################################
#' getGeneM
#' Helper function to retrieve the relevant portion of the hdf5 file for a given gene
#'
#' @param obj The object containing the path to the hdf5 file
#' @param gene Which gene to fetch hdf5 indexes
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
#' @param cells Optional parameter to retrieve indexes for a subset of cells
#'
#' @return Returns a data frame listing the start and count for each barcode with reads covering the requested gene
#' @export
#'
#' @examples getGeneM(obj, gene = "RBFOX3", type = "CH")
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
  gene_m <- map(.x = ids,
                .f = function(x) {
                  dataset_name <- paste0(type, "/", x)
                  start <- indexes %>% filter(cell_id == x) %>% pull(start)
                  count <- indexes %>% filter(cell_id == x) %>% pull(count)
                  h5read(obj@h5path, name = dataset_name, start = start, count = count, stride = 1)
                })

  # The result will be a list. Name the list with the corresponding cell barcodes.
  names(gene_m) <- ids

  # Row bind the list into one big list
  gene_m <- do.call(rbind, gene_m)

  # Make a new column containing the barcodes and remove unnecessary characters added automatically by R
  gene_m <- gene_m %>% rownames_to_column(var = "cell_id") %>% mutate(cell_id = sub("\\..*$", "", cell_id))
  output <- gene_m
}


############################################################################################################################
#' makeWindows
#' Calculate percent methylation across a genomic window of given size or a gene body.
#' If desired, aggregate average percent methylation according to a categorical value in the metadata.
#' @param obj Object for which to calculate percent methylation windows
#' @param groupBy Optional metadata parameter to group cells by. If included, the function will calculate mean percent methylation of a given window over the group.
#' @param genes If specified, the % methylation will be calculated over the indicated gene bodies.
#' @param stepsize If specified, the % methylation will be calculated over genomic windows of the indicated size.
#' @param nmin Optional. Minimum number of cytosines recorded in the window/gene to be included in the output.
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
#'
#' @return Returns a data frame with the genomic region or gene name as column names, cell barcode or group as row names, and values as percent methylated.
#' @export
#'
#' @examples makeWindows(obj, groupBy = cluster_id, genes = c("TBR2", "TBR1"), nmin = 10, type = "CH")
#' @examples makeWindows(obj, stepsize = 100000, type = "CH")
makeWindows <- function(
  obj,
  groupBy = NULL, # optional metadata parameter to group cells by
  genes = NULL, # provide only genes, bed, or step size
  stepsize = FALSE, # provide only genes, bed, or step size
  nmin = NULL,  # nmin is the minimum number of cytosines with data to include the cell in the analysis
  type, # What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
  threads = 1) {

  if (threads > 1) {
    plan(multicore, workers = threads)
  }

  # If stepsize is specified, make fixed size genomic windows
  if (stepsize) {
    # Optional aggregation over a categorical variable in obj@metadata
    if (!is.null(groupBy)) {
      groups <- as.list(unique(obj@metadata[[groupBy]]))
      aggregated <- list()
      for (i in groups) {
        barcodes <- rownames(obj@metadata %>% dplyr::filter(.data[[groupBy]] == i)) # Retrieve barcodes for the cells belonging to the appropriate metadata
        windows <- future_map(.x = barcodes, .f = function(x) {h5read(obj@h5path, name = paste0(type, "/", x)) %>%
            mutate(window = paste0(chr, "_", round_any(pos, stepsize, floor), "_", round_any(pos, stepsize, ceiling))) %>%
            dplyr::filter(pct %in% c(0, 100)) %>% dplyr::group_by(window) %>% dplyr::summarise(pct_m = sum(c != 0)*100/(sum(c != 0) + sum(t != 0)), .groups = "keep")},
            .progress = TRUE) # read in each barcode and calculate % methylation over 100kb genomic windows for each cell
        windows <- bind_rows(windows)
        windows <- windows %>% dplyr::group_by(window) %>% dplyr::summarise(pct_m = mean(pct_m)) # summarise average methylation across the group for each window
        windows <- pivot_wider(windows, names_from = window, values_from = pct_m) # change the format so genomic windows are columns
        rownames(windows) <- paste(as.character(i)) # set the row name to be the group
        aggregated[[i]] <- windows # add windows aggregated across the group to a list
      }
      windows <- bind_rows(aggregated) # aggregate the data frame after windows for each group are calculated. Columns will be genomic windows; rows will be the group name.
    } else {
      barcodes <- h5ls(obj@h5path)
      barcodes <- as.list(unique(barcodes %>% dplyr::filter(otype == "H5I_DATASET") %>% dplyr::pull(name))) # extract unique barcodes from the h5 file
      windows <- future_map(barcodes, .f = function(x) {h5read(obj@h5path, name = paste0(type, "/", x)) %>%
          mutate(window = paste0(chr, "_", round_any(pos, stepsize, floor), "_", round_any(pos, stepsize, ceiling))) %>%
          dplyr::filter(pct %in% c(0, 100)) %>% dplyr::group_by(window) %>% dplyr::summarise(pct_m = sum(c != 0)*100/(sum(c != 0) + sum(t != 0)), .groups = "keep")},
          .progress = TRUE) # read in each barcode and calculate % methylation over 100kb genomic windows for each cell
      names(windows) <- barcodes
      windows <- bind_rows(windows, .id = "cell_id")
      windows <- pivot_wider(windows, id_cols = cell_id, names_from = window, values_from = pct_m) %>% column_to_rownames(var = "cell_id") # change format so columns are genomic windows and rows are cells
    }
  }

  # If genes are specified, calculate percent methylation over the indicated gene bodies
  if(!is.null(genes)) {  # based on protein coding genes from gtf. If genes with the same name are present, the longest is taken.
    # Get genomic locations (chr, start, and end) from the reference
    ranges <- obj@ref %>% dplyr::filter(type == "gene") %>% dplyr::mutate(length = end - start) %>% dplyr::group_by(gene_name) %>% dplyr::mutate(replicate = n()) %>%
      dplyr::group_by(gene_name) %>% arrange(desc(length)) %>% dplyr::filter(row_number()==1) %>% select("seqid", "start", "end", "gene_name", "location") %>% arrange(seqid, start) %>%
      filter(gene_name %in% genes)
    if (!is.null(groupBy)) {
      groups <- as.list(unique(obj@metadata[[groupBy]]))
      aggregated <- list()
      for (i in groups) {
        barcodes <- rownames(obj@metadata %>% dplyr::filter(.data[[groupBy]] == i))
        windows <- lapply(X = ranges$gene_name, FUN = getGeneM, cells = barcodes, obj = {{obj}}, type = {{type}})
        names(windows) <- ranges$gene_name
        windows <- dplyr::bind_rows(windows, .id = "gene")
        windows <- windows %>% dplyr::group_by(gene, cell_id) %>% dplyr::summarise(n = n(), pct_m = sum(c != 0)*100/(sum(c != 0) + sum(t != 0)), .groups = "keep")
        if (!is.null(nmin)) {
          windows <- windows %>% dplyr::filter(n >= nmin)
        }
        windows <- windows %>% dplyr::group_by(gene) %>% dplyr::summarise(pct_m = mean(pct_m))
        windows <- pivot_wider(windows, names_from = gene, values_from = pct_m)
        rownames(windows) <- paste(as.character(i))
        aggregated[[i]] <- windows
      }
      windows <- bind_rows(aggregated)
    } else {
      windows <- lapply(X = ranges$gene_name, FUN = getGeneM, obj = {{obj}}, type = {{type}})
      names(windows) <- ranges$gene_name
      windows <- future_map(.x = windows, .f = function(x) {x %>% dplyr::group_by(cell_id) %>% dplyr::summarise(pct_m = sum(c != 0)*100/(sum(c != 0) + sum(t != 0)), n = n(), .groups = "keep")},
                            .progress = TRUE)
      windows <- dplyr::bind_rows(windows, .id = "gene")
      if (!is.null(nmin)) {
        windows <- windows %>% dplyr::filter(n >= nmin)
      }
      windows <- pivot_wider(windows, id_cols = cell_id, names_from = gene, values_from = pct_m) %>% column_to_rownames(var = "cell_id")
    }
  }

  if (threads > 1) {
    future::plan(NULL)
    gc()
  }

  output <- windows
}


############################################################################################################################
#' addMetadata
#' Add new cell metadata to the metadata slot of a pre-existing object. Make sure row names are cell IDs.
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
    obj@metadata <- left_join(obj@metadata %>% rownames_to_column(var = "cell_id"), metadata %>% rownames_to_column(var = "cell_id"), by = "cell_id") %>% column_to_rownames(var = "cell_id")
    output <- obj
  }
}
