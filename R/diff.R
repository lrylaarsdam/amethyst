# Authors: Lauren Rylaarsdam, PhD
# Note: for testDMR function, used fast fisher's exact test developed by @zellerivo; see https://github.com/al2na/methylKit/issues/96
# 2024-2025
############################################################################################################################
#' @title Find cluster markers
#' @description Find variably methylated features across clusters
#'
#' @param obj The amethyst object to analyze
#' @param genes Character vector of genes to test
#' @param nmin Minimum threshold of observations in either the member or non-member group to proceed with testing
#' @param groupBy Metadata column name of the variable to group by
#' @param method Multiple testing correction method. Options are c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param threads Optional; number of threads to use if parallelizing (recommended)
#' @param matrix Name of the pre-computed methylation values stored in the genomeMatrices slot for which to calculate variably methylated features
#'
#' @returns A data.frame containing the p.value of methylation levels over each gene in each cluster relative to all other clusters
#' @importFrom dplyr select filter mutate group_by rename
#' @importFrom stats wilcox.test p.adjust
#' @importFrom future plan multisession sequential
#' @importFrom furrr future_map
#' @export
#' @examples
#' \dontrun{
#'   gene_ch_cluster_markers <- findClusterMarkers(obj = obj, matrix = "gene_ch", genes = protein_coding, nmin = 10, groupBy = "cluster_id")
#' }

findClusterMarkers <- function(
    obj,
    matrix,
    genes,
    nmin = 10,
    groupBy = "cluster_id",
    method = "bonferroni",
    threads = 1) {

  options(scipen = 3)

  obj@metadata$cluster_id <- obj@metadata[[groupBy]]

  genematrix <- as.matrix(obj@genomeMatrices[[matrix]])
  membership <- obj@metadata |> dplyr::select("cluster_id")

  # Set up multithreading
  if (threads > 1) {
    future::plan(future::multicore, workers = threads)
  }

  results <- furrr::future_map(.x = genes, .f = function(gene) {
    gene_results <- list()  # Initialize outside the loop
    for (id in unique(membership$cluster_id)) {
      members <- rownames(membership |> dplyr::filter(cluster_id == id))
      nonmembers <- rownames(membership |> dplyr::filter(cluster_id != id))
      tryCatch(
        {
          if (length(members) < nmin || length(nonmembers) < nmin) {
            next
          } else {
            gene_results[[id]] <- data.frame(
              "p.val" = stats::wilcox.test(
                x = genematrix[gene, members],
                y = genematrix[gene, nonmembers])$p.value,
              "gene" = gene,
              "cluster_id" = id,
              mean_1 = mean(genematrix[gene, members], na.rm = TRUE),
              mean_2 = mean(genematrix[gene, nonmembers], na.rm = TRUE)
            ) |> dplyr::mutate(
              logFC = log2(mean_1 / mean_2),
              direction = ifelse(mean_1 > mean_2, "hypermethylated", "hypomethylated")
            )
          }
        },
        error = function(e) {
          cat("Warning:", gene, "and cluster", id, "not tested\n")
        }
      )
    }
    gene_results
  }, .progress = TRUE)

  if (threads > 1) {
    future::plan(future::sequential)
    gc()
  }

  results <- do.call(rbind, lapply(results, function(x) do.call(rbind, x)))
  results <- results |> dplyr::group_by(cluster_id) |> dplyr::mutate(p.adj = stats::p.adjust(p.val, method = method)) |> dplyr::select(p.val, p.adj, everything())

  return(results)

}

######################################################################
#' @title calcSmoothedWindows
#' @description Aggregate c and t observations over smoothed genomic windows of a user-defined size.
#' Optionally calculate the % methylation for each window.
#'
#' @param obj Amethyst object containing the h5 file paths
#' @param type Type of methylation - e.g. "CG" or "CH" - to calculate. Typically will be "CG"
#' @param threads Optional number of threads for parallelization (multithreading strongly recommended)
#' @param index Name of chr index in the index slot
#' @param step Width of the genomic window to calculate. Default is 500 bp
#' @param smooth Number of windows to include surrounding the target region; i.e. produces a sliding window matrix
#' Default parameter is 3, resulting in a 1500 x 500 bp sliding window matrix
#' @param genome Genome build of the organism(s) being analyzed. Options are currently "hg19", "hg38", "mm10", or "mm39".
#' @param futureType If using parallelization, should R multithread using "multicore" or "multisession"
#' @param groupBy Parameter contained in the metadata over which to aggregate observations. Default is "cluster_id"
#' @param returnSumMatrix Whether or not the function should return the matrix of summed c and t observations. Required for testDMR input
#' @param subset Optional; only calculate for portion of genome. Specify ranges like: "chr15_77563027_77870900"
#' @param save Optional; save each chromosome output intermediate. Helpful if you are working with large datasets
#' @param returnPctMatrix Whether or not the function should calculate % methylation over each genomic window. Required for heatMap input
#' @param chrList Optional; chromosome whitelist (only use if necessary)
#' @param chrSizes Optional; matched size vector of chromosomes provided in the whitelist (MUST include if chrList is provided)
#' @return if returnSumMatrix = TRUE, returns a data.table with the genomic location as chr, start, and end columns,
#' plus aggregated c and t observations for each groupBy value. If returnPctMatrix = TRUE, returns a data.table with the
#' genomic location as chr, start, and end columns, plus one column of the mean % methylation for each window per groupBy value.
#' If both are true, returns a list containing each result.
#' @importFrom dplyr filter pull mutate select left_join
#' @importFrom stringr str_count
#' @importFrom rhdf5 h5read h5ls
#' @importFrom future plan multicore sequential
#' @importFrom furrr future_map future_pmap
#' @importFrom plyr round_any .
#' @importFrom tidyr separate
#' @importFrom data.table data.table setnames setDT tstrsplit := setorder data.table frollapply frollsum frollmean shift copy
#' @importFrom purrr reduce
#' @importFrom utils read.table
#' @export
#' @examples
#' \dontrun{
#'   output <- calcSmoothedWindows(obj, groupBy = "cluster_id")
#'   output <- calcSmoothedWindows(obj, groupBy = "louvain_vmr_cluster", chrList = c("chr1", "chr2"), chrSizes = c(248956422, 242193529))
#' }
calcSmoothedWindows <- function(
    obj,
    type = "CG",
    threads = 1,
    step = 500,
    smooth = 3,
    genome = "hg38",
    subset = NULL,
    index = "chr_cg",
    futureType = "multicore",
    groupBy = "cluster_id",
    returnSumMatrix = TRUE,
    returnPctMatrix = TRUE,
    save = FALSE,
    chrList = NULL,
    chrSizes = NULL) {

  # get barcodes and paths from amethyst object
  if (is.null(obj@h5paths)) {
    stop("\nPlease generate the path list for each barcode and store in the obj@h5paths slot.")
  } else {
    h5paths <- obj@h5paths
  }

  if (!is.null(chrList) && is.null(chrSizes)) {
    stop("\nIf a chromosome whitelist is used, you also must provide a matched vector of chromosome sizes.")
  }

  # Generate genome in step size windows
  # Function to generate windows for a given chromosome
  generate_windows <- function(chromosome, size) {
    starts <- seq(0, size - 1, by = step)
    ends <- pmin(starts + step, size)
    data.frame(chr = chromosome, start = starts, end = ends)
  }


  if (is.null(chrList)) {
    if (genome %in% c("hg19", "hg38")) {
      chromosome = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                     "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                     "chr20", "chr21", "chr22", "chrX", "chrY")
      if (genome == "hg19") {
        size = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
                 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
      } else if (genome == "hg38") {
        size = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
                 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
                 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
      }

    } else if (genome %in% c("mm10", "mm39")) {
      chromosome = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                     "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                     "chrX", "chrY")
      if (genome == "mm10") {
        size = c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213,
                 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768,
                 94987271, 90702639, 61431566, 171031299, 91744698)
      } else if (genome == "mm39") {
        size = c(195154279, 181755017, 159745316, 156860686, 151758149, 149588044, 144995196, 130127694,
                 124359700, 130530862, 121973369, 120092757, 120883175, 125139656, 104073951, 98008968,
                 95294699, 90720763, 61420004, 169476592, 91455967)
      }

    } else if (!(genome %in% c("hg19", "hg38", "mm10", "mm39"))) {
      stop("Only hg19, hg38, mm10, or mm39 can be accommodated at this time.")
    }

    chromosome_sizes <- data.frame(chromosome, size)

  } else if (!is.null(chrList) && !is.null(chrSizes)) {
    chromosome_sizes <- data.frame(chromosome = chrList,
                                   size = chrSizes)
  }

  if (!is.null(subset)) {
    subset_coords <- data.table::tstrsplit(subset, "_", fixed = TRUE)
    subset_coords <- data.table(subset_chr = subset_coords[[1]], subset_start = as.integer(subset_coords[[2]]), subset_end = as.integer(subset_coords[[3]]))
    data.table::setkey(subset_coords, subset_chr, subset_start, subset_end)

    chromosome_sizes <- chromosome_sizes[chromosome_sizes$chromosome == subset_coords$subset_chr, ]
  }

  # Apply the function to each chromosome and combine results
  genomechunks <- do.call(rbind, lapply(1:nrow(chromosome_sizes), function(i) {
    generate_windows(chromosome_sizes$chromosome[i], chromosome_sizes$size[i])
  }))

  data.table::setDT(genomechunks)
  genomechunks <- genomechunks[, window := paste0(chr, "_", start, "_", end)]
  data.table::setkey(genomechunks, chr, start, end)

  if (!is.null(subset)) {
    filtered_chunks <- foverlaps(
      x = genomechunks,
      y = subset_coords,
      by.x = c("chr", "start", "end"),
      by.y = c("subset_chr", "subset_start", "subset_end"),
      nomatch = 0L,  # drops non-overlapping rows
      type = "any"   # any overlap is sufficient
    )
    genomechunks <- filtered_chunks[, .(chr, start, end, window)]
  }

  #define chr groups
  if (is.null(chrList)) {
    chr_groups <- as.list(unique(genomechunks$chr)) # store chr groups list to loop over
  } else {
    chr_groups <- as.list(chrList)
  }

  # get number of clusters
  membership <- obj@metadata |> dplyr::select(groupBy)
  colnames(membership) <- "membership"
  groups <- as.list(unique(membership$membership))


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

  by_chr <- list()
  for (chr in chr_groups) {

    # get index positions
    sites <- obj@index[[index]][[chr]] # get chr index for h5 file

    chr_group_results <- furrr::future_map(.x = groups, .f = function(gr) {
      members <- rownames(membership |> dplyr::filter(membership == gr))

      if (is.null(obj@h5paths$prefix)) {
        member_paths <- h5paths[match(members, h5paths$barcode), "path"]
        member_prefixes <- rep("", length(members))
      } else {
        h5paths$cell_ids <- paste0(obj@h5paths$prefix, obj@h5paths$barcode)
        member_paths <- h5paths[match(members, h5paths$cell_ids), "path"]
        member_prefixes <- h5paths[match(members, h5paths$cell_ids), "prefix"]
        members <- h5paths[match(members, h5paths$cell_ids), "barcode"]
      }

      # add up sum c and sum t in member cells
      member_results <- furrr::future_pmap(.l = list(member_paths, members, member_prefixes), .f = function(path, barcode, prefix) {
        tryCatch({
          barcode_name <- sub("\\..*$", "", barcode)
          data <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode, "/1"),
                                                       start = sites$start[sites$cell_id == paste0(prefix, barcode)],
                                                       count = sites$count[sites$cell_id == paste0(prefix, barcode)])) # read in 1 chr at a time
          data <- data[pos %% step == 0, pos := pos + 1] # add 1 to anything exactly divisible by window size otherwise it will be its own window
          data <- data[, window := paste0(chr, "_", plyr::round_any(pos, step, floor), "_", plyr::round_any(pos, step, ceiling))][, c("chr", "pos", "pct") := NULL]
          data <- data[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = window]
        }, error = function(e) {
          cat("Error processing data for barcode", barcode, ":", conditionMessage(e), "\n")
          return(NULL)  # Return NULL or any other value indicating failure
        })
      }, .progress = TRUE)

      # aggregate in groups so as to not exceed R row limit
      member_results <- split(member_results, ceiling(seq_along(member_results) / 100))
      member_results <- lapply(member_results, function(x) {
        x <- do.call(rbind, x)
        x <- x[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = window]
      })
      member_results <- do.call(rbind, member_results)
      member_results <- member_results[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = window]
      data.table::setnames(member_results, c("window", paste0(gr, "_c"), paste0(gr, "_t")))
    }, .progress = TRUE)

    by_chr[[chr]] <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE), chr_group_results)

    if (save) {
      saveRDS(by_chr[[chr]], paste0("tmp_", type, "_", step, "bp_", chr, ".RData"))
    }

    cat("\nCompleted ", chr,"\n")
  }

  if (threads > 1) {
    future::plan(future::sequential)
    gc()
  }

  count_matrix <- do.call(rbind, by_chr)
  count_matrix <- dplyr::left_join(genomechunks, count_matrix, by = "window")

  # Calculate the rolling mean for each column
  count_matrix <- count_matrix[, c("window") := NULL]
  data.table::setorder(count_matrix, chr, start)

  rolling_sums <- count_matrix[, lapply(.SD, function(x) as.integer(data.table::frollsum(x, n = smooth, align = "center", fill = NA, na.rm = TRUE))), .SDcols = grep("_c$|_t$", names(count_matrix))]

  # Add the rolling sums as new columns to your original data.table
  count_matrix <- count_matrix[, smooth_start := lapply(.SD, function(x) data.table::frollapply(x, n = smooth, FUN = function(y) y[1], align = "center", fill = NA)), .SDcols = "start"]
  count_matrix <- count_matrix[, smooth_end := lapply(.SD, function(x) data.table::frollapply(x, n = smooth, FUN = function(y) y[smooth], align = "center", fill = NA)), .SDcols = "end"]
  count_matrix <- count_matrix[, names(rolling_sums) := rolling_sums]
  count_matrix <- count_matrix[data.table::shift(chr, 1) == data.table::shift(chr, -2)]
  count_matrix <- count_matrix[, c("smooth_start", "smooth_end") := NULL]
  count_matrix <- count_matrix[rowSums(count_matrix[, .SD, .SDcols = -c("chr", "start", "end")]) != 0]

  # Calculate cluster track from count_matrix
  if (returnPctMatrix) {
    pctm_matrix <- data.table::copy(count_matrix)
    for (gr in groups) {
      group_c_cols <- paste0(gr, "_c")
      group_t_cols <- paste0(gr, "_t")
      pctm_matrix[, m := round(.SD[[group_c_cols]] * 100 / (.SD[[group_c_cols]] + .SD[[group_t_cols]]), 2), .SDcols = c(group_c_cols, group_t_cols)]
      data.table::setnames(pctm_matrix, "m", as.character(gr))
    }

    pctm_matrix <- pctm_matrix[, c(paste0(groups, "_c"), paste0(groups, "_t")) := NULL]
  }

  if (returnSumMatrix & returnPctMatrix) {
    return(list(sum_matrix = count_matrix, pct_matrix = pctm_matrix))
  } else if (!returnSumMatrix & returnPctMatrix) {
    return(pctm_matrix)
  } else if (returnSumMatrix & ! returnPctMatrix) {
    return(count_matrix)
  }
}

######################################################################
#' @title testDMR
#' @description calculate differentially methylated regions from a preconstructed matrix of summed observations
#'
#' @param sumMatrix Sum of c and t observations in each genomic window per group. Output of calcSmoothedWindows
#' @param eachVsAll If TRUE, each group found in the sumMatrix will be tested against all others
#' @param comparisons If eachVsAll is not desired, provide a data frame describing which tests to run.
#' The data.frame should have three columns with rows describing conditions of each test
#' "name" determines the name of the test in the output; "A" lists group members, and "B" lists group nonmembers
#' @param nminTotal Minimum number of observations across all groups to include the genome region in calculations
#' @param nminGroup Minimum number of observations across both members and nonmembers to include the genome region in calculations
#'
#' @return Returns a data.table containing the test results for each condition as an appended column
#' @export
#' @importFrom data.table copy setDT
#' @examples
#' \dontrun{
#'   dmrs <- testDMR(sumMatrix = calcSmoothedWindows(sumMatrix = smoothed_windows[["sum_matrix"]], eachVsAll = TRUE, nminTotal = 10, nminGroup = 10))
#' }
testDMR <- function(
    sumMatrix,
    eachVsAll = TRUE,
    comparisons = NULL,
    nminTotal = 3,
    nminGroup = 3) {

  if (!eachVsAll && is.null(comparisons)) {
    stop("Please either specify eachVsAll = TRUE or provide a data frame of comparisons to make.")
  }

  # filter counts table
  data.table::setDT(sumMatrix)
  counts <- data.table::copy(sumMatrix)
  counts <- counts[rowSums(counts[, .SD, .SDcols = patterns("_c$|_t$")], na.rm = TRUE) >= nminTotal]

  # fast fisher's exact test developed by @zellerivo; see https://github.com/al2na/methylKit/issues/96
  fast.fisher <- function (
    cntg_table) {
    q <- cntg_table[1, 1]
    m <- cntg_table[1, 1] + cntg_table[2, 1]
    n <- cntg_table[1, 2] + cntg_table[2, 2]
    k <- cntg_table[1, 1] + cntg_table[1, 2]
    pval_right <- phyper(q = q, m = m, n = n, k = k, lower.tail = FALSE) +
      (0.5 * dhyper(q, m, n, k))
    pval_left <- phyper(q = q - 1, m = m, n = n, k = k, lower.tail = TRUE) +
      (0.5 * dhyper(q, m, n, k))
    return(ifelse(test = pval_right > pval_left, yes = pval_left *
                    2, no = pval_right * 2))
  }

  if (is.null(comparisons)) {
    # get unique groups
    groups <- as.list(sub("_c$", "", colnames(sumMatrix)[grep("_c$", colnames(sumMatrix))]))

    for (gr in groups) {
      m_c <- paste0(gr, "_c") # m = member
      m_t <- paste0(gr, "_t")

      nm_c <- setdiff(grep("_c$", colnames(counts), value = TRUE), m_c) # nm = nonmember
      nm_t <- setdiff(grep("_t$", colnames(counts), value = TRUE), m_t)

      counts <- counts[, `:=`(
        member_c = get(paste0(gr, "_c")),
        member_t = get(paste0(gr, "_t")),
        nonmember_c = rowSums(.SD[, mget(nm_c)]),
        nonmember_t = rowSums(.SD[, mget(nm_t)])
      )]

      # don't test where the minimum observations per group is not met
      counts <- counts[member_c + member_t <= nminGroup | nonmember_c + nonmember_t <= nminGroup, c("member_c", "member_t", "nonmember_c", "nonmember_t") := .(NA, NA, NA, NA)]

      # apply fast fishers exact test
      counts <- counts[, paste0(gr, "_all_pval") := apply(.SD, 1, function(x) fast.fisher(matrix(x, nrow = 2, byrow = TRUE))), .SDcols = c("member_c", "member_t", "nonmember_c", "nonmember_t")]
      counts <- counts[, paste0(gr, "_all_logFC") := round(log2((member_c / (member_c + member_t)) / (nonmember_c / (nonmember_c + nonmember_t))), 4)]
      counts <- counts[, c("member_c", "member_t", "nonmember_c", "nonmember_t") := NULL] # this line used to be outside the loop in < v1.0.2, causing nonmember variable buildup :(
      cat(paste0("Finished group ", gr, "\n"))
    }

  } else if (!is.null(comparisons)) {
    for (i in 1:nrow(comparisons)) {
      m <- unlist(strsplit(comparisons[i, "A"], ','))
      nm <- unlist(strsplit(comparisons[i, "B"], ',', fixed = FALSE))
      name <- comparisons[i, "name"]

      m_c <- paste0(m, "_c") # m = member
      m_t <- paste0(m, "_t")

      nm_c <- paste0(nm, "_c") # n = nonmember
      nm_t <- paste0(nm, "_t")

      counts <- counts[, `:=`(
        member_c = rowSums(.SD[, mget(m_c)]),
        member_t = rowSums(.SD[, mget(m_t)]),
        nonmember_c = rowSums(.SD[, mget(nm_c)]),
        nonmember_t = rowSums(.SD[, mget(nm_t)])
      )]

      # don't test where the minimum observations per group is not met
      counts <- counts[member_c + member_t <= nminGroup | nonmember_c + nonmember_t <= nminGroup, c("member_c", "member_t", "nonmember_c", "nonmember_t") := .(NA, NA, NA, NA)]

      # apply fast fishers exact test
      counts <- counts[, paste0(name, "_pval") := apply(.SD, 1, function(x) fast.fisher(matrix(x, nrow = 2, byrow = TRUE))), .SDcols = c("member_c", "member_t", "nonmember_c", "nonmember_t")]
      counts <- counts[, paste0(name, "_logFC") := round(log2((member_c / (member_c + member_t)) / (nonmember_c / (nonmember_c + nonmember_t))), 4)]
      counts <- counts[, c("member_c", "member_t", "nonmember_c", "nonmember_t") := NULL] # this line used to be outside the loop in < v1.0.2, causing nonmember variable buildup :(
      cat(paste0("Finished testing ", name, ": ", paste0(m, collapse = ", "), " vs. ", paste0(nm, collapse = ", ")))
    }
  }
  return(counts)
}

######################################################################
#' @title filterDMR
#' @description Expand and filter the testDMR output
#'
#' @param dmrMatrix Results from testDMR
#' Should contain a data.table of summed observations per test group and genomic region with unadjusted p values and logFC for each condition.
#' @param method Method for p value adjustment. Bonferroni is the default. See ?p.adjust for all options.
#' @param filter If TRUE, removes insignificant results.
#' @param keepSums If TRUE, does not remove summed c and t observations per group from the resulting data.table.
#' @param pThreshold Maxmimum adjusted p value to allow if filter = TRUE.
#' @param logThreshold Minimum absolute value of the log2FC to allow if filter = TRUE.
#' @return Returns an expanded data table with results from testDMR.
#' @export
#' @importFrom data.table data.table copy := melt
#' @importFrom stats p.adjust
#' @examples
#' \dontrun{
#'   results <- filterDMR(dmrMatrix = dmrs, method = "BH", filter = TRUE, keepSums = FALSE, pThreshold = 0.05, logThreshold = 2)
#' }
filterDMR <- function(
    dmrMatrix,
    method = "bonferroni",
    filter = TRUE,
    keepSums = FALSE,
    pThreshold = 0.01,
    logThreshold = 1) {

  options(scipen = 3)
  ids <- sub("_pval$", "", colnames(dmrMatrix)[grep("_pval$", colnames(dmrMatrix))])
  results <- data.table::copy(dmrMatrix)
  results <- results[!rowSums(is.na(results[, .SD, .SDcols = -c("chr", "start", "end")])) == ncol(results[, .SD, .SDcols = -c("chr", "start", "end")])]
  if (keepSums) {
    results <- data.table::melt(results, id.vars = c("chr", "start", "end", paste0(ids, "_c"), paste0(ids, "_t")), measure.vars = patterns("_pval$", "_logFC$"), variable.name = "test_order", value.name = c("pval", "logFC"), na.rm = TRUE)
  } else if (!(keepSums)) {
    results <- data.table::melt(results, id.vars = c("chr", "start", "end"), measure.vars = patterns("_pval$", "_logFC$"), variable.name = "test_order", value.name = c("pval", "logFC"), na.rm = TRUE)
  }
  results <- results[, padj := stats::p.adjust(pval, method = method)]

  if (filter) {
    results <- results[padj < pThreshold & abs(logFC) > logThreshold]
  }
  results <- results[, direction := ifelse(logFC < 0, "hypo", "hyper")]

  key <- data.table::data.table(
    test_order = as.factor(1:length(ids)),  # assuming test is integer in your main table
    test = as.character(ids)  # or any corresponding label
  )

  results <- merge(results, key, by = "test_order", all.x = TRUE)
  results <- results[, c("test_order") := NULL]

  return(results)
}

############################################################################################################################
#' @title collapseDMR
#' @description Identify adjacent DMRs
#'
#' @param obj Amethyst object containing the h5 file paths
#' @param dmrMatrix Table of differentially methylated regions identified with testDMR and filterDMR
#' @param maxDist Numeric; Maximum allowable overlap between DMRs to be considered adjacent
#' @param reduce Boolean; Reduce results to unique observations (recommended)
#' @param mergeDirections Boolean; If TRUE, merges adjacent DMRs regardless of direction (hypo/hypermethylated).
#' Recommended to keep FALSE. Can be useful for annotating large loci of interest
#' @param annotate Boolean; Add column with overlapping gene names. Requires annotation file in @ref slot
#' @param minLength Minimum length of collapsed DMR window to include in the output
#'
#' @return Returns a data.table of collapsed DMRs
#' @export
#' @importFrom data.table copy setorder as.data.table setnames
#' @importFrom GenomicRanges reduce
#' @importFrom IRanges IRanges
#' @examples
#' \dontrun{
#'   results <- collapseDMR(dmrs, maxDist = 2000, minLength = 2000, reduce = T, annotate = T)
#' }
collapseDMR <- function(
    obj,
    dmrMatrix,
    maxDist = 0,
    minLength = 1000,
    reduce = TRUE,
    mergeDirections = FALSE,
    annotate = TRUE) {

  x <- data.table::copy(dmrMatrix)
  x <- x[, c("start", "end") := list(start - maxDist, end + maxDist)]
  data.table::setorder(x, test, direction, chr, start)

  if (mergeDirections) {
    collapsed <- x[, data.table::as.data.table(GenomicRanges::reduce(IRanges::IRanges(start, end))), by = .(test, chr)]
  } else {
    collapsed <- x[, data.table::as.data.table(GenomicRanges::reduce(IRanges::IRanges(start, end))), by = .(test, chr, direction)]
  }

  data.table::setnames(collapsed, old = c("start", "end"), new = c("dmr_start", "dmr_end"))

  if (mergeDirections) {
    output <- dmrMatrix[collapsed, on = .(chr = chr, test = test, start >= dmr_start, end <= dmr_end),
                        .(chr,
                          start = x.start,
                          end = x.end,
                          test,
                          pval = x.pval,
                          padj = x.padj,
                          logFC = x.logFC,
                          dmr_start = (dmr_start + maxDist),
                          dmr_end = (dmr_end - maxDist),
                          dmr_length = (width - (2 * maxDist) - 1))]
    output <- output[, c("dmr_padj", "dmr_logFC") := list(mean(padj), mean(abs(logFC))),
                     by = .(test, chr, dmr_start, dmr_end)]
    output <- output[dmr_length >= minLength]

    if (reduce) {
      output <- unique(output[, .(chr, test, dmr_start, dmr_end, dmr_length, dmr_padj, dmr_logFC), nomatch = 0 ])
    }

  } else {
    output <- dmrMatrix[collapsed, on = .(chr = chr, direction = direction, test = test, start >= dmr_start, end <= dmr_end),
                        .(chr,
                          start = x.start,
                          end = x.end,
                          test,
                          pval = x.pval,
                          padj = x.padj,
                          logFC = x.logFC,
                          direction = direction,
                          dmr_start = (dmr_start + maxDist),
                          dmr_end = (dmr_end - maxDist),
                          dmr_length = (width - (2 * maxDist) - 1))]
    output <- output[, c("dmr_padj", "dmr_logFC") := list(mean(padj), mean(logFC)),
                     by = .(test, direction, chr, dmr_start, dmr_end)]
    output <- output[dmr_length >= minLength]

    if (reduce) {
      output <- unique(output[, .(chr, test, direction, dmr_start, dmr_end, dmr_length, dmr_padj, dmr_logFC), nomatch= 0 ])
    }

  }

  if (annotate) {
    if (is.null(obj@ref)) {
      stop("\nReference slot must be filled to annotate.")
    }
    genes <- obj@ref |> dplyr::filter(type == "gene") |>
      dplyr::select(seqid, start, end, gene_name) |>
      dplyr::distinct(gene_name, .keep_all = TRUE) |>
      dplyr::rename("chr" = "seqid", "gene" = "gene_name") |>
      data.table::data.table()

    setkey(output, chr, dmr_start, dmr_end)
    setkey(genes, chr, start, end)
    overlaps <- foverlaps(output, genes, by.x = c("chr", "dmr_start", "dmr_end"), by.y = c("chr", "start", "end"), type="any")
    overlaps[, gene_names := paste(unique(gene), collapse = ", "), by = .(chr, dmr_start, dmr_end)] # Handle multiple overlaps by concatenating gene names
    overlaps_unique <- unique(overlaps[, .(chr, dmr_start, dmr_end, gene_names)]) # Select necessary columns and remove duplicates
    output <- merge(output, overlaps_unique, by = c("chr", "dmr_start", "dmr_end"), all.x = TRUE)  # Merge the gene names back into the original output data.table

  }

  return(output)
}

