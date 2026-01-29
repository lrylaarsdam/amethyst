#####################################################################################################
#                                                                                                   #
# Amethyst metacells: utilities to create and work with metacell-level Amethyst objects.            #
# Includes:                                                                                         #
#  • makeMetacells(): build metacells from a reduction                                              #
#  • metacellMeanReduction(): recompute a reduction on metacells from original cell embeddings      #
#  • makeMetacellWindows(): R-side metacell tiling / smoothing for modest datasets                  #
#  • calcMetacellSmoothedWindows(): streaming smoothed windows per metacell                         #
#  • writeMetacellJobJSON(): JSON job spec for Python metacell window H5 generator                  #
#  • H5 helpers + populateMetacellTrackFromH5(): lazy loading of metacell windows                   #
#  • metacellDims(): visualize metacell–cell relationships in embedding space                       #
#  • metaMap(): gene/region-centric metacell heatmaps with gene models                              #
#                                                                                                   #
#####################################################################################################


# =========================
# Internal utilities
# =========================

.safe_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(FALSE)
  TRUE
}

.get_embedding <- function(obj, reduction, dims) {
  if (!methods::is(obj, class(obj))) stop("obj must be an S4 object")
  if (!"reductions" %in% slotNames(obj)) stop("obj@reductions not found")
  reds <- methods::slot(obj, "reductions")
  if (is.null(reds[[reduction]])) stop(sprintf("reduction '%s' not found in obj@reductions", reduction))
  emb <- reds[[reduction]]
  if (!is.matrix(emb)) emb <- as.matrix(emb)
  if (is.null(rownames(emb))) stop("reduction matrix must have rownames as cell IDs")
  dims <- dims[dims >= 1 & dims <= ncol(emb)]
  if (length(dims) == 0) stop("no valid dims selected")
  emb[, dims, drop = FALSE]
}

.get_metadata <- function(obj) {
  if (!"metadata" %in% slotNames(obj)) stop("obj@metadata not found")
  md <- methods::slot(obj, "metadata")
  if (!is.data.frame(md)) md <- as.data.frame(md)
  if (is.null(rownames(md))) stop("obj@metadata must have rownames as cell IDs")
  md
}

.numeric_cols <- function(df) {
  vapply(df, is.numeric, logical(1))
}

# Farthest-Point Sampling (FPS) in Euclidean space for seed selection
.fps_seeds <- function(X, n_seeds, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  if (n_seeds >= n) return(seq_len(n))
  # start with point closest to mean
  mu <- colMeans(X)
  d0 <- rowSums((X - matrix(mu, n, ncol(X), byrow = TRUE))^2)
  first <- which.min(d0)
  seeds <- integer(n_seeds)
  seeds[1] <- first
  # maintain distance to nearest seed
  dmin <- rowSums((X - matrix(X[first, ], n, ncol(X), byrow = TRUE))^2)
  for (i in 2:n_seeds) {
    next_i <- which.max(dmin)
    seeds[i] <- next_i
    dnew <- rowSums((X - matrix(X[next_i, ], n, ncol(X), byrow = TRUE))^2)
    dmin <- pmin(dmin, dnew)
  }
  seeds
}

# kNN indices using RANN if available, else brute-force
.knn_indices <- function(X, k) {
  n <- nrow(X)
  k <- max(min(k, n), 1)
  
  if (.safe_require("RANN")) {
    res <- RANN::nn2(X, query = X, k = k)
    idx <- res$nn.idx %||% res$nn_idx  # support either name
    if (is.null(idx)) stop("RANN::nn2 did not return nn.idx")
    idx <- as.matrix(idx)
    if (!is.integer(idx)) idx <- matrix(as.integer(idx), nrow = nrow(idx))
    return(idx)
  }
  
  # fallback: brute-force distances (O(n^2))
  D <- as.matrix(stats::dist(X, method = "euclidean", diag = TRUE, upper = TRUE))
  ord <- apply(D, 1, function(d) order(d, decreasing = FALSE)[1:k])
  idx <- t(ord)
  if (!is.integer(idx)) idx <- matrix(as.integer(idx), nrow = nrow(idx))
  idx
}

# Compute eps from kNN distances targeting ~cellN neighbors
.estimate_eps <- function(X, cellN) {
  k <- max(ceiling(cellN), 4)
  if (.safe_require("dbscan")) {
    kd <- dbscan::kNNdist(X, k = k)
    # heuristic: median k-dist with a small inflation
    return(stats::median(kd) * 1.05)
  }
  # fallback: use average distance to k-th neighbor via brute-force
  idx <- .knn_indices(X, k)
  kth <- vapply(seq_len(nrow(idx)), function(i) {
    j <- idx[i, k]
    sum((X[i, ] - X[j, ])^2)^0.5
  }, numeric(1))
  stats::median(kth) * 1.05
}

# Aggregate numeric metadata and produce originating_cellIDs vector
.aggregate_metadata <- function(md, groups, metacell_ids) {
  # Drop any empty groups defensively
  groups <- Filter(function(ix) length(ix) > 0, groups)
  if (length(groups) == 0) stop("No non-empty metacell groups were formed")
  
  num_mask <- .numeric_cols(md)
  cell_ids <- rownames(md)
  
  out_num <- if (any(num_mask)) {
    do.call(rbind, lapply(groups, function(ix) {
      xi <- md[ix, num_mask, drop = FALSE]
      if (nrow(xi) == 0) return(rep(NA_real_, sum(num_mask)))
      suppressWarnings(colMeans(xi, na.rm = TRUE))
    }))
  } else {
    NULL
  }
  
  origin <- vapply(groups, function(ix) {
    if (length(ix) == 0) return(NA_character_)
    paste(cell_ids[ix], collapse = ",")
  }, character(1))
  
  out_md <- data.frame(originating_cellIDs = origin, stringsAsFactors = FALSE, check.names = FALSE)
  if (!is.null(out_num)) {
    out_md <- cbind(out_md, as.data.frame(out_num, check.names = FALSE))
  }
  rownames(out_md) <- metacell_ids[seq_len(nrow(out_md))]
  
  return(out_md)
}

# Build metacell reductions as group means
.aggregate_reduction <- function(X, groups, metacell_ids) {
  M <- do.call(rbind, lapply(groups, function(ix) colMeans(X[ix, , drop = FALSE])))
  rownames(M) <- metacell_ids
  M
}

.groups_knn <- function(X, cellN, oversample, seed = NULL, k_extra = 10) {
  n <- nrow(X)
  if (n < 2) stop("Reduction matrix must have at least 2 cells")
  
  if (!is.null(seed)) set.seed(seed)
  
  # Defensive: ensure numeric matrix
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  finite_mask <- apply(X, 1, function(x) all(is.finite(x)))
  if (!all(finite_mask)) {
    warning(sprintf("Removing %d non-finite rows", sum(!finite_mask)))
    X <- X[finite_mask, , drop = FALSE]
    n <- nrow(X)
  }
  if (n == 0) stop("No usable cells in reduction matrix after filtering")
  
  # Determine number of metacells (centers)
  M <- max(1L, ceiling(n * oversample / cellN))
  
  # ---- 1. Select seed centers via farthest point sampling ----
  seeds <- .fps_seeds(X, n_seeds = M, seed = seed)
  seeds <- as.integer(seeds[seeds >= 1 & seeds <= n])
  centers <- X[seeds, , drop = FALSE]
  
  # ---- 2. Assign each cell to its nearest seed ----
  if (.safe_require("RANN")) {
    assign_idx <- RANN::nn2(data = centers, query = X, k = 1)$nn.idx[, 1]
  } else {
    D <- as.matrix(stats::dist(rbind(centers, X)))[seq_len(M), (M + 1):(M + n)]
    assign_idx <- max.col(-t(D))
  }
  
  # ---- 3. Group cells by nearest seed ----
  groups <- split(seq_len(n), assign_idx)
  
  # ---- 4. Expand each group to include local neighbors (~cellN size) ----
  k <- max(min(cellN + k_extra, n), 2)
  idx <- .knn_indices(X, k = k)
  
  groups <- lapply(groups, function(ix) {
    if (length(ix) >= cellN) {
      return(ix)
    }
    # pull neighbors of members until reaching target size
    extra <- unique(as.vector(idx[ix, , drop = FALSE]))
    unique(c(ix, extra))[1:min(k, length(unique(c(ix, extra))))]  # allow <cellN if sparse
  })
  
  # ---- 5. Defensive cleanup ----
  groups <- Filter(function(ix) length(ix) > 0, groups)
  all_cells <- seq_len(n)
  used_cells <- unique(unlist(groups))
  missing <- setdiff(all_cells, used_cells)
  if (length(missing) > 0) {
    message(sprintf("Adding %d unassigned cells as singleton metacells.", length(missing)))
    groups <- c(groups, as.list(missing))
  }
  
  message(sprintf(
    "Generated %d metacell groups (median group size = %.1f, mean = %.1f)",
    length(groups),
    median(vapply(groups, length, numeric(1))),
    mean(vapply(groups, length, numeric(1)))
  ))
  
  groups
}

# Create metacell groups via DBSCAN; refine to ~cellN members per metacell
.groups_dbscan <- function(X, cellN, oversample, eps = NULL, minPts = NULL, seed = NULL) {
  if (!.safe_require("dbscan")) stop("dbscan package required for mode='dbscan'")
  if (!is.null(seed)) set.seed(seed)
  if (is.null(eps)) eps <- .estimate_eps(X, cellN)
  if (is.null(minPts)) minPts <- max(ceiling(cellN * 0.6), 4)
  fit <- dbscan::dbscan(X, eps = eps, minPts = minPts, borderPoints = TRUE)
  labs <- fit$cluster
  cl_ids <- sort(unique(labs[labs > 0]))
  # Split/expand clusters into ~cellN using kNN seeds inside each cluster
  k <- max(cellN, 2)
  idx_all <- .knn_indices(X, k = k)
  groups <- list()
  for (cl in cl_ids) {
    ix <- which(labs == cl)
    if (length(ix) < ceiling(cellN/2)) next
    # number of metacells to carve from this cluster proportional to its size and oversample
    m_cl <- max(1L, round(length(ix) * oversample / cellN))
    seeds <- .fps_seeds(X[ix, , drop = FALSE], n_seeds = m_cl, seed = seed)
    for (s in seeds) {
      si <- ix[s]
      nn <- idx_all[si, ]
      nn <- nn[nn %in% ix]
      groups[[length(groups) + 1L]] <- nn[1:min(length(nn), k)]
    }
  }
  # If still short (e.g., many noise points), backfill using KNN seeding globally
  n <- nrow(X)
  target_M <- ceiling(n * oversample / cellN)
  if (length(groups) < target_M) {
    add <- target_M - length(groups)
    extra <- .groups_knn(X, cellN = cellN, oversample = add * cellN / n, seed = seed)
    groups <- c(groups, extra)
  }
  groups
}

# Create metacell groups via adaptive grid on first 2 dims, then split/merge to ~cellN
.groups_adaptive <- function(X, cellN, oversample, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  # bandwidth from median kNN radius
  k <- max(cellN, 4)
  idx <- .knn_indices(X, k)
  kth_dist <- vapply(seq_len(n), function(i) {
    j <- idx[i, k]
    sqrt(sum((X[i, ] - X[j, ])^2))
  }, numeric(1))
  h <- stats::median(kth_dist) / sqrt(2)
  # grid on first two dims (fallback to more dims by hashing)
  X2 <- if (ncol(X) >= 2) X[, 1:2, drop = FALSE] else cbind(X[, 1], 0)
  mins <- apply(X2, 2, min)
  bin <- function(v) floor((v - mins) / h)
  keys <- paste(bin(X2[, 1]), bin(X2[, 2]), sep = ":")
  bins <- split(seq_len(n), keys)
  groups <- list()
  # split large bins with kmeans; merge tiny bins to nearest neighbor bin centroid
  for (b in names(bins)) {
    ix <- bins[[b]]
    if (length(ix) <= 0) next
    if (length(ix) > 2 * cellN) {
      n_parts <- ceiling(length(ix) / cellN)
      km <- stats::kmeans(X[ix, , drop = FALSE], centers = n_parts, iter.max = 50)
      for (cc in seq_len(n_parts)) {
        members <- ix[km$cluster == cc]
        groups[[length(groups) + 1L]] <- members
      }
    } else {
      groups[[length(groups) + 1L]] <- ix
    }
  }
  # merge small groups
  cent <- function(ix) colMeans(X[ix, , drop = FALSE])
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    sizes <- vapply(groups, length, integer(1))
    small <- which(sizes < ceiling(0.5 * cellN))
    if (length(small) == 0) break
    for (si in small) {
      if (length(groups[[si]]) == 0) next
      c_si <- cent(groups[[si]])
      # find nearest other group
      d <- sapply(groups, function(g) if (length(g) == 0) Inf else sum((c_si - cent(g))^2))
      d[si] <- Inf
      ni <- which.min(d)
      if (is.finite(d[ni])) {
        groups[[ni]] <- c(groups[[ni]], groups[[si]])
        groups[[si]] <- integer(0)
        changed <- TRUE
      }
    }
    groups <- groups[vapply(groups, length, integer(1)) > 0]
  }
  # within-group trimming to ~cellN via seed + nearest neighbors
  idx_all <- .knn_indices(X, k = max(cellN, 2))
  out <- list()
  for (g in groups) {
    if (length(g) <= cellN) {
      out[[length(out) + 1L]] <- g
    } else {
      m_g <- max(1L, round(length(g) * oversample / cellN))
      seeds <- .fps_seeds(X[g, , drop = FALSE], n_seeds = m_g, seed = seed)
      for (s in seeds) {
        si <- g[s]
        nn <- idx_all[si, ]
        nn <- nn[nn %in% g]
        out[[length(out) + 1L]] <- nn[1:cellN]
      }
    }
  }
  # If we need more, backfill using global KNN
  target_M <- ceiling(n * oversample / cellN)
  if (length(out) < target_M) {
    add <- target_M - length(out)
    extra <- .groups_knn(X, cellN = cellN, oversample = add * cellN / n, seed = seed)
    out <- c(out, extra)
  }
  out
}


# =========================
# Public API
# =========================

#' @title makeMetacells
#' @description Generate metacells from a reduction embedding using kNN, DBSCAN, or an adaptive scheme.
#'
#' @param obj Amethyst single-cell object (S4) with @reductions and @metadata slots.
#' @param reduction Name of reduction in obj@reductions to use (e.g. "irlba").
#' @param cellN Target number of cells per metacell (default 50).
#' @param oversample Oversampling factor controlling number of metacells (default 1).
#' @param mode Metacell construction mode: "knn", "dbscan", or "adaptive".
#' @param dims Dimensions of the reduction to use (default 1:20).
#' @param threads Unused (reserved for future parallelization; currently must be 1).
#' @param seed Optional random seed for reproducibility.
#' @param k_extra Extra neighbors to pull in during kNN group expansion (mode="knn").
#' @param eps DBSCAN eps parameter (mode="dbscan"; if NULL, estimated from kNN distances).
#' @param minPts DBSCAN minPts parameter (mode="dbscan"; if NULL, derived from cellN).
#'
#' @return A new Amethyst object of the same class, with:
#'   \itemize{
#'     \item \code{metadata}: rows = metacells, including \code{originating_cellIDs}
#'     \item \code{reductions}: only a single reduction \code{mean_<reduction>} with metacell coordinates
#'     \item \code{genomeMatrices}: cleared for downstream metacell-specific matrices
#'   }
#' @export
#' @importFrom methods slot slotNames
#' @importFrom stats dist kmeans median as.dist hclust
makeMetacells <- function(
    obj,
    reduction,
    cellN = 50,
    oversample = 1,
    mode = c("knn", "dbscan", "adaptive"),
    dims = 1:20,
    threads = 1,
    seed = NULL,
    # mode-specific tuning knobs
    k_extra = 10,
    eps = NULL,
    minPts = NULL
) {
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)
  if (cellN < 2) stop("cellN must be >= 2")
  if (oversample <= 0) stop("oversample must be > 0")
  
  X <- .get_embedding(obj, reduction = reduction, dims = dims)
  md <- .get_metadata(obj)
  
  # align rows
  common <- intersect(rownames(X), rownames(md))
  if (length(common) == 0) stop("No overlapping cell IDs between reductions and metadata")
  X <- X[common, , drop = FALSE]
  md <- md[common, , drop = FALSE]
  
  # generate groups (as integer index vectors into X/md)
  groups <- switch(
    mode,
    knn      = .groups_knn(X, cellN = cellN, oversample = oversample, seed = seed, k_extra = k_extra),
    dbscan   = .groups_dbscan(X, cellN = cellN, oversample = oversample, eps = eps, minPts = minPts, seed = seed),
    adaptive = .groups_adaptive(X, cellN = cellN, oversample = oversample, seed = seed)
  )
  # Final defensive filter against empties
  groups <- Filter(function(ix) length(ix) > 0, groups)
  if (length(groups) == 0) stop("No metacell groups were formed; try adjusting parameters or mode")
  
  # Build metacell IDs
  M <- length(groups)
  metacell_ids <- sprintf("MC%04d", seq_len(M))
  
  # Aggregate metadata and reductions
  metacell_md  <- .aggregate_metadata(md, groups, metacell_ids)
  metacell_red <- .aggregate_reduction(X, groups, metacell_ids)
  
  # Construct new object of the same class
  new_obj <- obj
  methods::slot(new_obj, "metadata") <- metacell_md
  
  reds <- methods::slot(new_obj, "reductions")
  reds_new <- list()
  mean_reduction_name <- paste0("mean_", reduction)
  reds_new[[mean_reduction_name]] <- metacell_red
  methods::slot(new_obj, "reductions") <- reds_new
  
  # clear genomeMatrices if present (metacell-level matrices will be added downstream)
  if ("genomeMatrices" %in% slotNames(new_obj)) {
    methods::slot(new_obj, "genomeMatrices") <- list()
  }
  
  new_obj
}


#' @title metacellMeanReduction
#' @description Recompute a reduction on the metacell object by averaging the original reduction
#'   across member cells of each metacell.
#'
#' @param metacell_obj Metacell Amethyst object (from makeMetacells()).
#' @param original_obj Original single-cell Amethyst object (with reduction to average).
#' @param original_reduction Name of reduction in original_obj@reductions to average.
#' @param metacell_reduction Name to give the new reduction in metacell_obj@reductions.
#'
#' @return The metacell_obj with a new reduction added.
#' @export
#' @importFrom methods slot
metacellMeanReduction <- function(
    metacell_obj,
    original_obj,
    original_reduction,
    metacell_reduction
) {
  # sanity checks
  if (!"originating_cellIDs" %in% colnames(metacell_obj@metadata)) {
    stop("metacell_obj@metadata must contain 'originating_cellIDs'")
  }
  if (!"reductions" %in% slotNames(original_obj)) {
    stop("original_obj has no @reductions slot")
  }
  reds <- methods::slot(original_obj, "reductions")
  if (is.null(reds[[original_reduction]])) {
    stop(sprintf("Reduction '%s' not found in original_obj@reductions", original_reduction))
  }
  
  # get original reduction matrix
  red_mat <- as.matrix(reds[[original_reduction]])
  if (is.null(rownames(red_mat))) {
    stop("Original reduction matrix must have rownames as cell IDs")
  }
  
  # parse originating cell lists for each metacell
  cell_lists <- strsplit(
    metacell_obj@metadata$originating_cellIDs,
    ",",
    fixed = TRUE
  )
  
  # average coordinates for each metacell
  metacell_coords <- do.call(rbind, lapply(cell_lists, function(cells) {
    cells <- intersect(cells, rownames(red_mat))
    if (length(cells) == 0) {
      return(rep(NA_real_, ncol(red_mat)))
    }
    colMeans(red_mat[cells, , drop = FALSE], na.rm = TRUE)
  }))
  
  rownames(metacell_coords) <- rownames(metacell_obj@metadata)
  
  # attach new reduction matrix to metacell object
  metareds <- methods::slot(metacell_obj, "reductions")
  if (!is.list(metareds)) metareds <- list()
  metareds[[metacell_reduction]] <- metacell_coords
  methods::slot(metacell_obj, "reductions") <- metareds
  
  message(sprintf(
    "Added reduction '%s' to metacell object (%d metacells, %d dims)",
    metacell_reduction, nrow(metacell_coords), ncol(metacell_coords)
  ))
  
  metacell_obj
}


#' @title makeMetacellWindows
#' @description Build metacell-level genomic windows by aggregating raw H5 reads
#'   across the originating cells of each metacell (R-based; best for smaller datasets).
#'
#' @param original_obj Original Amethyst single-cell object with @h5paths, @index, and @metadata.
#' @param metacell_obj Metacell Amethyst object from makeMetacells().
#' @param stepsize Window size in bp (default 50000).
#' @param type Methylation context ("CG" or "CH"; default "CG").
#' @param metric Summary metric: "percent", "score", or "ratio" (default "score").
#' @param index Name of chromosome index in original_obj@index (default "chr_cg").
#' @param threads Parallel workers (via future/furrr) for per-metacell computation.
#' @param futureType Future backend: "multicore" or "multisession".
#' @param nmin Minimum number of total observations in a window to retain it.
#' @param chrList Optional vector of chromosomes to restrict to.
#' @param save Deprecated/unused (kept for backward compatibility).
#' @param matrix_name Name for the resulting genome matrix in metacell_obj@genomeMatrices.
#'
#' @return The metacell_obj with a new matrix added to @genomeMatrices.
#' @export
#' @importFrom rhdf5 h5read
#' @importFrom future plan sequential multicore multisession
#' @importFrom furrr future_map future_pmap furrr_options
#' @importFrom plyr round_any
#' @importFrom data.table data.table setorder setnames tstrsplit :=
#' @importFrom tibble column_to_rownames
#' @importFrom purrr map reduce
makeMetacellWindows <- function(
    original_obj,
    metacell_obj,
    stepsize = 50000,
    type = "CG",
    metric = "score",
    index = paste0("chr_", tolower(type)),
    threads = 1,
    futureType = "multicore",
    nmin = 2,
    chrList = NULL,
    save = FALSE,
    matrix_name = NULL
) {
  # ---------- sanity & setup ----------
  if (is.null(original_obj@h5paths$path) | is.null(original_obj@h5paths$barcode)) {
    stop("Please make sure original_obj@h5paths has valid 'path' and 'barcode' columns.")
  }
  if (is.null(original_obj@h5paths$prefix)) {
    original_obj@h5paths$prefix <- ""
  }
  
  if (is.null(original_obj@index[[index]])) {
    stop("Please index which rows in the H5 files correspond to each chromosome using indexChr.")
  }
  if (!metric %in% c("percent", "score", "ratio")) {
    stop("Metric must be one of: 'percent', 'score', 'ratio'.")
  }
  if (!"originating_cellIDs" %in% colnames(metacell_obj@metadata)) {
    stop("metacell_obj@metadata must contain 'originating_cellIDs'.")
  }
  
  # default output matrix name
  if (is.null(matrix_name)) {
    matrix_name <- paste0(tolower(type), "_", round(stepsize/1000), "k_", tolower(metric))
  }
  
  # threads
  if (threads > 1) {
    if (futureType == "multicore") {
      future::plan(future::multicore, workers = threads)
    } else if (futureType == "multisession") {
      future::plan(future::multisession, workers = threads)
    } else {
      stop("futureType must be 'multicore' or 'multisession' when threads > 1.")
    }
  }
  
  # ensure cell_id in h5paths matches index expectation
  if (!is.null(original_obj@h5paths$prefix)) {
    original_obj@h5paths$cell_id <- paste0(original_obj@h5paths$prefix, original_obj@h5paths$barcode)
  } else {
    original_obj@h5paths$cell_id <- original_obj@h5paths$barcode
  }
  
  # handy lookup
  h5tab <- original_obj@h5paths
  h5tab$cell_name <- paste0(ifelse(is.null(h5tab$prefix), "", h5tab$prefix),
                            sub("\\..*$", "", h5tab$barcode))
  
  # parse metacell membership
  mc_ids <- rownames(metacell_obj@metadata)
  members_list <- strsplit(metacell_obj@metadata$originating_cellIDs, ",", fixed = TRUE)
  names(members_list) <- mc_ids
  
  # helper: fetch global methylation for metacell (for score/ratio)
  get_metacell_global_m <- function(mc_id, context_type) {
    colname <- paste0("m", tolower(context_type), "_pct")
    # prefer metacell-level numeric metadata if present
    if (colname %in% colnames(metacell_obj@metadata)) {
      val <- metacell_obj@metadata[mc_id, colname]
      if (is.finite(val)) return(as.numeric(val)/100)
    }
    # fallback: mean of member-cell globals from original metadata
    mem <- members_list[[mc_id]]
    mem <- intersect(mem, rownames(original_obj@metadata))
    if (length(mem) == 0 || !(colname %in% colnames(original_obj@metadata))) return(NA_real_)
    mean(as.numeric(original_obj@metadata[mem, colname]), na.rm = TRUE)/100
  }
  
  # define chromosomes to iterate
  if (is.null(chrList)) {
    chr_groups <- as.list(names(original_obj@index[[index]]))
  } else {
    chr_groups <- chrList
  }
  
  by_chr <- list()
  
  # ---------- main per-chromosome loop ----------
  for (chr in chr_groups) {
    sites <- original_obj@index[[index]][[chr]]           # chr index for H5 file, keyed by cell_id
    if (is.null(sites) || nrow(sites) == 0) next
    
    # per-metacell aggregation
    mc_cols <- furrr::future_map(mc_ids, function(mc_id) {
      members <- members_list[[mc_id]]
      if (length(members) == 0) return(NULL)
      
      # map members to h5tab
      sel <- h5tab$cell_id %in% members | h5tab$cell_name %in% members | rownames(original_obj@metadata) %in% members
      subtab <- h5tab[sel, , drop = FALSE]
      if (nrow(subtab) == 0) return(NULL)
      
      mem_windows <- furrr::future_pmap(
        .l = list(subtab$path, subtab$barcode, subtab$prefix),
        .f = function(path, barcode, prefix) {
          tryCatch({
            cell_name <- paste0(ifelse(is.null(prefix), "", prefix), sub("\\..*$", "", barcode))
            
            # positions for this cell in this chromosome
            st <- sites$start[sites$cell_id == cell_name]
            ct <- sites$count[sites$cell_id == cell_name]
            if (length(st) == 0 || length(ct) == 0) return(NULL)
            
            h5 <- data.table::data.table(
              rhdf5::h5read(path, name = paste0(type, "/", barcode, "/1"),
                            start = st, count = ct)
            )
            h5[pos %% stepsize == 0, pos := pos + 1]
            h5[, window := paste0(chr, "_",
                                  plyr::round_any(pos, stepsize, floor), "_",
                                  plyr::round_any(pos, stepsize, ceiling))]
            h5[, .(m = sum(c != 0),
                   u = sum(t != 0),
                   n = sum(c + t, na.rm = TRUE)),
               by = window][n >= nmin, .(window, m, u)]
          }, error = function(e) {
            cat("Error processing H5 for barcode", barcode, ":", conditionMessage(e), "\n")
            return(NULL)
          })
        },
        .progress = FALSE
      )
      
      mem_windows <- Filter(Negate(is.null), mem_windows)
      if (length(mem_windows) == 0) return(NULL)
      
      for (i in seq_along(mem_windows)) {
        cn <- paste0("cell", i)
        data.table::setnames(mem_windows[[i]],
                             c("m", "u"),
                             c(paste0("m_", cn), paste0("u_", cn)))
      }
      
      mem_sum <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE), mem_windows)
      
      m_cols <- grep("^m_", colnames(mem_sum), value = TRUE)
      u_cols <- grep("^u_", colnames(mem_sum), value = TRUE)
      mem_sum$m <- rowSums(mem_sum[, ..m_cols], na.rm = TRUE)
      mem_sum$u <- rowSums(mem_sum[, ..u_cols], na.rm = TRUE)
      mem_sum <- mem_sum[, .(window, m, u)]
      
      denom <- mem_sum$m + mem_sum$u
      val <- ifelse(denom > 0, mem_sum$m / denom, NA_real_)
      
      if (metric != "percent") {
        meth_cell <- get_metacell_global_m(mc_id, type)
        if (!is.na(meth_cell)) {
          if (metric == "score") {
            val <- round(ifelse(val - meth_cell > 0,
                                (val - meth_cell) / (1 - meth_cell),
                                (val - meth_cell) / meth_cell), 3)
          } else if (metric == "ratio") {
            val <- round(val / meth_cell, 3)
          }
        } else {
          warning(sprintf("Global methylation (%s) missing for metacell %s; leaving raw proportion", tolower(type), mc_id))
        }
      } else {
        val <- val * 100
      }
      
      out <- data.table::data.table(window = mem_sum$window)
      out[[mc_id]] <- val
      out
    }, .options = furrr::furrr_options(seed = TRUE))
    
    mc_cols <- Filter(Negate(is.null), mc_cols)
    if (length(mc_cols) == 0) next
    mc_cols <- split(mc_cols, ceiling(seq_along(mc_cols)/1000))
    windows_merged <- Reduce(
      function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE),
      furrr::future_map(mc_cols, ~ Reduce(function(a, b)
        merge(a, b, by = "window", all = TRUE, sort = FALSE), .x))
    )
    
    by_chr[[chr]] <- windows_merged
    cat("\nCompleted ", chr, "\n")
  }
  
  if (length(by_chr) == 0) stop("No per-chromosome results produced")
  windows_merged <- data.table::rbindlist(by_chr, fill = TRUE, use.names = TRUE)
  
  windows_merged <- windows_merged[, c("chr", "start", "end") := {
    split_parts <- data.table::tstrsplit(window, "_", fixed = TRUE)
    list(split_parts[[1]], as.numeric(split_parts[[2]]), as.numeric(split_parts[[3]]))
  }]
  windows_merged <- windows_merged[(end - start) == stepsize]
  data.table::setorder(windows_merged, chr, start)
  windows_merged[, c("chr", "start", "end") := NULL]
  
  windows_merged <- windows_merged |>
    tibble::column_to_rownames(var = "window")
  windows_merged <- windows_merged[!sapply(rownames(windows_merged), function(name)
    length(strsplit(name, "_")[[1]]) > 3 || grepl("chrEBV|chrM|KI", name)), ]
  
  if (threads > 1) {
    future::plan(future::sequential)
    gc()
  }
  
  gm <- methods::slot(metacell_obj, "genomeMatrices")
  if (!is.list(gm)) gm <- list()
  gm[[matrix_name]] <- windows_merged
  methods::slot(metacell_obj, "genomeMatrices") <- gm
  
  message(sprintf("Added genome matrix '%s' (%d windows × %d metacells)",
                  matrix_name, nrow(windows_merged), ncol(windows_merged)))
  
  metacell_obj
}


#' @title metacellDims
#' @description Visualize metacell–cell relationships in UMAP / t-SNE space by
#'   connecting each cell to the centroid of its metacell.
#'
#' @param original_obj The original single-cell Amethyst object.
#' @param metacell_obj The metacell Amethyst object (from makeMetacells()).
#' @param reduction_orig Name of reduction in original_obj@reductions (e.g. "umap_irlba_100k_score").
#' @param reduction_meta Name of reduction in metacell_obj@reductions (e.g. "mean_umap_irlba_100k_score").
#' @param point_size Size of metacell points.
#' @param line_alpha Transparency of connecting lines.
#' @param line_width Width of connecting lines.
#' @param cell_alpha Transparency of original cell points.
#' @param cell_size Size of original cell points.
#' @param bg_color Background color for panel and plot.
#'
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_segment theme_void coord_cartesian theme element_rect
metacellDims <- function(
    original_obj,
    metacell_obj,
    reduction_orig,
    reduction_meta,
    point_size = 2.5,
    line_alpha = 0.1,
    line_width = 0.3,
    cell_alpha = 0.4,
    cell_size = 0.3,
    bg_color = "transparent"
) {
  if (is.null(original_obj@reductions[[reduction_orig]]))
    stop(paste0("Reduction '", reduction_orig, "' not found in original_obj."))
  if (is.null(metacell_obj@reductions[[reduction_meta]]))
    stop(paste0("Reduction '", reduction_meta, "' not found in metacell_obj."))
  
  orig_red <- as.data.frame(original_obj@reductions[[reduction_orig]])
  meta_red <- as.data.frame(metacell_obj@reductions[[reduction_meta]])
  orig_red$cell_id <- rownames(orig_red)
  meta_red$metacell_id <- rownames(meta_red)
  colnames(orig_red)[1:2] <- c("x", "y")
  colnames(meta_red)[1:2] <- c("x", "y")
  
  members_list <- strsplit(metacell_obj@metadata$originating_cellIDs, ",", fixed = TRUE)
  names(members_list) <- rownames(metacell_obj@metadata)
  
  map_df <- data.frame(
    cell_id = unlist(members_list),
    metacell_id = rep(names(members_list), lengths(members_list)),
    stringsAsFactors = FALSE
  )
  
  merged <- merge(orig_red, map_df, by = "cell_id")
  merged <- merge(merged, meta_red, by = "metacell_id",
                  suffixes = c("_orig", "_meta"))
  
  x_range <- range(c(merged$x_orig, merged$x_meta), na.rm = TRUE)
  y_range <- range(c(merged$y_orig, merged$y_meta), na.rm = TRUE)
  
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = merged,
      ggplot2::aes(x = x_orig, y = y_orig, xend = x_meta, yend = y_meta),
      color = "black", alpha = line_alpha, linewidth = line_width
    ) +
    ggplot2::geom_point(
      data = merged,
      ggplot2::aes(x = x_orig, y = y_orig),
      size = cell_size, alpha = cell_alpha, color = "grey60"
    ) +
    ggplot2::geom_point(
      data = meta_red,
      ggplot2::aes(x = x, y = y),
      size = point_size, color = "#A45EE5"
    ) +
    ggplot2::coord_cartesian(xlim = x_range, ylim = y_range, expand = TRUE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = bg_color, color = NA),
      plot.background  = ggplot2::element_rect(fill = bg_color, color = NA)
    )
}


#' @title calcMetacellSmoothedWindows
#' @description Aggregate c and t observations for metacells across smoothed genomic windows,
#'   computing one metacell at a time (memory-efficient).
#'
#' @param original_obj Amethyst single-cell object with valid @h5paths and @index.
#' @param metacell_obj Metacell object from makeMetacells() with metadata$originating_cellIDs.
#' @param type Methylation context ("CG" or "CH"; default "CG").
#' @param threads Parallel workers for H5 reading.
#' @param step Window width in bp.
#' @param smooth Rolling window size (odd integer; centered).
#' @param genome Genome identifier (currently unused; kept for compatibility).
#' @param index Name of chr index in @index (default "chr_cg").
#' @param futureType "multicore" or "multisession" for future.
#' @param returnSumMatrix Return c/t matrices (TRUE/FALSE).
#' @param returnPctMatrix Return %m matrix (TRUE/FALSE).
#' @param chrList Optional whitelist of chromosomes.
#' @param verbose Print progress.
#'
#' @return Either a list(sum_matrix, pct_matrix), or a single data.table depending on flags.
#' @export
#' @importFrom rhdf5 h5read
#' @importFrom future plan sequential multicore multisession
#' @importFrom purrr pmap
#' @importFrom plyr round_any
#' @importFrom data.table data.table as.data.table setorder tstrsplit :=
calcMetacellSmoothedWindows <- function(
    original_obj,
    metacell_obj,
    type = "CG",
    threads = 1,
    step = 500,
    smooth = 3,
    genome = "hg38",
    index = "chr_cg",
    futureType = "multicore",
    returnSumMatrix = TRUE,
    returnPctMatrix = TRUE,
    chrList = NULL,
    verbose = TRUE
) {
  Sys.setenv(OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1, HDF5_USE_FILE_LOCKING = "FALSE")
  
  if (is.null(original_obj@h5paths$path) || is.null(original_obj@h5paths$barcode))
    stop("original_obj@h5paths must have 'path' and 'barcode'.")
  
  if (is.null(original_obj@index[[index]]))
    stop("Chromosome index missing in original_obj@index[[index]]. Run indexChr() first.")
  
  h5paths <- data.table::as.data.table(original_obj@h5paths)
  if (is.null(h5paths$prefix)) h5paths[, prefix := ""]
  if (is.null(h5paths$cell_id)) h5paths[, cell_id := paste0(prefix, barcode)]
  
  if (is.null(metacell_obj@metadata$originating_cellIDs))
    stop("metacell_obj@metadata$originating_cellIDs is missing.")
  members_list <- strsplit(metacell_obj@metadata$originating_cellIDs, ",", fixed = TRUE)
  names(members_list) <- rownames(metacell_obj@metadata)
  metacell_ids <- names(members_list)
  
  chr_groups <- if (is.null(chrList)) names(original_obj@index[[index]]) else chrList
  
  if (threads > 1) {
    if (futureType == "multicore") future::plan(future::multicore, workers = threads)
    else if (futureType == "multisession") future::plan(future::multisession, workers = threads)
    else stop("futureType must be 'multicore' or 'multisession'.")
  }
  
  combined_counts <- NULL
  
  for (i in seq_along(metacell_ids)) {
    mc_id <- metacell_ids[i]
    if (verbose) cat("\n[Metacell", mc_id, "|", i, "of", length(metacell_ids), "] Starting...\n")
    
    members <- trimws(members_list[[mc_id]])
    members <- members[nzchar(members)]
    if (length(members) == 0) next
    
    m_idx <- match(members, h5paths$cell_id)
    valid <- which(!is.na(m_idx))
    if (length(valid) == 0) next
    
    m_paths   <- h5paths$path[m_idx[valid]]
    m_barc    <- h5paths$barcode[m_idx[valid]]
    m_prefix  <- h5paths$prefix[m_idx[valid]]
    m_cellids <- h5paths$cell_id[m_idx[valid]]
    
    by_chr <- list()
    
    for (chr in chr_groups) {
      sites <- data.table::as.data.table(original_obj@index[[index]][[chr]])
      sites <- sites[cell_id %in% m_cellids]
      if (nrow(sites) == 0) next
      
      if (verbose) cat("  [", chr, "] Reading...", sep = "")
      
      member_results <- purrr::pmap(
        .l = list(m_paths, m_barc, m_prefix, m_cellids),
        .f = function(path, barcode, prefix, cell_id) {
          tryCatch({
            st <- sites$start[sites$cell_id == cell_id]
            ct <- sites$count[sites$cell_id == cell_id]
            if (length(st) == 0L || length(ct) == 0L || is.na(st) || is.na(ct) || ct <= 0L)
              return(data.table::data.table(window = character(), c = numeric(), t = numeric()))
            
            dt <- rhdf5::h5read(path, name = paste0(type, "/", barcode, "/1"),
                                start = st, count = ct)
            dt <- data.table::as.data.table(dt)
            if (!all(c("chr", "pos", "c", "t") %in% names(dt)))
              return(data.table::data.table(window = character(), c = numeric(), t = numeric()))
            
            dt[pos %% step == 0L, pos := pos + 1L]
            dt[, window := paste0(chr, "_",
                                  plyr::round_any(pos, step, f = floor), "_",
                                  plyr::round_any(pos, step, f = ceiling))]
            dt <- dt[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = window]
            data.table::as.data.table(dt)
          }, error = function(e) {
            if (verbose) cat("\n    Error in cell", barcode, ":", conditionMessage(e), "\n")
            data.table::data.table(window = character(), c = numeric(), t = numeric())
          })
        },
        .progress = FALSE
      )
      
      member_results <- lapply(member_results, data.table::as.data.table)
      member_results <- member_results[vapply(member_results, nrow, integer(1)) > 0]
      if (length(member_results) == 0) next
      
      merged_chr <- data.table::rbindlist(member_results, use.names = TRUE, fill = TRUE)
      merged_chr <- merged_chr[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = window]
      data.table::setnames(merged_chr, c("c", "t"), c(paste0(mc_id, "_c"), paste0(mc_id, "_t")))
      by_chr[[chr]] <- merged_chr
      if (verbose) cat(" Done.\n")
    }
    
    if (length(by_chr) == 0) next
    mc_counts <- data.table::rbindlist(by_chr, fill = TRUE, use.names = TRUE)
    
    sp <- data.table::tstrsplit(mc_counts$window, "_", fixed = TRUE)
    mc_counts[, chr := sp[[1]]]
    mc_counts[, start := suppressWarnings(as.numeric(sp[[2]]))]
    mc_counts[, end := suppressWarnings(as.numeric(sp[[3]]))]
    data.table::setorder(mc_counts, chr, start)
    
    if (is.null(combined_counts)) {
      combined_counts <- mc_counts
    } else {
      combined_counts <- merge(combined_counts, mc_counts, by = c("window", "chr", "start", "end"), all = TRUE, sort = FALSE)
    }
    
    if (verbose) cat("  [Metacell", mc_id, "] completed and appended.\n")
  }
  
  if (threads > 1) future::plan(future::sequential)
  
  ct_cols <- grep("(_c|_t)$", names(combined_counts), value = TRUE)
  if (length(ct_cols) > 0 && smooth > 1) {
    combined_counts[, (ct_cols) := lapply(.SD, function(x) {
      data.table::frollsum(x, n = smooth, align = "center", na.rm = TRUE)
    }), by = chr, .SDcols = ct_cols]
  }
  
  pct_dt <- NULL
  if (returnPctMatrix) {
    pct_dt <- data.table::copy(combined_counts)
    mc_ids2 <- unique(sub("_(c|t)$", "", ct_cols))
    mc_ids2 <- mc_ids2[sapply(mc_ids2, function(m) all(paste0(m, c("_c", "_t")) %in% names(pct_dt)))]
    for (m in mc_ids2) {
      c_col <- paste0(m, "_c")
      t_col <- paste0(m, "_t")
      pct_dt[, (m) := {
        num <- get(c_col); den <- get(c_col) + get(t_col)
        ifelse(den > 0, round(100 * num / den, 2), NA_real_)
      }]
    }
    pct_keep <- c("window", "chr", "start", "end", mc_ids2)
    pct_dt <- pct_dt[, ..pct_keep]
  }
  
  if (returnSumMatrix & returnPctMatrix) {
    sum_keep <- c("window", "chr", "start", "end", ct_cols)
    list(sum_matrix = combined_counts[, ..sum_keep],
         pct_matrix = pct_dt)
  } else if (returnPctMatrix) {
    pct_dt
  } else {
    combined_counts
  }
}


#' @title writeMetacellJobJSON
#' @description Write a JSON job description for Python metacell window generation.
#'
#' @param original_obj Amethyst object with @h5paths and @ref.
#' @param metacell_obj Metacell Amethyst object (from makeMetacells()) with
#'   metadata$originating_cellIDs (Amethyst cell IDs, i.e. prefix + barcode).
#' @param out_json Path to JSON file to write.
#' @param type Methylation type, e.g. "CG" (default "CG").
#' @param step Window size in bp (default 1000).
#' @param smooth Smoothing width in windows (default 3).
#' @param output_h5 Path to output HDF5 file that Python will write.
#'
#' @return Invisibly returns the job list; writes JSON to disk.
#' @export
#' @importFrom data.table as.data.table
#' @importFrom jsonlite write_json
writeMetacellJobJSON <- function(
    original_obj,
    metacell_obj,
    out_json,
    type      = "CG",
    step      = 1000,
    smooth    = 3,
    output_h5 = "metacell_windows.h5"
) {
  if (is.null(original_obj@ref))
    stop("@ref slot is missing; load gene/region reference first.")
  
  ref <- data.table::as.data.table(original_obj@ref)
  
  chr_col <- NULL
  if ("seqid" %in% names(ref)) {
    chr_col <- "seqid"
  } else if ("chr" %in% names(ref)) {
    chr_col <- "chr"
  } else {
    stop("@ref slot must contain either 'seqid' or 'chr' column.")
  }
  
  if (!"end" %in% names(ref)) {
    stop("@ref slot must contain 'end' column for chromosome sizes.")
  }
  
  chr_sizes_dt <- ref[, .(max_end = max(end, na.rm = TRUE)), by = chr_col]
  chrom_sizes_list <- as.list(chr_sizes_dt$max_end)
  names(chrom_sizes_list) <- chr_sizes_dt[[chr_col]]
  
  h5paths_dt <- data.table::as.data.table(original_obj@h5paths)
  
  if (is.null(h5paths_dt$barcode) || is.null(h5paths_dt$path)) {
    stop("@h5paths must contain 'barcode' and 'path' columns.")
  }
  
  if (is.null(h5paths_dt$prefix)) {
    h5paths_dt[, prefix := ""]
  }
  
  if (is.null(h5paths_dt$cell_ids)) {
    h5paths_dt[, cell_ids := paste0(prefix, barcode)]
  }
  
  if (is.null(metacell_obj@metadata$originating_cellIDs)) {
    stop("metacell_obj@metadata$originating_cellIDs is missing.")
  }
  
  members_str <- metacell_obj@metadata$originating_cellIDs
  if (is.null(names(members_str))) {
    names(members_str) <- rownames(metacell_obj@metadata)
  }
  
  metacells_list <- lapply(names(members_str), function(mc_id) {
    cell_ids_am <- trimws(unlist(strsplit(members_str[[mc_id]], ",", fixed = TRUE)))
    cell_ids_am <- cell_ids_am[nzchar(cell_ids_am)]
    
    if (length(cell_ids_am) == 0L) {
      return(list())
    }
    
    m_idx <- match(cell_ids_am, h5paths_dt$cell_ids)
    
    if (any(is.na(m_idx))) {
      missing_ids <- cell_ids_am[is.na(m_idx)]
      warning(sprintf(
        "Some originating_cellIDs for metacell '%s' were not found in original_obj@h5paths$cell_ids. Example: %s",
        mc_id,
        paste(head(missing_ids, 5), collapse = ", ")
      ))
    }
    sub_dt <- h5paths_dt[m_idx[!is.na(m_idx)], ]
    
    if (nrow(sub_dt) == 0L) {
      warning(sprintf("Metacell '%s' has no valid members with H5 mapping; it will be empty.", mc_id))
      return(list())
    }
    
    sub_dt <- unique(sub_dt, by = c("barcode", "path"))
    
    lapply(seq_len(nrow(sub_dt)), function(i) {
      list(
        barcode = sub_dt$barcode[i],
        path    = sub_dt$path[i]
      )
    })
  })
  
  names(metacells_list) <- names(members_str)
  
  job <- list(
    type        = type,
    step        = step,
    smooth      = smooth,
    chrom_sizes = chrom_sizes_list,
    output_h5   = output_h5,
    metacells   = metacells_list
  )
  
  jsonlite::write_json(job, out_json, pretty = TRUE, auto_unbox = TRUE)
  message("✓ Wrote job JSON: ", out_json)
  
  invisible(job)
}


# =========================
# H5 helper utilities
# =========================

.metacell_h5_index_cache <- new.env(parent = emptyenv())

#' @title buildMetacellH5Index
#' @description Build or retrieve a cached index for a metacell H5 file.
#'
#' @param h5file Path to metacell H5 (e.g. "metacell_cg_1kbp.h5").
#' @param type H5 group name (e.g. "CG").
#' @param use_cache Logical; reuse cached index if available (default TRUE).
#'
#' @return data.table with columns: idx, chr, start, end (keyed by chr,start,end).
#' @export
#' @importFrom rhdf5 h5read
#' @importFrom data.table data.table setkey
buildMetacellH5Index <- function(h5file, type = "CG", use_cache = TRUE) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package 'rhdf5' is required for buildMetacellH5Index().")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for buildMetacellH5Index().")
  }
  
  key <- paste0(normalizePath(h5file, winslash = "/", mustWork = FALSE), "::", type)
  if (use_cache && exists(key, envir = .metacell_h5_index_cache, inherits = FALSE)) {
    return(get(key, envir = .metacell_h5_index_cache))
  }
  
  chr   <- rhdf5::h5read(h5file, paste0("/", type, "/chr"))
  start <- rhdf5::h5read(h5file, paste0("/", type, "/start"))
  end   <- rhdf5::h5read(h5file, paste0("/", type, "/end"))
  
  idx_dt <- data.table::data.table(
    idx   = seq_along(chr),
    chr   = as.character(chr),
    start = as.numeric(start),
    end   = as.numeric(end)
  )
  data.table::setkey(idx_dt, chr, start, end)
  
  if (use_cache) {
    assign(key, idx_dt, envir = .metacell_h5_index_cache)
  }
  idx_dt
}

#' @title .get_ref_dt_and_chr_col
#' @description Internal helper to coerce obj@ref to data.table and detect chromosome column.
#' @keywords internal
#' @importFrom data.table as.data.table
.get_ref_dt_and_chr_col <- function(obj) {
  if (is.null(obj@ref)) {
    stop("obj@ref is NULL; please run makeRef() on this object first.")
  }
  ref <- obj@ref
  if (!"data.table" %in% class(ref)) {
    ref <- data.table::as.data.table(ref)
  }
  
  if ("chr" %in% names(ref)) {
    chr_col <- "chr"
  } else if ("seqid" %in% names(ref)) {
    chr_col <- "seqid"
  } else {
    stop("obj@ref must contain either a 'chr' or 'seqid' column.")
  }
  
  list(ref = ref, chr_col = chr_col)
}

#' @title getRangesFromGenes
#' @description Convert a set of gene names into genomic ranges using obj@ref.
#'
#' @param ref data.table version of obj@ref.
#' @param chr_col Name of chromosome column in ref ("chr" or "seqid").
#' @param genes Character vector of gene symbols.
#' @param trackOverhang bp to extend on each side of gene coordinates.
#'
#' @return data.table with columns: gene_name, chrom, range_start, range_end.
#' @export
#' @importFrom data.table data.table
getRangesFromGenes <- function(ref, chr_col, genes, trackOverhang = 5000) {
  if (!"gene_name" %in% names(ref)) {
    stop("obj@ref must contain 'gene_name' to use genes=.") 
  }
  if (!"type" %in% names(ref) || !"start" %in% names(ref) || !"end" %in% names(ref)) {
    stop("obj@ref must contain 'type', 'start', and 'end' columns.")
  }
  
  genes <- unique(genes)
  ref_subset <- ref[gene_name %in% genes & type == "gene"]
  if (nrow(ref_subset) == 0) {
    stop("None of the requested genes were found in obj@ref$type == 'gene'.")
  }
  
  ref_best <- ref_subset[, .SD[which.max(end - start)], by = gene_name]
  ref_best[, chrom       := get(chr_col)]
  ref_best[, range_start := pmax(0L, as.integer(start) - trackOverhang)]
  ref_best[, range_end   := as.integer(end) + trackOverhang]
  
  ref_best[, .(gene_name, chrom, range_start, range_end)]
}

#' @title getRangesFromRegions
#' @description Parse region strings "chr_start_end" into genomic ranges.
#'
#' @param regions Character vector of "chr_start_end".
#'
#' @return data.table with columns: region_id, chrom, range_start, range_end.
#' @export
#' @importFrom data.table data.table tstrsplit
getRangesFromRegions <- function(regions) {
  if (length(regions) == 0) {
    stop("regions vector is empty.")
  }
  parts <- data.table::tstrsplit(regions, "_", fixed = TRUE)
  if (length(parts) != 3) {
    stop("Each region must be of the form 'chr_start_end'.")
  }
  dt <- data.table::data.table(
    region_id   = regions,
    chrom       = parts[[1]],
    range_start = as.numeric(parts[[2]]),
    range_end   = as.numeric(parts[[3]])
  )
  dt
}

#' @title getWindowIdxForRanges
#' @description Determine which window indices from an H5 index overlap requested ranges.
#'
#' @param idx_dt data.table from buildMetacellH5Index().
#' @param ranges_dt data.table with columns: chrom, range_start, range_end.
#'
#' @return Integer vector of unique window indices.
#' @export
getWindowIdxForRanges <- function(idx_dt, ranges_dt) {
  keep_idx <- integer(0L)
  for (i in seq_len(nrow(ranges_dt))) {
    chrom <- ranges_dt$chrom[i]
    rs    <- ranges_dt$range_start[i]
    re    <- ranges_dt$range_end[i]
    this  <- idx_dt[chr == chrom & start >= rs & end <= re, idx]
    if (length(this)) {
      keep_idx <- c(keep_idx, this)
    }
  }
  sort(unique(keep_idx))
}

#' @title readH5WindowsSubset
#' @description Read a subset of windows from the H5 dataset and return a matrix [n_windows x n_metacells].
#'
#' @param h5file Path to H5.
#' @param type H5 group name (e.g. "CG").
#' @param which Dataset name inside group ("pct", "c", or "t").
#' @param keep_idx Integer vector of window indices to read (1-based).
#' @param n_metacells Number of metacells (nrow(obj@metadata)).
#'
#' @return Numeric matrix with nrow = length(keep_idx), ncol = n_metacells.
#' @export
#' @importFrom rhdf5 h5read
readH5WindowsSubset <- function(h5file, type, which, keep_idx, n_metacells) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package 'rhdf5' is required for readH5WindowsSubset().")
  }
  if (length(keep_idx) == 0) {
    stop("keep_idx is empty; nothing to read from H5.")
  }
  keep_idx <- sort(unique(as.integer(keep_idx)))
  
  path <- paste0("/", type, "/", which)
  
  probe1 <- rhdf5::h5read(h5file, path, index = list(1L, NULL))
  len1   <- length(probe1)
  if (len1 == n_metacells) {
    orientation <- "windows_first"
  } else {
    probe2 <- rhdf5::h5read(h5file, path, index = list(NULL, 1L))
    len2   <- length(probe2)
    if (len2 == n_metacells) {
      orientation <- "metacells_first"
    } else {
      stop(
        sprintf(
          "Cannot infer H5 orientation: length(probe1)=%d, length(probe2)=%d, n_metacells=%d",
          len1, len2, n_metacells
        )
      )
    }
  }
  
  if (orientation == "windows_first") {
    sub <- rhdf5::h5read(h5file, path, index = list(keep_idx, NULL))
    mat <- as.matrix(sub)
  } else {
    sub <- rhdf5::h5read(h5file, path, index = list(NULL, keep_idx))
    mat <- t(as.matrix(sub))
  }
  
  mat
}


#' @title populateMetacellTrackFromH5
#' @description Populate or extend a metacell track from a Python-generated H5 file,
#'   without loading the entire H5 matrix into memory.
#'
#' @param obj Amethyst metacell object (from makeMetacells()).
#' @param h5file Path to metacell H5 file (from Python script).
#' @param type H5 group name (default "CG").
#' @param which Matrix name inside group ("pct" for %m, or "c"/"t").
#' @param track Name of track in obj@tracks to populate (e.g. "cg_1kbp_pct").
#' @param genes Optional character vector of gene names.
#' @param regions Optional character vector of "chr_start_end" strings.
#' @param trackOverhang bp to extend on each side of gene coordinates.
#' @param verbose Print progress messages.
#' @param force Logical; if TRUE, overwrite any existing windows in the requested regions.
#'
#' @return Modified obj with obj@tracks[[track]] containing the requested windows.
#' @export
#' @importFrom data.table data.table as.data.table setorder
populateMetacellTrackFromH5 <- function(
    obj,
    h5file,
    type         = "CG",
    which        = "pct",
    track,
    genes        = NULL,
    regions      = NULL,
    trackOverhang = 5000,
    verbose      = TRUE,
    force        = FALSE
) {
  if (missing(track) || is.null(track) || !nzchar(track)) {
    stop("Please provide a non-empty 'track' name.")
  }
  if (is.null(genes) && is.null(regions)) {
    stop("Specify at least one of 'genes' or 'regions'.")
  }
  
  ref_info <- .get_ref_dt_and_chr_col(obj)
  ref     <- ref_info$ref
  chr_col <- ref_info$chr_col
  
  if (!is.null(genes)) {
    if (verbose) {
      message("Converting genes to genomic ranges...")
    }
    gr_genes <- getRangesFromGenes(ref, chr_col, genes, trackOverhang = trackOverhang)
    ranges_dt <- data.table::data.table(
      chrom       = gr_genes$chrom,
      range_start = gr_genes$range_start,
      range_end   = gr_genes$range_end
    )
  } else {
    if (verbose) {
      message("Using provided regions as genomic ranges...")
    }
    gr_regions <- getRangesFromRegions(regions)
    ranges_dt <- gr_regions[, .(chrom, range_start, range_end)]
  }
  
  if (verbose) {
    message("Building (or retrieving cached) H5 index...")
  }
  idx_dt <- buildMetacellH5Index(h5file, type = type, use_cache = TRUE)
  
  if (verbose) {
    message("Determining windows overlapping requested ranges...")
  }
  keep_idx <- getWindowIdxForRanges(idx_dt, ranges_dt)
  if (length(keep_idx) == 0) {
    warning("No windows in the H5 index overlap the requested ranges.")
    return(obj)
  }
  if (verbose) {
    message(sprintf("Will read %d windows from H5.", length(keep_idx)))
  }
  
  existing <- NULL
  if (!is.null(obj@tracks[[track]]) && !force) {
    existing <- obj@tracks[[track]]
    if (!"data.table" %in% class(existing)) {
      existing <- data.table::as.data.table(existing)
    }
    
    if (!all(c("chr", "start", "end") %in% names(existing))) {
      warning(sprintf("Existing track '%s' is missing chr/start/end; treating as empty.", track))
      existing <- NULL
      
    } else {
      existing[, window := paste(chr, start, end, sep = "_")]
      idx_dt[, window := paste(chr, start, end, sep = "_")]
      
      add_windows <- idx_dt[idx %in% keep_idx, unique(window)]
      already_have <- existing$window
      new_windows <- setdiff(add_windows, already_have)
      
      if (length(new_windows) == 0) {
        if (verbose) {
          message("All requested windows already present in existing track; nothing to add.")
        }
        existing[, window := NULL]
        idx_dt[, window := NULL]
        return(obj)
      }
      
      keep_idx <- idx_dt[window %in% new_windows, idx]
      keep_idx <- sort(unique(keep_idx))
      
      if (verbose) {
        message(sprintf("After excluding existing windows, %d new windows remain to fetch.", length(keep_idx)))
      }
      
      existing[, window := NULL]
      idx_dt[, window := NULL]
    }
    
  } else if (force && !is.null(obj@tracks[[track]]) && verbose) {
    message("force = TRUE: existing windows will be overwritten for this region.")
    existing <- obj@tracks[[track]]
    if (!"data.table" %in% class(existing)) {
      existing <- data.table::as.data.table(existing)
    }
  }
  
  n_metacells <- nrow(obj@metadata)
  if (verbose) {
    message(sprintf("Reading %s/%s subset from H5...", type, which))
  }
  mat_subset <- readH5WindowsSubset(
    h5file      = h5file,
    type        = type,
    which       = which,
    keep_idx    = keep_idx,
    n_metacells = n_metacells
  )
  
  idx_sub <- idx_dt[idx %in% keep_idx]
  data.table::setorder(idx_sub, idx)
  
  if (nrow(idx_sub) != nrow(mat_subset)) {
    stop("Row count mismatch between idx_sub and mat_subset; indexing error.")
  }
  
  track_new <- data.table::data.table(
    window = paste(idx_sub$chr, idx_sub$start, idx_sub$end, sep = "_"),
    chr    = idx_sub$chr,
    start  = idx_sub$start,
    end    = idx_sub$end
  )
  
  mc_ids <- rownames(obj@metadata)
  if (is.null(mc_ids) || length(mc_ids) != ncol(mat_subset)) {
    warning("rownames(obj@metadata) mismatch; using generic MC IDs.")
    mc_ids <- paste0("MC", sprintf("%04d", seq_len(ncol(mat_subset))))
  }
  mat_dt <- data.table::as.data.table(mat_subset)
  data.table::setnames(mat_dt, mc_ids)
  
  track_new <- cbind(track_new, mat_dt)
  
  if (!is.null(existing)) {
    existing <- data.table::as.data.table(existing)
    
    all_cols <- union(names(existing), names(track_new))
    for (cc in setdiff(all_cols, names(existing))) existing[[cc]]     <- NA_real_
    for (cc in setdiff(all_cols, names(track_new))) track_new[[cc]] <- NA_real_
    
    existing  <- existing[, ..all_cols]
    track_new <- track_new[, ..all_cols]
    
    combined <- data.table::rbindlist(list(existing, track_new), use.names = TRUE, fill = TRUE)
    combined <- combined[!duplicated(combined[, .(chr, start, end)])]
    data.table::setorder(combined, chr, start)
    
    obj@tracks[[track]] <- combined
  } else {
    data.table::setorder(track_new, chr, start)
    obj@tracks[[track]] <- track_new
  }
  
  if (verbose) {
    message(sprintf("Track '%s' now has %d windows and %d metacell columns.",
                    track,
                    nrow(obj@tracks[[track]]),
                    ncol(obj@tracks[[track]]) - 4L))
  }
  
  obj
}


#' @title metaMap
#' @description Plot methylation levels for metacells aggregated across genes or regions,
#'   with optional hierarchical clustering of metacell rows (local or global) and gene models.
#'
#' @param obj Amethyst object containing metacell methylation matrices in @tracks and a gene annotation in @ref.
#' @param track Name of the matrix in @tracks containing methylation values.
#' @param genes Gene list to plot.
#' @param regions Optional genomic regions to plot (e.g. "chr1_1000_2000").
#' @param colors Optional custom color palette.
#' @param trackOverhang Number of bp to extend beyond the gene on each side.
#' @param arrowOverhang Number of bp arrow extends beyond gene body.
#' @param legend Logical; include legend.
#' @param removeNA Logical; drop non-finite values before plotting.
#' @param trackScale Scale for gene track height.
#' @param arrowScale Manual arrow size (NULL = auto).
#' @param colorMax Upper bound on color scale (NULL = auto).
#' @param order Optional manual order of metacell rows.
#' @param nrow Number of rows in multi-gene layout.
#' @param clusterRows Logical; cluster metacell rows.
#' @param clusterMethod Linkage method for hclust.
#' @param clusterScope "local" (per gene) or "global" (shared order).
#' @param extraTracks Optional named list of coordinate tracks to show below gene body.
#' @param extraTrackColors Optional colors for extraTracks.
#' @param remove Optional regex to remove genes from obj@ref by name.
#'
#' @return A grid of ggplot objects arranged with gridExtra::grid.arrange().
#' @export
#' @importFrom data.table as.data.table tstrsplit
#' @importFrom dplyr filter distinct group_by mutate arrange slice_head rowwise
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_tile geom_rect geom_segment geom_text scale_fill_gradientn ylab xlab theme element_blank arrow
#' @importFrom grid unit
#' @importFrom gridExtra grid.arrange
#' @importFrom scales squish
metaMap <- function(obj,
                    track = NULL,
                    genes = NULL,
                    regions = NULL,
                    colors = NULL,
                    trackOverhang = 5000,
                    arrowOverhang = 3000,
                    legend = TRUE,
                    removeNA = TRUE,
                    trackScale = .07,
                    arrowScale = NULL,
                    colorMax = NULL,
                    order = NULL,
                    nrow = max(length(genes), length(regions)),
                    clusterRows = TRUE,
                    clusterMethod = "ward.D2",
                    clusterScope = "local",
                    extraTracks = NULL,
                    extraTrackColors = NULL,
                    remove = NULL) {
  
  if (is.null(obj@ref)) {
    stop("Please make sure a genome annotation file has been added to the obj@ref slot with makeRef().")
  }
  if (!"data.table" %in% class(obj@ref)) {
    obj@ref <- data.table::as.data.table(obj@ref)
  }
  
  if (!is.null(remove)) {
    obj@ref <- obj@ref[!grepl(x = obj@ref$gene_name, pattern = remove), ]
  }
  
  if (!is.null(extraTracks)) {
    for (i in seq_along(extraTracks)) {
      if (!(all(c("chr", "start", "end") %in% names(extraTracks[[i]])))) {
        stop("Each extraTrack must contain chr, start, and end columns.")
      }
      if (!(is.numeric(extraTracks[[i]]$start) && is.numeric(extraTracks[[i]]$end))) {
        stop("All start/end columns in extraTracks must be numeric.")
      }
    }
  }
  
  if (!is.null(regions)) {
    ranges <- regions
  }
  if (!is.null(genes)) {
    genes <- unique(genes)
    ref_subset <- obj@ref[gene_name %in% genes & type %in% c("gene", "exon")]
    ranges <- ref_subset[type == "gene", .SD[which.max(end - start)], by = gene_name][order(match(gene_name, genes))]$location
  }
  
  if (clusterRows && clusterScope == "global" && ncol(obj@tracks[[track]]) > 4) {
    all_values <- obj@tracks[[track]]
    numeric_cols <- setdiff(colnames(all_values), c("window", "chr", "start", "end"))
    mat_for_clust <- as.matrix(all_values[, ..numeric_cols])
    mat_for_clust[is.na(mat_for_clust)] <- 0
    d <- stats::as.dist(1 - stats::cor(mat_for_clust, use = "pairwise.complete.obs"))
    hc_global <- stats::hclust(d, method = clusterMethod)
    clustered_order_global <- hc_global$labels[hc_global$order]
  } else {
    clustered_order_global <- NULL
  }
  
  p <- vector("list", length(ranges))
  for (i in seq_along(ranges)) {
    split_parts <- data.table::tstrsplit(ranges[[i]], "_", fixed = TRUE)
    chrom <- split_parts[[1]]
    min <- as.numeric(split_parts[[2]]) - trackOverhang
    max <- as.numeric(split_parts[[3]]) + trackOverhang
    
    values <- obj@tracks[[track]][(chr == chrom & start >= min & end <= max)]
    ngroups <- ncol(values) - 3
    
    if (clusterRows && clusterScope == "local" && ncol(values) > 4) {
      numeric_cols <- setdiff(colnames(values), c("window", "chr", "start", "end"))
      mat_for_clust <- as.matrix(values[, ..numeric_cols])
      mat_for_clust[is.na(mat_for_clust)] <- 0
      d <- stats::as.dist(1 - stats::cor(mat_for_clust, use = "pairwise.complete.obs"))
      hc <- stats::hclust(d, method = clusterMethod)
      clustered_order <- hc$labels[hc$order]
      values <- values[, c(c("window", "chr", "start", "end"), clustered_order), with = FALSE]
    }
    
    if (clusterRows && clusterScope == "global" && !is.null(clustered_order_global)) {
      numeric_cols <- setdiff(colnames(values), c("window", "chr", "start", "end"))
      matched_cols <- intersect(clustered_order_global, numeric_cols)
      values <- values[, c(c("window", "chr", "start", "end"), matched_cols), with = FALSE]
    }
    
    if (!is.null(regions)) {
      ref <- obj@ref |>
        dplyr::filter(seqid == chrom & start >= min & end <= max & type %in% c("gene", "exon")) |>
        dplyr::distinct(gene_name, start, end, type, strand) |>
        dplyr::group_by(gene_name) |>
        dplyr::mutate(label = dplyr::cur_group_id()) |>
        dplyr::mutate(trackHeight = ngroups * label * trackScale,
                      ymin = -(trackHeight - (ngroups * 0.03)),
                      ymax = -(trackHeight + (ngroups * 0.03)))
    }
    
    if (!is.null(genes)) {
      ref <- ref_subset |>
        dplyr::filter(gene_name == genes[i] & type %in% c("gene", "exon")) |>
        dplyr::distinct(gene_name, start, end, type, strand) |>
        dplyr::mutate(trackHeight = ngroups * trackScale,
                      ymin = -(trackHeight - (ngroups * 0.03)),
                      ymax = -(trackHeight + (ngroups * 0.03)))
    }
    
    genenames <- ref |> dplyr::group_by(gene_name) |> 
      dplyr::arrange(dplyr::desc(end)) |> dplyr::slice_head(n = 1) |>
      dplyr::mutate(text_start = ifelse(strand == "+", (end + 2500), (end + 1000)))
    genebody <- ref |> dplyr::filter(type == "gene") |>
      dplyr::mutate(x = ifelse(strand == "+", start, end),
                    xend = ifelse(strand == "+", end + arrowOverhang, start - arrowOverhang))
    promoters <- ref |> dplyr::filter(type == "gene") |>
      dplyr::mutate(promoter_start = ifelse(strand == "+", (start - 1500), (end - 1500)),
                    promoter_end = ifelse(strand == "+", (start + 1500), (end + 1500)))
    exons <- ref |> dplyr::filter(type == "exon")
    
    if (is.null(arrowScale)) arrowHeight <- (5 / (ngroups * nrow(genenames)))
    else arrowHeight <- arrowScale
    
    values <- tidyr::pivot_longer(
      values,
      cols = setdiff(colnames(values), c("window", "chr", "start", "end")),
      names_to = "group",
      values_to = "pct_m"
    ) |>
      dplyr::rowwise() |>
      dplyr::mutate(middle = mean(c(start, end), na.rm = TRUE))
    
    if (removeNA) values <- values |> dplyr::filter(is.finite(pct_m))
    
    if (!is.null(order)) values$group <- factor(values$group, levels = order)
    else values$group <- factor(values$group, levels = unique(values$group))
    
    trackHeight <- mean(ref$trackHeight, na.rm = TRUE)
    width <- mean(values$end - values$start, na.rm = TRUE)
    
    if (is.null(colors)) {
      if (mean(values$pct_m, na.rm = TRUE) > 10) {
        pal <- c("#005eff", "#9daec9", "#cccccc", "#dbdbdb")
      } else {
        pal <- c("black", "red", "yellow")
      }
    } else pal <- colors
    
    if (is.null(colorMax)) {
      if (mean(values$pct_m, na.rm = TRUE) > 10) colorMax <- 100
      else colorMax <- as.numeric(stats::quantile(values$pct_m, probs = 0.999, na.rm = TRUE))
    }
    
    p[[i]] <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = values, ggplot2::aes(x = middle, y = group, fill = pct_m), width = width) +
      ggplot2::geom_rect(data = promoters, fill = "pink", ggplot2::aes(xmin = promoter_start, xmax = promoter_end, ymin = ymin, ymax = ymax)) +
      ggplot2::geom_rect(data = exons, fill = "black", ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax)) +
      ggplot2::geom_segment(data = genebody,
                            ggplot2::aes(x = x, xend = xend, y = -trackHeight, yend = -trackHeight),
                            arrow = ggplot2::arrow(length = grid::unit(trackHeight * arrowHeight, "cm"))) +
      ggplot2::geom_text(data = genenames,
                         ggplot2::aes(x = text_start, y = -trackHeight, label = gene_name),
                         hjust = 0) +
      ggplot2::scale_fill_gradientn(colors = pal, limits = c(0, colorMax), oob = scales::squish) +
      ggplot2::ylab(NULL) +
      ggplot2::xlab(paste0(genes[i], " (", chrom, ")")) +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        legend.position = if (legend) "right" else "none"
      )
  }
  
  gridExtra::grid.arrange(grobs = p, nrow = nrow)
}
