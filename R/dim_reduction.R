# Authors: Lauren Rylaarsdam, PhD; Andrew Adey, PhD
# 2024-2025
############################################################################################################################
#' @title dimEstimate
#' @description Estimate the nv value needed for singular value decomposition with irlba
#'
#' @param obj Amethyst object containing the matrix to calculate, which should be in the genomeMatrices slot
#' @param genomeMatrices Name of the matrix in the genomeMatrices slot to calculate nv
#' @param threshold Amount of variance that must be explained by the nv value, with 1 being 100% of variance
#' @param dims Number of singular values to test for each matrix
#'
#' @return Integer indicating the number of principal components required to meet the variance threshold for each matrix
#' @export
#'
#' @examples
#' \dontrun{
#'   dimEstimate(obj = combined, genomeMatrices = c("cg_100k_score", "ch_100k_pct"), dims = c(50, 50), threshold = 0.98)
#' }
#'
dimEstimate <- function(
    obj,
    genomeMatrices,
    dims,
    threshold = 0.98) {

  if (length(genomeMatrices) != length(dims)) {
    stop("Number of input matrices must equal the length of the dimension list")
  }

  svd_output <- list()
  for (i in 1:length(genomeMatrices)) {
    svd_output[[i]] <- obj@genomeMatrices[[genomeMatrices[i]]]
    windows <- rownames(svd_output[[i]])
    cells <- colnames(svd_output[[i]])
    svd_output[[i]][is.na(svd_output[[i]])] <- 0
    svd_output[[i]] <- irlba::irlba(as.matrix(svd_output[[i]]), dims[[i]])
    svd_output[[i]] <- as.data.frame(svd_output[[i]]$d)
    colnames(svd_output[[i]]) <- paste(genomeMatrices[[i]])
    rownames(svd_output[[i]]) <- paste("DIM", 1:dims[[i]])
  }
  svd_output <- do.call(cbind, svd_output)

  dims_to_use <- apply(svd_output, 2, function(x) {
    cumulative_variance_explained <- cumsum((x^2) / sum(x^2))
    if (any(cumulative_variance_explained >= threshold)) {
      min(which(cumulative_variance_explained >= threshold))
    } else {
      stop(paste("No nv value explains ", (threshold*100), "% of variance."))
    }
  })

  dims_to_use

}

############################################################################################################################
#' @title runIrlba
#' @description Perform dimensionality reduction based on methylation levels over a matrix stored in the @genomeMatrices slot
#'
#' @param obj Object for which to run irlba
#' @param genomeMatrices list of matrices in the genomeMatrices slot to use for irlba
#' @param replaceNA IRLBA can't accept NA values. Replace NA values with 0, 1, "mch_pct", or "mcg_pct"
#' @param dims list of how many dimensions to output for each matrix
#' @return Returns a matrix of appended irlba dimensions as columns and cells as rows
#' @importFrom irlba irlba
#' @export
#' @examples
#' \dontrun{
#'   obj@reductions[["irlba"]] <- runIrlba(obj, genomeMatrices = "cg_100k_score", dims = 10, replaceNA = 0)
#'   obj@reductions[["irlba"]] <- runIrlba(obj, genomeMatrices = c("ch_2M_pct", "ch_2M_score", "cg_2M_score"), dims = c(10, 10, 10))
#' }
runIrlba <- function(
    obj,
    genomeMatrices,
    dims,
    replaceNA = rep(0, length(dims))) {

  if (length(genomeMatrices) != length(dims)) {
    stop("Number of input matrices must equal the length of the dimension list")

  } else {

    matrix <- list()

    for (i in 1:length(genomeMatrices)) {
      matrix[[i]] <- obj@genomeMatrices[[genomeMatrices[i]]]
      default <- replaceNA[[i]]

      if (default == 0) {
        matrix[[i]][is.na(matrix[[i]])] <- 0
      } else if (default == 1) {
        matrix[[i]][is.na(matrix[[i]])] <- 1
      } else if (default %in% c("mch_pct", "mcg_pct")) {
        matrix[[i]] <- matrix[[i]][, colnames(matrix[[i]]) %in% rownames(obj@metadata)]
        if (default == "mch_pct") {
          replacement_vector <- obj@metadata$mch_pct[match(colnames(matrix[[i]]), rownames(obj@metadata))]
        }
        if (default == "mcg_pct") {
          replacement_vector <- obj@metadata$mcg_pct[match(colnames(matrix[[i]]), rownames(obj@metadata))]
        }
        # Replace NA values in the matrix with corresponding mch_pct values
        na_indices <- which(is.na(matrix[[i]]), arr.ind = TRUE)
        matrix[[i]][na_indices] <- replacement_vector[na_indices[, 2]]
      } else {
        stop("Default NA replacement must be 0, 1, 'mch_pct', or 'mcg_pct'.")
      }

      windows <- rownames(matrix[[i]])
      cells <- colnames(matrix[[i]])

      matrix[[i]] <- irlba::irlba(as.matrix(matrix[[i]]), dims[[i]])
      matrix[[i]] <- as.data.frame(matrix[[i]]$v)
      colnames(matrix[[i]]) <- paste("DIM", 1:dims[[i]], "_", genomeMatrices[i], sep = "")
      rownames(matrix[[i]]) <- cells
    }

    result <- do.call(cbind, matrix)
    return(result)
  }
}

############################################################################################################################
#' @title runTsne
#' @description This function runs t-Distributed Stochastic Neighbor Embedding (t-SNE) on single-cell methylation data contained within an amethyst object.
#' @param obj An amethyst object. A dimensionality reduction method, such as runIrlba, must have been performed
#' @param perplexity Numeric; the perplexity parameter for t-SNE. Typical values range from 5 to 50. Default is 30
#' @param method Character; the distance metric to use. Common choices are "euclidean", "cosine", etc. Default is "euclidean"
#' @param theta Numeric; speed/accuracy trade-off parameter for t-SNE. Values range from 0.0 (exact) to 1.0 (fast). Default is 0.5
#' @param reduction Character; the name of the dimensionality reduction to use from the object. Default is "irlba"
#' @return t-SNE coordinates labeled "dim_x" and "dim_y"
#' @export
#' @importFrom Rtsne Rtsne
#' @importFrom dplyr rename select
#' @importFrom tibble column_to_rownames
#' @examples
#' \dontrun{
#'   obj@reductions[["tsne"]] <- runTsne(obj, perplexity = 5, method = "euclidean", theta = 0.2, reduction = "irlba")
#' }
runTsne <- function(obj,
                    perplexity = 30,
                    method = "euclidean",
                    theta = 0.5,
                    reduction = "irlba") {

  # Ensure data is in a suitable format (data frame)
  cells <- rownames(obj@reductions[[reduction]])
  data_matrix <- as.matrix(obj@reductions[[reduction]])

  # Execute t-SNE
  tsne_result <- Rtsne::Rtsne(data_matrix,
                              dims = 2,
                              perplexity = {{perplexity}},
                              theta = {{theta}},
                              check_duplicates = FALSE,
                              pca = FALSE)

  # Create a data frame from t-SNE output and align row names
  tsne_dims <- as.data.frame(tsne_result$Y)
  rownames(tsne_dims) <- cells
  colnames(tsne_dims) <- c("dim_x", "dim_y")

  return(tsne_dims)
}

############################################################################################################################
#' @title runUmap
#' @description Perform dimension reduction with Uniform Manifold APproximation and Projection for Dimension Reduction
#'
#' @param obj Amethyst object for which to determine umap coordinates
#' @param neighbors Number of closest points to factor into projection calculations. A higher number will capture more global structure
#' @param dist Distance between point pairs
#' @param reduction Character; the name of the dimensionality reduction to use from the object. Default is "irlba"
#' @param method Distance metric to utilize. Default is euclidean
#'
#' @return A data frame of UMAP x and y coordinates labeled "dim_x" and "dim_y"
#' @importFrom umap umap
#' @export
#' @examples
#' \dontrun{
#'   obj@reductions[["umap"]] <- runUmap(obj, neighbors = 30, dist = 0.1, method = "euclidean", reduction = "irlba")
#' }
runUmap <- function(obj,
                    neighbors = 30,
                    dist = 0.1,
                    method = "euclidean",
                    reduction = "irlba") {

  umap_dims <- as.data.frame(umap::umap(as.data.frame(obj@reductions[[reduction]]),
                                        method = "naive",
                                        dims = 2,
                                        n_components = 2,
                                        n_neighbors = neighbors,
                                        min_dist = dist,
                                        metric = method)$layout)

  colnames(umap_dims) <- c("dim_x", "dim_y")
  return(umap_dims)

}
