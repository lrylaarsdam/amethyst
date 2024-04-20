############################################################################################################################
#' @title Find cluster markers
#' @description Find variably methylated features across clusters
#'
#' @param obj The amethyst object to analyze. Clustering should be done prior to this step.
#' @param matrix Name of the pre-computed methylation values stored genomeMatrices slot for which to calculate variably methylated features
#' @returns A data.frame containing the p.value of methylation levels over each gene in each cluster relative to all other clusters
#' @export
#' @examples clusterMarkers <- findClusterMarkers(obj = obj, matrix = "gene_ch")
#' @importFrom dplyr select filter mutate group_by
#' @importFrom stats wilcox.test p.adjust
findClusterMarkers <- function(
    obj,
    matrix,
    genes) {

  genematrix <- as.matrix(obj@genomeMatrices[[matrix]])
  membership <- obj@metadata |> dplyr::select("cluster_id")

  results <- list()
  for (gene in genes) {
    results[[gene]] <- list()
    for (id in unique(membership$cluster_id)) {
      members <- rownames(membership |> dplyr::filter(cluster_id == id))
      nonmembers <- rownames(membership |> dplyr::filter(cluster_id != id))
      tryCatch(
        {
          results[[gene]][[id]] <- data.frame(
            "p.val" = stats::wilcox.test(
              x = genematrix[gene, members],
              y = genematrix[gene, nonmembers])$p.value,
            "gene" = gene,
            "cluster_id" = id,
            mean_1 = mean(genematrix[gene, members]),
            mean_2 = mean(genematrix[gene, nonmembers])
          ) |>
            dplyr::mutate(
              logFC = log2(mean_2 / mean_1),
              direction = ifelse(mean_1 > mean_2, "hypermethylated", "hypomethylated"))
        },
        error = function(e) {
          cat("Error processing gene:", gene, "and cluster:", id, "\n")
          results[[gene]][[id]] <- NA
        }
      )
    }
  }
  results <- do.call(rbind, lapply(results, function(x) do.call(rbind, x)))
  results <- results |> dplyr::group_by(cluster_id) |> dplyr::mutate(p.adj = stats::p.adjust(p.val, method = "bonferroni"))
  results <- results |> dplyr::select(p.val, p.adj, everything())
}
