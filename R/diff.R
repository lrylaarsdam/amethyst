#' @title Find cluster markers
#' @description Find variably methylated features across clusters
#'
#' @param obj The amethyst object to analyze. Clustering should be done prior to this step.
#' @param impute Optional parameter to address coverage issues. If TRUE, Rmagic will be used to impute values based on available data.
#' @param genes Which genes to test. Genes must be indexed.
#' @param nmin Optional minimum threshold
#' @param type
#'
#' @returns A data.frame of ...
#' @export
#'
#' @examples
#' @importFrom dplyr mutate rename filter everything select group_by
#' @importFrom Rmagic magic
findClusterMarkers <- function(
    obj,
    impute = FALSE,
    genes = names(obj@index[["CH"]]),
    nmin = NULL, # nmin is the minimum number of covered cytosines to include the cell in the analysis
    type = NULL,
    matrix = NULL) {
  if (is.null(matrix) && is.null(type)) {
    stop("If matrix is not provided, 'type' must be specified.")
  }

  if (!is.null(matrix)) {
    markers <- obj@genomeMatrices[[matrix]]
  } else if (is.null(matrix)) {
    markers <- makeWindows(obj, genes = genes, type = {{ type }}, nmin = {{ nmin }})
    if (impute) {
      markers[is.na(markers)] <- 0
      markers <- Rmagic::magic(markers, npca = 10)
      markers <- markers[["result"]]
    }
  }

  markers <- merge(obj@metadata, markers, by = 0) |>
    dplyr::rename("cell_id" = "Row.names")
  markers <- markers |>
    dplyr::mutate(cluster_id = paste0("cluster_", cluster_id))


  results <- list()
  for (gene in genes) {
    results[[gene]] <- list()
    for (id in unique(markers$cluster_id)) {
      cluster_markers <- markers %>%
        dplyr::filter(cluster_id == id)[[gene]]
      noncluster_markers <- markers %>%
        dplyr::filter(cluster_id != id)[[gene]]
      tryCatch(
        {
          results[[gene]][[id]] <- data.frame(
            "p.val" = wilcox.test(
              x = cluster_markers,
              y = noncluster_markers
            )$p.value,
            "gene" = gene,
            "id" = id,
            mean_1 = mean(cluster_markers),
            mean_2 = mean(noncluster_markers)
          ) |>
            dplyr::mutate(
              logFC = log2(mean_2 / mean_1),
              direction = ifelse(mean_1 > mean_2, "increase", "decrease")
            )
          # colnames(results[[gene]][[id]]) <- "p.val"
          # results[[gene]][[id]] <- results[[gene]][[id]] |>
          #   dplyr::mutate(
          #     gene = paste(gene),
          #     id = paste(id),
          #   mean_1 = mean(markers %>% dplyr::filter(cluster_id == id) %>% pull(paste(gene))),
          #   mean_2 = mean(markers %>% dplyr::filter(cluster_id != id) %>% pull(paste(gene))),
          #   logFC = log2(mean_2/mean_1),
          #   direction = ifelse(mean_1 > mean_2, "increase", "decrease"))
        },
        error = function(e) {
          cat("Error processing gene:", gene, "and cluster:", id, "\n")
          results[[gene]][[id]] <- NA
        }
      )
    }
  }
  results <- do.call(rbind, lapply(results, function(x) do.call(rbind, x)))
  results <- results %>% dplyr::select(-p.val, dplyr::everything(), p.val)
  results <- results %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(p.adj = p.adjust(p.val, method = "bonferroni"))
}
