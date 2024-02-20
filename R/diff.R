############################################################################################################################
#' findClusterMarkers
#' Find variably methylated features across clusters
#'
#' @param obj The amethyst object to analyze. Clustering should be done prior to this step.
#' @param impute Optional parameter to address coverage issues. If TRUE, Rmagic will be used to impute values based on available data.
#' @param genes Which genes to test. Genes must be indexed.
#' @param nmin Optional minimum threshold
#' @param type
#'
#' @return
#' @export
#'
#' @examples
findClusterMarkers <- function(
  obj,
  impute = FALSE,
  genes = names(obj@index[["CH"]]),
  nmin = NULL, # nmin is the minimum number of covered cytosines to include the cell in the analysis
  type) {

  markers <- makeWindows(obj, genes = genes, type = {{type}}, nmin = {{nmin}})
  if (impute) {
    markers [is.na(markers )] <- 0
    markers <- magic(markers, npca = 10)
    markers <- markers[["result"]]
  }
  markers <- merge(obj@metadata, markers, by = 0) %>% dplyr::rename("cell_id" = "Row.names")
  markers <- markers %>% mutate(cluster_id = paste0("cluster_", cluster_id))

  results <- list()
  for (gene in genes) {
    results[[gene]] <- list()
    for (id in unique(markers$cluster_id)) {
      tryCatch({
        results[[gene]][[id]] <-  wilcox.test(x = markers %>% dplyr::filter(cluster_id == id) %>% pull(paste(gene)),
                                              y = markers %>% dplyr::filter(cluster_id != id) %>% pull(paste(gene)))$p.value
        results[[gene]][[id]] <- results[[gene]][[id]] %>% dplyr::mutate(direction = ifelse(
          mean(markers %>% dplyr::filter(cluster_id == id) %>% pull(paste(gene))) >
            mean(markers %>% dplyr::filter(cluster_id != id) %>% pull(paste(gene))), "increase", "decrease"))
      }, error = function(e) {
        cat("Error processing gene:", gene, "and cluster:", id, "\n")
        results[[gene]][[id]] <- NA
      })
    }
  }
  results <- do.call(rbind, results)
  results <- pivot_longer(as.data.frame(results) %>% rownames_to_column("gene"), cols = starts_with("cluster"), names_to = "cluster", values_to = "p.val")
  results <- results %>% group_by(cluster) %>% mutate(p.adj = p.adjust(p.val, method = "bonferroni"))

}
