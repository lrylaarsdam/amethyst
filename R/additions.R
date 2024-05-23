############################################################################################################################
### Nearest Neighbor Label Xfer
#' @title transferLabelsNN
#' @description Transfer metadata labels based on nearest neighbor analysis.
#'
#' @param obj Amethyst object on which to perform label transfer. The object must contain metadata with the specified labels and an `irlba` slot.
#' @param labels A character string specifying the metadata column in the Amethyst object that contains the labels to be transferred. Defaults to 'labels'.
#' @param transferred A character string specifying the name of the new column in the metadata slot where the transferred labels will be stored. Defaults to 'transferred'.
#' @param nnn An integer specifying the maximum number of nearest neighbors to search. Defaults to 5.
#' @return The Amethyst object with a new column in the metadata containing the transferred labels.
#' @export
#' @importFrom FNN get.knnx
transferLabelsNN <- function(obj,
                             labels = 'labels',
                             transferred = 'transferred',
                             nnn = 5,
                             reduction = "irlba") {
  # Check if the label slot exists in the metadata
  if (!labels %in% names(obj@metadata)) {
    stop(paste("The label slot", labels, "does not exist in obj@metadata."))
  }

  # Check if the reduction slot exists
  if (!reduction %in% names(obj@reductions)) {
    stop(paste0("The reduction ", reduction, " does not exist in the object."))
  }

  labels_vector <- obj@metadata[[labels]]
  dimred_matrix <- obj@reductions[[reduction]]

  # Initialize transferred labels with original labels
  transferred_labels <- labels_vector

  # Exclude rows with any NA values in dimred_matrix and corresponding labels
  valid_indices <- which(rowSums(is.na(dimred_matrix)) == 0)
  dimred_matrix_valid <- dimred_matrix[valid_indices, ]
  labels_vector_valid <- labels_vector[valid_indices]

  # Update indices for NA and non-NA labels after filtering
  na_indices_valid <- which(is.na(labels_vector_valid))
  non_na_indices_valid <- which(!is.na(labels_vector_valid))

  if (length(na_indices_valid) > 0 && length(non_na_indices_valid) > 0) {
    # Calculate nearest neighbors only for NA labeled rows using valid non-NA rows
    knn <- FNN::get.knnx(data = dimred_matrix_valid[non_na_indices_valid, , drop = FALSE],
                    query = dimred_matrix_valid[na_indices_valid, , drop = FALSE],
                    k = nnn)

    # Propagate labels based on mode of nearest neighbors
    for (i in seq_along(na_indices_valid)) {
      nearest_labels <- labels_vector_valid[non_na_indices_valid][knn$nn.index[i, ]]
      most_common_label <- names(which.max(table(nearest_labels)))
      transferred_labels[valid_indices[na_indices_valid[i]]] <- most_common_label
    }
  }

  # Store the transferred labels in the specified new slot
  obj@metadata[[transferred]] <- transferred_labels

  return(obj)
}

############################################################################################################################
### Cluster comp label xfer
#' @title clusterLabelTransfer
#' @description Transfers cluster labels from one set of annotations to another within an Amethyst object.
#' This function calculates the most common label within each cluster and assigns it as the representative label for that cluster.
#' @param obj An Amethyst object to perform cluster label transfer on. The object must contain metadata with the specified labels and cluster IDs.
#' @param labels A character string specifying the metadata column in the Amethyst object that contains the labels to be transferred. Defaults to 'type'.
#' @param cluster_ids A character string specifying the metadata column in the Amethyst object that contains the cluster IDs. Defaults to 'cluster_id'.
#' @return The Amethyst object with updated metadata containing the transferred cluster labels.
#' @export
#' @importFrom stats setNames ave
#' @examples
clusterLabelTransfer <- function(obj,
                                 labels = 'type',
                                 cluster_ids = 'cluster_id') {
  # Check if the necessary slots exist
  if (!labels %in% names(obj@metadata)) {
    stop(paste("The label slot", labels, "does not exist in obj@metadata."))
  }
  if (!cluster_ids %in% names(obj@metadata)) {
    stop(paste("The cluster_id slot", cluster_ids, "does not exist in obj@metadata."))
  }

  # Extract labels and cluster IDs
  label_vector <- obj@metadata[[labels]]
  cluster_vector <- obj@metadata[[cluster_ids]]

  # Calculate the frequency of each label within each cluster
  table_df <- as.data.frame(table(cluster_vector, label_vector))
  names(table_df) <- c("ClusterID", "Label", "Freq")

  # Calculate totals for each cluster to determine proportions
  total_per_cluster <- aggregate(Freq ~ ClusterID, data = table_df, FUN = sum)
  names(total_per_cluster)[2] <- "Total"

  # Merge totals back to main table and calculate proportions
  table_df <- merge(table_df, total_per_cluster, by = "ClusterID")
  table_df$Proportion <- with(table_df, Freq / Total)

  # Determine the label with the maximum proportion for each cluster
  idx_max_prop <- table_df$Proportion == stats::ave(table_df$Proportion, table_df$ClusterID, FUN = max)
  max_labels <- table_df[idx_max_prop, ]

  # If ties exist, select the first occurrence (This could be adapted if ties should be handled differently)
  max_labels <- max_labels[!duplicated(max_labels$ClusterID), ]

  # Map cluster IDs to their most common labels
  cluster_map <- stats::setNames(as.character(max_labels$Label), max_labels$ClusterID)

  # Create unique labels if necessary to resolve duplicates
  unique_labels <- make.unique(as.vector(cluster_map))
  cluster_map <- stats::setNames(unique_labels, names(cluster_map))

  # Assign this new mapping to a new slot
  obj@metadata[[paste0(cluster_ids, "_label_transfer")]] <- cluster_map[as.character(cluster_vector)]

  return(obj)
}
