############################################################################################################################
#' @title makeDoubletObject
#' @description This function makes the new object with true cells as well as artificial doublets
#'
#' @param obj An amethyst object containing the original cell data.
#' @param simFraction A numeric value specifying the fraction of cells to simulate as doublets. Default is 0.25.
#' @param threads An integer specifying the number of threads to use for parallel processing. Default is 10.
#' @param genomeMatrices A character vector specifying the names of the genome matrices to use. Default is c('cg_100k_score','ch_100k_pct').
#' @return An amethyst object containing both true cells and artificial doublets.
#' @export
#' @importFrom dplyr left_join mutate
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom caret createDataPartition train
#' @importFrom future multicore plan sequential
#' @importFrom future.apply future_lapply
#' @examples dbobj <- makeDoubletObject(obj, simFraction=0.25, threads = 10, genomeMatrices=c("cg_100kb_score", "ch_100k_pct"))
makeDoubletObject <- function(
    obj,
    simFraction = 0.25,
    threads = 10,
    genomeMatrices = c('cg_100k_score','ch_100k_pct')) {

  # Ensure the specified genome matrices exist
  for (i in 1:length(genomeMatrices)) {
    if (!any(names(slot(obj, "genomeMatrices")) == genomeMatrices[i])) {
      stop(paste("The genome matrix", sQuote(genomeMatrices[i]), "does not exist within 'genomeMatrices'."))
    }
  }

  # Set up parallel processing with a specified number of threads
  future::plan(future::multicore, workers = threads)  # Configures the number of threads to use

  # Number of artificial doublets to create
  cell_names <- colnames(as.matrix(obj@genomeMatrices[[genomeMatrices[[1]]]]))
  num_cells <- length(cell_names)
  num_doublets <- ceiling(num_cells * simFraction)

  # Set up the doublet pairings
  pairings <- matrix(sample(cell_names, num_doublets * 2, replace = TRUE), nrow = 2, ncol = num_doublets)

  # Create the new amethyst object
  dbobj <- amethyst::createObject()

  # Add in doublet status info
  cat("Adding doublet status info\n")
  doublet_info <- data.frame(cell_id = c(cell_names, paste0("Doublet_", seq_len(num_doublets))),
                             doublet_info = c(rep("True_Cell", num_cells), rep("Artif_Doub", num_doublets)),
                             stringsAsFactors = FALSE)
  rownames(doublet_info) <- doublet_info$cell_id

  dbobj@metadata <- doublet_info

  # Add in additional relevant cellInfo metadata, doublets will have NA values
  cat("Adding additional metadata\n")
  if ("cov" %in% colnames(obj@metadata)) {
    additional_cov <- data.frame(cell_id = rownames(obj@metadata), cov = obj@metadata[['cov']], stringsAsFactors = FALSE)
    rownames(additional_cov) <- additional_cov$cell_id
    dbobj@metadata <- dplyr::left_join(dbobj@metadata, additional_cov, by = "cell_id")
  }

  if ("mCG" %in% colnames(obj@metadata)) {
    additional_mcg <- data.frame(cell_id = rownames(obj@metadata), mCG = obj@metadata[['mCG']], stringsAsFactors = FALSE)
    rownames(additional_mcg) <- additional_mcg$cell_id
    dbobj@metadata <- dplyr::left_join(dbobj@metadata, additional_mcg, by = "cell_id")
  }

  if ("mCH" %in% colnames(obj@metadata)) {
    additional_mch <- data.frame(cell_id = rownames(obj@metadata), mCH = obj@metadata[['mCH']], stringsAsFactors = FALSE)
    rownames(additional_mch) <- additional_mch$cell_id
    dbobj@metadata <- dplyr::left_join(dbobj@metadata, additional_mch, by = "cell_id")
  }

  rownames(dbobj@metadata) <- dbobj@metadata$cell_id
  dbobj@metadata$cell_id <- NULL

  # Process each genome matrix
  newMatrices <- lapply(seq_along(genomeMatrices), function(i) {
    matrix_name <- genomeMatrices[i]
    pct_matrix <- as.matrix(obj@genomeMatrices[[matrix_name]])

    cat("Processing matrix:", matrix_name, "\n")

    # Generate the doublet_matrix directly using lapply for clarity
    doublet_matrix <- matrix(nrow = nrow(pct_matrix), ncol = num_doublets)

    # Use future_lapply for parallel processing
    doublets <- future.apply::future_lapply(1:ncol(pairings), function(col_idx) {
      cell1 <- pairings[1, col_idx]
      cell2 <- pairings[2, col_idx]
      sapply(seq_len(nrow(pct_matrix)), function(j) {
        vals <- c(pct_matrix[j, cell1], pct_matrix[j, cell2])
        if (all(is.na(vals))) {
          NA
        } else if (any(is.na(vals))) {
          # 50% chance to keep as NA or report the non-NA
          if (runif(1) < 0.5) {
            mean(vals, na.rm = TRUE)
          } else {
            NA
          }
        } else {
          mean(vals, na.rm = TRUE)
        }
      })
    }, future.seed = TRUE)

    # Bind the results back into a matrix
    doublet_matrix <- do.call(cbind, doublets)

    # Prepare combined matrix
    combined_matrix <- as.data.frame(cbind(pct_matrix, doublet_matrix))

    cat("Completed matrix:", matrix_name, "\n")
    return(combined_matrix)
  })

  if (threads > 1) {
    future::plan(future::sequential)
    gc()
  }

  # Name the slots in dbobj@genomeMatrices as the original slot
  for (i in seq_along(newMatrices)) {
    matrix_name <- genomeMatrices[i]
    cat("Adding matrix to new object:", matrix_name, "\n")

    dbobj@genomeMatrices[[matrix_name]] <- newMatrices[[i]]
	colnames(dbobj@genomeMatrices[[matrix_name]]) <- rownames(dbobj@metadata)

  }

  cat("Completed creating doublet object\n")
  return(dbobj)
}

############################################################################################################################
#' @title addMatrixMeans
#' @description This function calculates the mean values for each cell (column) in the specified genome matrix, excluding NAs,
#' and adds these mean values as a new column in the metadata of the provided amethyst object.
#' If you want to process irlba where mean is used for the NA values, this function can be used.
#'
#' @param obj An amethyst object containing the original cell data.
#' @param matrix_name A character string specifying the name of the genome matrix to process.
#' @param name A character string specifying the name of the new column to be added to the metadata.
#' @return An amethyst object with updated metadata including the mean values from the specified genome matrix.
#' @export
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column column_to_rownames
#' @examples
addMatrixMeans <- function(obj, matrix_name, name) {
  # Calculate the mean for each cell (column), excluding NAs
  cell_means <- apply(obj@genomeMatrices[[matrix_name]], 2, function(x) mean(x, na.rm = TRUE))

  # Create a data frame for the new metadata slot
  matrix_means <- data.frame(cell_id = colnames(obj@genomeMatrices[[matrix_name]]), stringsAsFactors = FALSE)
  matrix_means[[name]] <- cell_means
  rownames(matrix_means) <- matrix_means$cell_id
  matrix_means$cell_id <- NULL

  # Add the new metadata to the obj
  obj@metadata <- dplyr::left_join(
    obj@metadata |> tibble::rownames_to_column(var = "cell_id"),
    matrix_means |> tibble::rownames_to_column(var = "cell_id"),
    by = "cell_id"
  ) |> tibble::column_to_rownames(var = "cell_id")

  return(obj)
}

############################################################################################################################
#' @title buildDoubletModel
#' @description This function builds a model to identify doublets using a specified machine learning method.
#' The function splits the data into training and test sets, trains the model, and predicts probabilities on the test set.
#'
#' @param dbobj An amethyst object containing the data with doublets.
#' @param method A character string specifying the machine learning method to use for model training. Default is 'rf' (random forest).
#' @param reduction A character string specifying the name of the reduction slot to use. Default is 'irlba'.
#' @return A list containing the trained model, the predictions on the test set, and the test data.
#' @importFrom caret createDataPartition train
#' @importFrom stats predict
#' @importFrom randomForest randomForest
#' @export
#' @examples result <- buildDoubletModel(dbobj, reduction = "irlba", method="rf")
buildDoubletModel <- function(dbobj, method = 'rf', reduction = 'irlba') {
  # Ensure the irlba slot exists
  if (is.null(dbobj@reductions[[reduction]])) {
    stop("The reduction slot does not exist in the provided dbobj.")
  }

  # Ensure the doublet_info metadata exists
  if (!"doublet_info" %in% colnames(dbobj@metadata)) {
    stop("The 'doublet_info' metadata slot does not exist in the provided dbobj.")
  }

  # Extract the reduced dimension matrix from the irlba slot
  irlba_matrix <- dbobj@reductions[[reduction]]

  # Prepare the dataset
  labels <- dbobj@metadata$doublet_info
  data <- as.data.frame(irlba_matrix)
  data$label <- as.factor(labels)

  # Split into training and test sets
  set.seed(1)
  trainIndex <- caret::createDataPartition(data$label, p = 0.8, list = FALSE)
  trainData <- data[trainIndex, ]
  testData <- data[-trainIndex, ]

  # Train a logistic regression model
  model <- caret::train(label ~ ., data = trainData, method = method, family = "binomial")

  # Predict probabilities on the test set
  predictions <- stats::predict(model, newdata = testData, type = "prob")

  # Return the trained model and the predictions
  list(model = model, predictions = predictions, testData = testData)
}


############################################################################################################################
# this acores the cells and makes the doublet_score metadata slot
#' @title predictDoubletScores
#' @description This function calculates doublet scores for cells and adds them to the metadata of the provided doublet object.
#'
#' @param dbobj An amethyst doublet object containing the data with artificial doublets; output of makeDoubletObject
#' @param model A trained model used to predict doublet probabilities; output of buildDoubletMode
#' @param reduction A character string specifying the name of the reduction slot to use.
#' @return The amethyst doublet object with updated metadata including doublet scores.
#' @export
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom stats predict
#' @examples dbobj <- predictDoubletScores(dbobj, reduction = "irlba", model = result$model)
predictDoubletScores <- function(dbobj, model, reduction) {
  #Ensure the irlba slot exists
  if (is.null(dbobj@reductions[[reduction]])) {
    stop("The reduction slot does not exist in the provided dbobj.")
  }
  # Extract the reduced dimension matrix from the irlba slot
  irlba_matrix <- dbobj@reductions[[reduction]]

  # Predict probabilities on the entire dataset
  predictions <- predict(model, newdata = as.data.frame(irlba_matrix), type = "prob")

  # Store the doublet scores in the metadata slot
  dbobj@metadata <- dbobj@metadata %>%
    rownames_to_column(var = "cell_id") %>%
    mutate(doublet_score = predictions$`Artif_Doub`) %>%
    column_to_rownames(var = "cell_id")

  return(dbobj)
}

############################################################################################################################
#' @title addDoubletScores
#' @description This function adds doublet scores from a doublet object to the metadata of the original amethyst object.
#'
#' @param obj An amethyst object containing the original cell data.
#' @param dbobj An amethyst object containing the doublet scores in its metadata.
#' @return The original amethyst object with updated metadata including doublet scores.
#' @export
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column column_to_rownames
#' @examples obj <- addDoubletScores(obj, dbobj)
addDoubletScores <- function(obj, dbobj) {
  if (is.null(dbobj@metadata[['doublet_score']])) {
    stop("The 'doublet_score' metadata slot does not exist in the provided dbobj.")
  }

  # Prepare the doublet scores for joining
  doublet_score <- data.frame(
    cell_id = rownames(dbobj@metadata),
    doublet_score = dbobj@metadata[['doublet_score']],
    stringsAsFactors = FALSE
  )
  rownames(doublet_score) <- doublet_score$cell_id

  # Add cell_id to obj@metadata for joining
  obj_metadata <- obj@metadata |> tibble::rownames_to_column(var = "cell_id")

  # Join the doublet scores
  obj_metadata <- dplyr::left_join(obj_metadata, doublet_score, by = "cell_id")

  # Convert rownames back to original
  obj@metadata <- obj_metadata |> tibble::column_to_rownames(var = "cell_id")

  return(obj)
}



