#' @title Sample comp
#' @description Visualize the distribution of a metadata variable across groups of cells
#'
#' @param obj The amethyst object to plot
#' @param groupBy What column in the metadata to group cells by. This will be on the x axis.
#' @param colorBy What column in the metadata to color the bar chart by.
#'
#' @return A bar graph showing the distribution of categorical metadata in user-defined groups of cells
#' @export
#'
#' @examples sampleComp(obj, groupBy = sample, colorBy = cluster_id)
sampleComp <- function(
    obj,
    groupBy,
    colorBy,
    colors = NULL) {
  # Use ggplot2 to generate bar plot with the groupBy variable on the x axis and the fill determined by colorBy.
  ggplot2::ggplot(obj@metadata, aes(x = {{ groupBy }})) +
    ggplot2::geom_bar(aes(fill = {{ colorBy }}), position = "fill") +
    {
      if (!is.null(colors)) ggplot2::scale_fill_manual(values = {{ colors }})
    } +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "sample identity", y = "percentage", title = "sample composition") +
    ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}

#' @title markerHeatmap
#' @description Plot the methylation levels of sample groups across many genes at once
#'
#' @param obj The amethyst object to plot
#' @param level Hierarchy of reference to plot. i.e., "cell_class" plots "Exc"  "Inh"  "NonN"
#' @param impute Logical indicating whether to impute values using Rmagic
#' @param rownames Logical indicating whether to show row names on the heatmap
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
#'
#' @return Returns a heatmap of average methylation levels of key marker genes
#' @export
#'
#' @examples
markerHeatmap <- function(
    obj,
    level = "CellClass",
    impute = TRUE,
    rownames = FALSE,
    type = "CH") {
  markers <- keygenestable %>%
    dplyr::filter(cluster_level == level & marker %in% names(obj@index[[type]])) %>%
    dplyr::group_by(cluster_level) %>%
    dplyr::arrange(marker, .by_group = TRUE) %>%
    dplyr::pull(marker)
  markers <- makeWindows(obj, genes = markers, type = {{ type }})
  markers[is.na(markers)] <- 0
  if (impute) {
    markers <- Rmagic::magic(markers, npca = 10)
    markers <- markers[["result"]]
  }
  rel_key_genes_for_join <- keygenestable %>%
    dplyr::filter(cluster_level == level) %>%
    dplyr::rename("gene" = "marker")

  rel_meta_data <- obj@metadata |>
    tibble::rownames_to_column(var = "cell_id")

  markers <- markers |>
    tibble::rownames_to_column(var = "cell_id") |>
    tidyr::pivot_longer(
      cols = c(2:(length(colnames(markers)) + 1)),
      names_to = "gene",
      values_to = "pct_m"
    ) |>
    dplyr::inner_join(rel_key_genes, by = "gene") |>
    dplyr::inner_join(rel_meta_data, by = "cell_id")


  # markers <- dplyr::inner_join(markers, rel_key_genes, by = "gene")
  # markers <- dplyr::inner_join(markers, rel_meta_data, by = "cell_id")
  markers <- markers |>
    dplyr::rowwise() |>
    dplyr::mutate(group = paste0(make_clean_names(description), "_", gene)) |>
    dplyr::arrange(cluster_id, description)
  toplot <- tidyr::pivot_wider(
    dplyr::select(markers, c(cell_id, pct_m, group)),
    names_from = "group",
    values_from = "pct_m"
  ) |> tibble::column_to_rownames(var = "cell_id")

  # plot pheatmap
  annot_col <- markers |>
    dplyr::select(group, name) |>
    dplyr::distinct(group, name) |>
    tibble::column_to_rownames("group")
  annot_row <- markers |>
    dplyr::select(cell_id, cluster_id) |>
    dplyr::distinct(cell_id, cluster_id) |>
    tibble::column_to_rownames("cell_id")

  labels <- markers |>
    dplyr::select(group, name, gene) |>
    dplyr::distinct(group, name, .keep_all = TRUE) |>
    dplyr::pull("gene")

  pheatmap::pheatmap(toplot,
    annotation_row = annot_row,
    annotation_col = annot_col,
    show_rownames = rownames,
    labels_col = labels,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
}


#' @title Cluster compare
#' @description Correlates average percent methylation over 100kb windows aggregated by metadata to a reference matrix to aid in cell type identification.
#'
#' @param matrix Average percent methylation aggregated over 100kb windows by a grouping factor (see makeWindows). Column names should be genomic positions.
#' @param ref The reference structure to correlate to. The default for brain used here is from PMC8494641.
#' @param level The grouping level in the reference to compare to
#' @param type The type of methylation in the reference to compare to (cg or ch). Make sure the input matrix reflects the same type.
#' @param method Correlation coefficient to be computed ("pearson", "kendall", or "spearman")
#' @param n Number of top correlations to plot.
#'
#' @return A ggplot dot plot showing the top 5 correlations of each group to the reference
#' @export
#'
#' @examples clusterCompare(cluster_100kb_cg, level = "major_type", type = "cg")
clusterCompare <- function(
    matrix,
    ref = celltyperefs,
    level,
    type,
    method = "pearson",
    n = 5,
    aggregate = NULL) {
  matrix <- obj@genomeMatrices[[matrix]]

  # append the reference matrix to the experimental
  withref <- dplyr::bind_rows(matrix, ref[[level]][[type]])

  # generate correlation matrix
  cor <- stats::cor(t(withref), method = {{ method }}, use = "pairwise.complete.obs")

  # determine top correlations
  annotation <- tidyr::pivot_longer(as.data.frame(cor) %>% rownames_to_column("group1"), cols = c(2:(length(colnames(cor)) + 1)), names_to = "group2", values_to = "cor") %>%
    dplyr::group_by(group1) %>%
    dplyr::arrange(desc(cor), .by_group = TRUE) %>%
    dplyr::filter((group1 %in% rownames(matrix)) & !(group2 %in% rownames(matrix))) %>%
    top_n(n = {{ n }})

  p <- ggplot2::ggplot(annotation, ggplot2::aes(x = group1, y = group2, size = cor, color = cor)) +
    ggplot2::geom_point() +
    ggplot2::theme_classic() +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(x = "cluster", y = paste0("Ref ", type, " methylation by ", level)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p
}


#' @title noAxes
#' @description Add to a ggplot object to remove axes
#' Developed by the Satija lab https://github.com/satijalab/seurat/blob/master/R/visualization.R
noAxes <- function(..., keep.text = FALSE, keep.ticks = FALSE) {
  blank <- element_blank()
  no.axes.theme <- ggplot2::theme(
    # Remove the axis elements
    axis.line.x = blank,
    axis.line.y = blank,
    # Validate the theme
    validate = TRUE,
    ...
  )
  if (!keep.text) {
    no.axes.theme <- no.axes.theme +
      ggplot2::theme(
        axis.text.x = blank,
        axis.text.y = blank,
        axis.title.x = blank,
        axis.title.y = blank,
        validate = TRUE,
        ...
      )
  }
  if (!keep.ticks) {
    no.axes.theme <- no.axes.theme +
      ggplot2::theme(
        axis.ticks.x = blank,
        axis.ticks.y = blank,
        validate = TRUE,
        ...
      )
  }
  no.axes.theme
}

#' @title UMAP feature
#' Plot any feature in the metadata in UMAP space
#'
#' @param obj The amethyst object to plot
#' @param colorBy The metadata feature to plot in UMAP space
#' @param colors Optional color scale
#' @param pointSize Optional adjustment of point size
#'
#' @return Returns a ggplot graph plotting xy coordinates of cells colored according to the specified feature
#' @export
#'
#' @examples umapFeature(obj, colorBy = cluster_id)
#'
umapFeature <- function(
    obj,
    colorBy,
    colors = NULL,
    pointSize = 0.5) {
  p <- ggplot2::ggplot(obj@metadata, ggplot2::aes(x = umap_x, y = umap_y, color = {{ colorBy }})) +
    ggplot2::geom_point(size = pointSize) +
    {
      if (!is.null(colors)) ggplot2::scale_color_manual(values = {{ colors }})
    } +
    ggplot2::theme_classic() +
    noAxes() +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3)))
  p
}


#' @title umapGeneM
#' @description Plot cumulative % methylation over a gene body on umap
#'
#' @param obj The amethyst object to plot
#' @param genes List of genes to plot % methylation in UMAP space
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
#' @param impute Logical indicating whether to impute values using Rmagic
#' @param blend Logical indicating whether to blend methylation levels of two or three genes
#' @param nrow Number of rows to distribute graphs over
#' @param matrix Pre-computed matrix to use in the genomeMatrices slot, if available. The type parameter is not needed if matrix is definied.
#'
#' @return Returns a plot of methylation over given genomic region in relation to UMAP
#' @export
#'
#' @examples
#' \dontrun{
#' umapGeneM(obj, genes = c("GAD1", "SATB2"), type = "CH")
#' umapGeneM(obj, genes = c("GAD1", "SATB2"), type = "CH", blend = TRUE, nrow = 1)
#' umapGeneM(both, genes = c("Cux1", "Lingo1"), matrix = "genemeans_cac_imputed", blend = T)
#' }
umapGeneM <- function(
    obj,
    genes,
    type = NULL,
    impute = FALSE,
    matrix = NULL,
    blend = FALSE,
    colors = NULL,
    nrow = round(sqrt(length(genes)))) {
  if (is.null(matrix) && is.null(type)) {
    stop("If matrix is not provided, 'type' must be specified.")
  }

  if (!blend) {
    # make empty plot list
    p <- vector("list", length(genes)) # empty plot list
    for (i in 1:length(genes)) {
      # get average methylation across a gene
      if (!is.null(matrix)) {
        genem <- obj@genomeMatrices[[matrix]]
        genem <- genem %>%
          dplyr::select(genes[i]) %>%
          tibble::rownames_to_column(var = "cell_id")
        colnames(genem) <- c("cell_id", "pct_m")
      } else if (is.null(matrix)) {
        genem <- getGeneM({{ obj }}, gene = genes[i], type = {{ type }})
        genem <- genem %>%
          dplyr::group_by(cell_id) %>%
          dplyr::summarise(n = n(), pct_m = (sum(c != 0) * 100 / (sum(c != 0) + sum(t != 0))), .groups = "keep") %>%
          dplyr::filter(n >= 2)
        if (impute) {
          imputed <- magic(genem$pct_m)
          genem$pct_m <- imputed[["result"]]$V1
        }
      }
      # merge with metadata
      plot <- dplyr::inner_join(genem, obj@metadata %>% rownames_to_column(var = "cell_id"), by = "cell_id")
      # plot
      p[[i]] <- ggplot2::ggplot(plot, ggplot2::aes(x = umap_x, y = umap_y, color = pct_m)) +
        ggplot2::geom_point(size = 0.5) +
        ggplot2::theme_classic() +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3))) +
        ggplot2::scale_color_viridis_c(name = "pct.genebody") +
        {
          if (!is.null(colors)) ggplot2::scale_color_gradientn(colors = {{ colors }})
        } +
        ggplot2::labs(title = paste0("%", type, " methylation across ", genes[i], " gene body")) +
        noAxes()
    }
    gridExtra::grid.arrange(grobs = p, nrow = nrow)
  } else if (blend) {
    if (!is.null(matrix)) {
      genem <- obj@genomeMatrices[[matrix]]
      genem <- genem %>% dplyr::select(genes)
    } else if (is.null(matrix)) {
      genem <- makeWindows(obj = {{ obj }}, genes = {{ genes }}, type = {{ type }})
      genem[is.na(genem)] <- 0
      if (impute) {
        genem <- Rmagic::magic(genem, npca = 10)
        genem <- genem[["result"]]
      }
    }
    genem <- genem * 0.01
    if (length(genes) == 2) {
      names(genem) <- c("gene1", "gene2")
      genem <- genem |>
        dplyr::mutate(mix = rgb(red = gene1, green = gene2, blue = 0, maxColorValue = max(genem[1:2])))
      plot <- merge(genem, obj@metadata, by = 0)
      legend <- expand.grid(red = seq(0, max(genem[1:2]), by = 0.02), green = seq(0, max(genem[1:2]), by = 0.02))
      legend <- within(legend, mix <- rgb(green = green, red = red, blue = 0, maxColorValue = max(genem[1:2])))
      p1 <- ggplot2::ggplot(plot, ggplot2::aes(x = umap_x, y = umap_y, color = mix)) +
        ggplot2::geom_point(size = 0.5) +
        ggplot2::theme_classic() +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3))) +
        ggplot2::scale_color_identity() +
        ggplot2::labs(title = paste0(type, " methylation across ", genes[1], " and ", genes[2])) +
        noAxes()
      p2 <- ggplot2::ggplot(legend, ggplot2::aes(x = red, y = green)) +
        ggplot2::geom_tile(aes(fill = mix), color = "white") +
        ggplot2::scale_fill_identity() +
        ggplot2::theme_classic() +
        ggplot2::labs(x = paste(genes[1]), y = paste(genes[2]), title = "Legend")
      p2 <- p2 + ggplot2::coord_fixed(ratio = 1 / 1)
      cowplot::plot_grid(p1, p2, rel_widths = c(3, 1))
    } else if (length(genes) == 3) {
      names(genem) <- c("gene1", "gene2", "gene3")
      genem <- genem %>% dplyr::mutate(mix = rgb(red = gene1, green = gene2, blue = gene3, maxColorValue = max(genem[1:3])))
      plot <- merge(genem, obj@metadata, by = 0)
      p1 <- ggplot(plot, aes(x = umap_x, y = umap_y, color = mix)) +
        geom_point(size = 0.5) +
        theme_classic() +
        guides(colour = guide_legend(override.aes = list(size = 3))) +
        scale_color_identity() +
        labs(title = paste0(type, " methylation across ", genes[1], ", ", genes[2], ", and ", genes[3])) +
        noAxes() +
        theme(legend.position = "none")
      legend <- expand.grid(red = seq(0, .1, by = 0.02), blue = seq(0, .1, by = 0.02), green = seq(0, .1, by = 0.02))
      legend <- within(legend, mix <- rgb(green = green, red = red, blue = blue, maxColorValue = 0.1))
      p2 <- plotly::plot_ly(x = legend$red, y = legend$green, z = legend$blue, color = ~ I(legend$mix)) %>%
        plotly::layout(scene = list(
          xaxis = list(title = paste(genes[[1]])),
          yaxis = list(title = paste(genes[[2]])),
          zaxis = list(title = paste(genes[[3]]))
        ))
      plotly::subplot(p1, p2, nrows = 1, widths = c(.8, 0.2))
    } else {
      print("Error: Blending only accommodates 2 or 3 genes")
    }
  }
}

############################################################################################################################
# plot cumulative % methylation over a gene body with violin plot

vlnGeneM <- function(
    obj,
    genes,
    type = NULL,
    matrix = NULL,
    groupBy,
    colors = NULL,
    nrow = round(sqrt(length(genes)))) {
  if (is.null(matrix) && is.null(type)) {
    stop("If matrix is not provided, 'type' must be specified.")
  }

  if (is.null(colors)) {
    pal <- c("#004A4A", "#F05252", "#419e9e", "#fcba2b", "#bd083e", "#FB9A99", "#75C8D2", "#FF8B73", "#B2DF8A", "#1F78B4", "#E31A1C", "#aae3e3", "#FFA976")
  } else if (!is.null(colors)) {
    pal <- colors
  }

  if (!is.null(matrix)) {
    genem <- obj@genomeMatrices[[matrix]]
    genem <- genem %>% dplyr::select(genes)
    genem <- merge(genem, obj@metadata, by = 0) %>% dplyr::rename("cell_id" = "Row.names")
    genem <- pivot_longer(genem, cols = genes, names_to = "gene", values_to = "pct_m")

    ggplot(genem, aes(x = {{ groupBy }}, y = pct_m, fill = {{ groupBy }}, color = {{ groupBy }})) +
      geom_violin(alpha = 0.3) +
      geom_jitter(aes(color = {{ groupBy }}), size = 0.2) +
      theme_classic() +
      scale_fill_manual(values = pal) +
      scale_color_manual(values = pal) +
      facet_wrap(vars(gene))
  } else if (is.null(matrix)) {
    p <- vector("list", length(genes)) # empty plot list
    for (i in 1:length(genes)) {
      genem <- getGeneM({{ obj }}, gene = genes[i], type = {{ type }})
      genem <- genem %>%
        dplyr::group_by(cell_id) %>%
        dplyr::summarise(n = n(), pct_m = (sum(c != 0) * 100 / (sum(c != 0) + sum(t != 0))), .groups = "keep") %>%
        dplyr::filter(n >= 5)
      genem <- inner_join(genem, obj@metadata %>% rownames_to_column(var = "cell_id"), by = "cell_id")
      p[[i]] <- ggplot(genem, aes(x = {{ groupBy }}, y = pct_m, fill = {{ groupBy }}, color = {{ groupBy }})) +
        geom_violin(alpha = 0.3) +
        geom_jitter(aes(color = {{ groupBy }}), size = 0.2) +
        theme_classic() +
        scale_fill_manual(values = pal) +
        scale_color_manual(values = pal) +
        ggtitle(paste0("Avg m", type, " across ", genes[i], " gene body")) +
        theme(legend.position = "none")
    }
    gridExtra::grid.arrange(grobs = p, nrow = nrow)
  }
}

############################################################################################################################

dotGeneM <- function(
    obj,
    genes,
    type,
    groupBy) {
  genem <- makeWindows({{ obj }}, gene = {{ genes }}, type = {{ type }}, nmin = 3)
  genem <- merge(genem, obj@metadata, by = 0) %>% column_to_rownames(var = "Row.names")
  genem <- pivot_longer(genem, cols = genes, names_to = "gene", values_to = "pct_m")
  genem <- genem %>%
    dplyr::group_by(gene, {{ groupBy }}) %>%
    dplyr::summarise(pct_nonzero = sum(pct_m != 0, na.rm = T) * 100 / (n()), n = n(), avg_m = mean(pct_m, na.rm = TRUE), .groups = "keep")

  ggplot(genem, aes(x = gene, y = {{ groupBy }})) +
    geom_point(aes(size = n, color = avg_m)) +
    theme_classic() +
    scale_color_viridis_c(name = paste0("avg_m", type)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}

############################################################################################################################
# plot pattern of % methylation over gene body

histGeneM <- function(
    obj,
    type, # "CG" or "CH" methylation, or "both"
    groupBy = NULL, # metadata to divide by
    genes,
    bins = 100, # number of columns to plot across gene body
    colors = pal,
    smooth = TRUE,
    track = TRUE,
    promoter = FALSE,
    track.color = "cornflowerblue",
    alpha = 0.2,
    nrow = length(genes),
    legend = F) {
  # make empty plot list
  p <- vector("list", length(genes)) # empty plot list

  if (is.null(obj@ref)) {
    print("No reference found. Make sure obj@ref slot contains the proper annotation file.")
  }
  for (i in 1:length(genes)) {
    # determine size of bins
    binwidth <- (obj@ref %>% dplyr::filter(type == "gene" & gene_type == "protein_coding") %>% dplyr::mutate(length = end - start) %>% dplyr::group_by(gene_name) %>% dplyr::mutate(replicate = n()) %>%
      dplyr::group_by(gene_name) %>% arrange(desc(length)) %>% dplyr::filter(row_number() == 1 & gene_name == genes[i]) %>% dplyr::pull(length)) / bins

    # get average methylation across a gene
    if (type == "CG" | type == "CH") {
      genem <- getGeneM({{ obj }}, gene = genes[i], type = {{ type }})
    } else if (type == "both") {
      both <- as.list(c("CH", "CG"))
      genem <- map(.x = both, .f = function(x) {
        getGeneM({{ obj }}, gene = genes[i], type = x)
      })
      names(genem) <- c("CH", "CG")
      genem <- rbindlist(genem, idcol = "type") %>% relocate(type, .after = last_col())
    }
    genem <- genem %>% dplyr::mutate(start = round_any(pos, binwidth, floor), end = round_any(pos, binwidth, ceiling), bin = paste0(chr, "_", round_any(pos, binwidth, floor), "_", round_any(pos, binwidth, ceiling)))
    if (type == "CG" | type == "CH") {
      genem <- genem %>%
        dplyr::group_by(cell_id, bin) %>%
        dplyr::summarise(n = n(), start = mean(start), end = mean(end), middle = mean(c(start, end)), pct_m = (sum(c != 0) * 100 / (sum(c != 0) + sum(t != 0))), .groups = "keep")
    } else if (type == "both") {
      genem <- genem %>%
        dplyr::group_by(cell_id, bin, type) %>%
        dplyr::summarise(n = n(), start = mean(start), end = mean(end), middle = mean(c(start, end)), pct_m = (sum(c != 0) * 100 / (sum(c != 0) + sum(t != 0))), .groups = "keep")
    }
    # merge with metadata for grouping information
    genem <- inner_join(genem, obj@metadata %>% rownames_to_column(var = "cell_id"), by = "cell_id")
    if (type == "CG" | type == "CH") {
      summary <- genem %>%
        dplyr::group_by(start, middle, end, bin, {{ groupBy }}) %>%
        dplyr::summarise(mean = mean(pct_m))
      summary$type <- as.character(type)
    } else if (type == "both") {
      summary <- genem %>%
        dplyr::group_by(start, middle, end, bin, {{ groupBy }}, type) %>%
        dplyr::summarise(mean = mean(pct_m))
    }
    if (promoter) {
      promoters <- getGeneM({{ obj }}, gene = genes[i], type = "promoters")
      promoters <- inner_join(promoters, obj@metadata %>% rownames_to_column(var = "cell_id"), by = "cell_id")
      promoters <- promoters %>%
        dplyr::group_by(cell_id, {{ groupBy }}) %>%
        dplyr::summarise(n = n(), pct_m = (sum(c != 0) * 100 / (sum(c != 0) + sum(t != 0))), type = "promoter", .groups = "keep")
      promoters$middle <- obj@ref %>%
        dplyr::filter(type == "gene" & gene_type == "protein_coding") %>%
        mutate(length = end - start) %>%
        group_by(gene_name) %>%
        dplyr::mutate(replicate = n()) %>%
        group_by(gene_name) %>%
        arrange(desc(length)) %>%
        dplyr::filter(row_number() == 1) %>%
        mutate(tss = ifelse(strand == "+", start, end)) %>%
        mutate(promoter_start = tss - 1500, promoter_end = tss + 1500) %>%
        distinct(seqid, start, end, .keep_all = T) %>%
        dplyr::filter(seqid != "chrM") %>%
        select("seqid", "promoter_start", "promoter_end", "gene_name", "location") %>%
        arrange(seqid, promoter_start) %>%
        dplyr::filter(gene_name %in% genes[i]) %>%
        dplyr::mutate(middle = mean(c(promoter_start, promoter_end))) %>%
        dplyr::pull(middle)
      promoters <- promoters %>%
        dplyr::mutate(start = (middle - 1500), end = (middle + 1500), bin = "promoter") %>%
        dplyr::group_by(start, middle, end, {{ groupBy }}, type) %>%
        dplyr::summarise(mean = mean(pct_m))
      summary <- do.call(rbind, list(summary, promoters))
    }
    if (track) {
      exon <- obj@ref %>%
        dplyr::filter(gene_name == genes[i] & type == "exon") %>%
        distinct(seqid, start, end)
    }
    # plot
    p[[i]] <- ggplot() +
      {
        if (track & promoter) geom_rect(data = promoters, aes(xmin = start, xmax = end, ymin = (0 - (max(summary$mean) * .06)), ymax = (0 - (max(summary$mean) * .03))), fill = "red", alpha = .4)
      } +
      {
        if (track) geom_rect(data = exon, aes(xmin = start, xmax = end, ymin = (0 - (max(summary$mean) * .06)), ymax = (0 - (max(summary$mean) * .03))), color = track.color, fill = track.color, alpha = 1)
      } + # if exons is true, plot exons track
      {
        if (track) geom_rect(data = exon, aes(xmin = min(start), xmax = max(end), ymin = (0 - (max(summary$mean) * .045)), ymax = (0 - (max(summary$mean) * .045))), color = track.color, alpha = 1)
      } + # also plot light line for intron track
      {
        if (type == "CG" | type == "CH") geom_rect(data = summary %>% dplyr::filter(type != "promoter"), aes(xmin = start, xmax = end, ymin = 0, ymax = mean, fill = {{ groupBy }}), alpha = 0.4)
      } +
      {
        if (type == "both") geom_rect(data = summary %>% dplyr::filter(type != "promoter"), aes(xmin = start, xmax = end, ymin = 0, ymax = mean, fill = {{ groupBy }}, alpha = type))
      } +
      {
        if (promoter) geom_rect(data = summary %>% dplyr::filter(type == "promoter"), aes(xmin = start, xmax = end, ymin = 0, ymax = mean), fill = "red", alpha = 0.2)
      } +
      {
        if (smooth & type == "CG") geom_smooth(data = genem, aes(x = middle, y = pct_m, color = {{ groupBy }}))
      } + # smoothed lines of mean % CH of exons and introns
      {
        if (smooth & type == "CH") geom_smooth(data = genem, aes(x = middle, y = pct_m, color = {{ groupBy }}))
      } +
      {
        if (smooth & type == "both") geom_smooth(data = genem, aes(x = middle, y = pct_m, color = {{ groupBy }}, alpha = type))
      } +
      facet_wrap(vars({{ groupBy }}), nrow = 1) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + # wrap by cell class; make graph look nice
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) + # colors
      {
        if (!legend) theme(legend.position = "none")
      } +
      scale_alpha_manual(values = c(.2, .6)) +
      {
        if (type != "both") ggtitle(paste0({{ type }}, " methylation along ", genes[i], " body"))
      } +
      {
        if (type == "both") ggtitle(paste0("CG+CH methylation along ", genes[i], " body"))
      } +
      labs(x = paste0("Location along ", unique(separate(data = genem, col = bin, into = "chr", sep = "_", extra = "drop")$chr)))
  }
  gridExtra::grid.arrange(grobs = p, nrow = nrow)
}
