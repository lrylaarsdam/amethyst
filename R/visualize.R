############################################################################################################################
#' @title sampleComp
#'
#' @description Visualize the distribution of a metadata variable across groups of cells
#' @param obj The amethyst object to plot
#' @param groupBy What column in the metadata to group cells by. This will be on the x axis.
#' @param colorBy What column in the metadata to color the bar chart by.
#' @return A bar graph showing the distribution of categorical metadata in user-defined groups of cells
#' @export
#' @examples sampleComp(obj = obj, groupBy = "sample", colorBy = "cluster_id")
#' @importFrom ggplot2 ggplot geom_bar scale_fill_manual theme_classic labs theme
sampleComp <- function(
    obj,
    groupBy,
    colorBy,
    colors = NULL) {
  # Use ggplot2 to generate bar plot with the groupBy variable on the x axis and the fill determined by colorBy.
  ggplot2::ggplot(obj@metadata, aes(x = .data[[groupBy]])) +
    ggplot2::geom_bar(aes(fill = .data[[colorBy]]), position = "fill") +
    {if (!is.null(colors)) ggplot2::scale_fill_manual(values = {{ colors }})} +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "sample identity", y = "percentage", title = "sample composition") +
    ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}

############################################################################################################################
#' @title noAxes
#' @description Add to a ggplot object to remove axes
#'
#' Developed by the Satija lab https://github.com/satijalab/seurat/blob/master/R/visualization.R
#' @export
#' @importFrom ggplot2 theme
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

############################################################################################################################
#' @title dimFeature
#'
#' @description Plot any feature in the metadata in UMAP or TSNE space
#'
#' @param obj The amethyst object to plot
#' @param colorBy The metadata feature to plot in UMAP or TSNE space
#' @param colors Optional color scale
#' @param reduction Whether to plot "umap" or "tsne" coordinates. Note runUmap and/or runTsne must have been run.
#' @param pointSize Optional adjustment of point size
#'
#' @return Returns a ggplot graph plotting xy coordinates of cells colored according to the specified feature
#' @export
#' @examples dimFeature(obj = obj, colorBy = "cluster_id", reduction = "umap")
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_classic guides guide_legend element_blank
dimFeature <- function(
    obj,
    colorBy,
    colors = NULL,
    pointSize = 0.5,
    reduction = "umap") {

  if (reduction == "umap") {
    dim_x <- "umap_x"
    dim_y <- "umap_y"
  }
  if (reduction == "tsne") {
    dim_x <- "tsne_x"
    dim_y <- "tsne_y"
  }

  # shuffle points so they aren't directly on top of each other
  set.seed(123)
  obj@metadata <- obj@metadata[sample(x = 1:nrow(x = obj@metadata)), ]

  p <- ggplot2::ggplot(obj@metadata, ggplot2::aes(x = .data[[dim_x]], y = .data[[dim_y]], color = {{colorBy}})) +
    ggplot2::geom_point(size = pointSize) +
    {if (!is.null(colors)) ggplot2::scale_color_manual(values = {{colors}})} +
    ggplot2::theme_classic() +
    noAxes() +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)))
  p

}

############################################################################################################################
#' @title dimM
#' @description Plot cumulative % methylation over a gene body on UMAP or TNSE dimensionality reductions
#'
#' @param obj The amethyst object to plot
#' @param genes List of genes to plot % methylation in UMAP or TSNE space
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
#' @param blend Logical indicating whether to blend methylation levels of two or three genes
#' @param nrow Number of rows to distribute graphs over
#' @param metric If matrix is not specified, what methylation parameter to calculate, i.e. "percent", "score", or "ratio"
#' @param colors Optional color gradient to include
#' @param matrix Pre-computed matrix to use in the genomeMatrices slot, if available. The type parameter is not needed if matrix is definied.
#' @return Returns a plot of methylation over given genomic region using UMAP or TSNE coordinates
#' @export
#' @examples dimM(obj, genes = c("GRM8", "SATB1"), matrix = "gene_ch", pointSize = 0.8)
#' @importFrom dplyr inner_join mutate
#' @importFrom ggplot2 ggplot aes geom_point theme_classic guides scale_color_viridis_c scale_color_gradientn labs scale_color_identity geom_tile scale_fill_identity theme
#' @importFrom grDevices rgb
#' @importFrom gridExtra grid.arrange
#' @importFrom plotly plot_ly subplot
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom scales squish
#' @importFrom cowplot plot_grid
dimM <-  function(
    obj,
    genes,
    matrix,
    blend = FALSE,
    colors = NULL,
    squish = NULL,
    reduction = "umap",
    pointSize = 0.1,
    nrow = round(sqrt(length(genes)))) {

  if (is.null(obj@genomeMatrices[[matrix]])) {
    stop("Please construct a matrix of methylation levels over the desired genes using makeWindows and store in the genomeMatrices slot.")
  }

  if (reduction == "umap") {
    dim_x <- "umap_x"
    dim_y <- "umap_y"
  }
  if (reduction == "tsne") {
    dim_x <- "tsne_x"
    dim_y <- "tsne_y"
  }

  if (!blend) {

    if (!is.null(colors)) {
      pal <- colors
    } else {
      pal <- c("black", "turquoise", "gold", "red")
    }

    # make empty plot list
    p <- vector("list", length(genes)) # empty plot list
    for (i in 1:length(genes)) {
      # get average methylation across a gene
      genem <- obj@genomeMatrices[[matrix]]
      genem <- genem[genes[i],] |> tibble::rownames_to_column(var = "gene")
      genem <- tidyr::pivot_longer(genem, cols = c(2:ncol(genem)), names_to = "cell_id", values_to = "pctm")
      genem <- genem |> dplyr::filter(!is.na(pctm))

      # merge with metadata
      plot <- dplyr::inner_join(genem, obj@metadata |> tibble::rownames_to_column(var = "cell_id"), by = "cell_id")
      # plot
      p[[i]] <- ggplot2::ggplot(plot, ggplot2::aes(x = .data[[dim_x]], y = .data[[dim_y]], color = pctm)) +
        ggplot2::geom_point(size = pointSize) + ggplot2::theme_classic() +
        ggplot2::scale_color_gradientn(colors = pal) +
        {if (!is.null(squish)) ggplot2::scale_color_gradientn(colors = pal, limits = c(0, squish), oob = scales::squish)} +
        ggplot2::labs(title = paste0(genes[i])) + noAxes()
    }
    gridExtra::grid.arrange(grobs = p, nrow = nrow)
  }

  if (blend) {
    genem <- obj@genomeMatrices[[matrix]]
    genem <- genem[genes,]
    genem <- as.data.frame(t(genem))
    genem <- genem * 0.01
    if (length(genes) == 2) {
      names(genem) <- c("gene1", "gene2")
      genem <- genem |> dplyr::mutate(mix = grDevices::rgb(red = gene1, green = gene2, blue = 0, maxColorValue = max(genem[1:2])))
      plot <- merge(genem, obj@metadata, by = 0)
      legend <- expand.grid(red = seq(0, max(genem[1:2]), by = 0.02), green = seq(0, max(genem[1:2]), by = 0.02))
      legend <- within(legend, mix <- grDevices::rgb(green = green, red = red, blue = 0, maxColorValue = max(genem[1:2])))
      p1 <- ggplot2::ggplot(plot, ggplot2::aes(x = .data[[dim_x]], y = .data[[dim_y]], color = mix)) +
        ggplot2::geom_point(size = pointSize) +
        ggplot2::theme_classic() +
        ggplot2::guides(colour = guide_legend(override.aes = list(size=3))) +
        ggplot2::scale_color_identity() +
        ggplot2::labs(title = paste0("Methylation across ", genes[1], " and ", genes[2])) + noAxes()
      p2 <- ggplot2::ggplot(legend, ggplot2::aes(x=red, y=green)) +
        ggplot2::geom_tile(aes(fill=mix), color="white") +
        ggplot2::scale_fill_identity() +
        ggplot2::theme_classic() +
        ggplot2::labs(x = paste(genes[1]), y = paste(genes[2]), title = "Legend")
      p2 <- p2 + coord_fixed(ratio = 1/1)
      cowplot::plot_grid(p1, p2, rel_widths = c(3,1))
    }
    else if (length(genes) == 3) {
      names(genem) <- c("gene1", "gene2", "gene3")
      genem <- genem |> dplyr::mutate(mix = grDevices::rgb(red = gene1, green = gene2, blue = gene3, maxColorValue = max(genem[1:3])))
      plot <- merge(genem, obj@metadata, by = 0)
      p1 <- ggplot2::ggplot(plot, ggplot2::aes(x = .data[[dim_x]], y = .data[[dim_y]], color = mix)) +
        ggplot2::geom_point(size = pointSize) +
        ggplot2::theme_classic() +
        ggplot2::guides(colour = guide_legend(override.aes = list(size=3))) +
        ggplot2::scale_color_identity() +
        ggplot2::labs(title = paste0("Methylation across ", genes[1], " (red), ", genes[2], " (green), and ", genes[3], " (blue)")) + noAxes() +
        ggplot2::theme(legend.position = "none")
      legend <- expand.grid(red = seq(0, .1, by = 0.02), blue = seq(0, .1, by = 0.02), green = seq(0, .1, by = 0.02))
      legend <- within(legend, mix <- grDevices::rgb(green = green, red = red, blue = blue, maxColorValue = 0.1))

      p1

    } else {
      stop("Error: Blending only accommodates 2 or 3 genes")
    }
  }
}

############################################################################################################################
#' @title violinM
#' @description Generates a violin plot of percent methylation over a gene body. Best for non-CpG methylation.
#'
#' @param obj The amethyst object to plot
#' @param genes Genes to plot
#' @param type Type of methylation - e.g. "CH" or "CG" - to calculate if matrix is not provided
#' @param matrix Optional name of a pre-calculated matrix of values contained in the genomeMatrices slot to use
#' @param groupBy Parameter to group by on the x axis
#' @param colors Optional color palette
#' @param nrow Number of rows to plot if visualizing many genes
#' @return Returns ggplot2 geom_violin plots of methylation levels faceted by the requested genes
#' @export
#' @examples violinM(obj, genes = c("GRM8"), matrix = "gene_ch", groupBy = "cluster_id")
#' @importFrom dplyr inner_join
#' @importFrom ggplot2 ggplot aes geom_violin geom_jitter theme_classic scale_fill_manual scale_color_manual facet_wrap
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
violinM <-  function(
    obj,
    genes,
    matrix,
    groupBy,
    colors = NULL,
    pointSize = 0.2,
    nrow = round(sqrt(length(genes)))) {

  if (is.null(obj@genomeMatrices[[matrix]])) {
    stop("Please construct a matrix of methylation levels over the desired genes using makeWindows and store in the genomeMatrices slot.")
  }

  genem <- obj@genomeMatrices[[matrix]]
  genem <- genem[genes,] |> tibble::rownames_to_column(var = "gene")
  genem <- tidyr::pivot_longer(genem, cols = c(2:ncol(genem)), names_to = "cell_id", values_to = "pctm")

  # merge with metadata
  plot <- dplyr::inner_join(genem, obj@metadata |> tibble::rownames_to_column(var = "cell_id"), by = "cell_id")

  ggplot2::ggplot(plot, ggplot2::aes(x = .data[[groupBy]], y = pctm, fill = .data[[groupBy]], color = .data[[groupBy]])) +
    ggplot2::geom_violin(alpha = 0.3) +
    ggplot2::geom_jitter(ggplot2::aes(color = .data[[groupBy]]), size = pointSize) +
    ggplot2::theme_classic() +
    {if (!is.null(colors)) scale_fill_manual(values = colors) } +
    {if (!is.null(colors)) scale_color_manual(values = colors) } +
    ggplot2::facet_wrap(vars(gene), nrow = nrow)
}

############################################################################################################################
#' @title dotM
#' @description Generates a dot plot of percent methylation over a gene body
#'
#' @param obj The amethyst object to plot
#' @param genes  List of genes to plot methylation on the x axis
#' @param matrix Name of a pre-calculated methylation levels over genes of interest contained in the genomeMatrices slot
#' @param groupBy A categorical variable contained in the metadata to group cells by on the y axis
#' @param splitBy Optional additional facet level
#' @param nrow If splitBy is specified, this option controls the number of rows the facets are distributed across
#' @param colors Option to override default viridis scale
#'
#' @return Returns a ggplot object displaying average methylation values for each gene by the grouping variable
#' @export
#' @examples dotM(brain, genes = c("SATB2", "GAD1", "LINGO1"), matrix = "gene_ch", groupBy = "type")
#' @examples dotM(brain, genes = c("SATB2", "GAD1", "LINGO1"), matrix = "gene_ch", groupBy = "type", splitBy = "batch", colors = c("black", "red", "yellow"), nrow = 1)
#' @importFrom dplyr inner_join summarise group_by
#' @importFrom ggplot2 ggplot aes geom_point theme_classic scale_color_viridis_c theme facet_wrap
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
dotM <-  function(
    obj,
    genes,
    groupBy,
    splitBy = NULL,
    nrow = NULL,
    colors = NULL,
    matrix) {

  if (is.null(obj@genomeMatrices[[matrix]])) {
    stop("Please construct a matrix of methylation levels over the desired genes using makeWindows and store in the genomeMatrices slot.")
  }
  genem <- obj@genomeMatrices[[matrix]]
  genem <- genem[genes,] |> tibble::rownames_to_column(var = "gene")
  genem <- tidyr::pivot_longer(genem, cols = c(2:ncol(genem)), names_to = "cell_id", values_to = "pctm")
  genem <- dplyr::inner_join(genem, obj@metadata |> tibble::rownames_to_column(var = "cell_id"), by = "cell_id")

  if (is.null(splitBy)) {
    genem <- genem |> dplyr::group_by(.data[[groupBy]], gene) |> dplyr::summarise(pctm = mean(pctm, na.rm = TRUE))
  } else if (!is.null(splitBy)) {
    genem <- genem |> dplyr::group_by(.data[[groupBy]], .data[[splitBy]], gene) |> dplyr::summarise(pctm = mean(pctm, na.rm = TRUE))
  }
  genem$gene <- factor(genem$gene, levels = genes)

  ggplot2::ggplot(genem, ggplot2::aes(x = gene, y = .data[[groupBy]])) +
    ggplot2::geom_point(aes(size = pctm, color = pctm)) +
    ggplot2::theme_classic() +
    {if (is.null(colors)) ggplot2::scale_color_viridis_c()} +
    {if (!is.null(colors)) ggplot2::scale_color_gradientn(colors = {{ colors }})} +
    {if (!is.null(splitBy)) ggplot2::facet_wrap(vars(.data[[splitBy]]), nrow = nrow)} +
    ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

}

############################################################################################################################
#' @title histograM
#' @description Plot methylation levels over a gene body in histogram format
#'
#' @param obj Object containing the matrix to plot
#' @param genes Gene list to plot
#' @param matrix Matrix contained in the genomeMatrices slot; pct matrix output of calcSmoothedWindows
#' @param colors Optional list of colors to include
#' @param trackOverhang Number of base pairs to extend beyond the gene
#' @param legend Boolean indicating whether to plot the color legend
#' @param removeNA Boolean indicating whether to remove NA values as opposed to plotting them grey
#' @param ncol Number of columns to distribute results across. Defaults to number of input genes
#' @param arrowOverhang Number of base pairs the track arrow should extend beyond the gene
#' @param trackScale Enables adjustment of gene track height
#' @param arrowScale Enables adjustment of gene directional arrow size
#' @param colorMax Set upper bound on color scale. Everything exceeding this threshold will be plotted at the max value
#' @param order Arrange the order of groups
#' @param baseline 0 or "mean". Define whether the histogram baseline starts at 0 or the mean global % methylation of the group.
#' @return A ggplot geom_tile object with colors indicating % methylation over tiled windows and the gene of interest beneath
#' @export
#' @examples histograM(brain, genes = "SATB2", matrix = "cg_type_tracks")
#' @examples histograM(brain, genes = "SATB2", matrix = "ch_type_tracks", baseline = "mean", arrowScale = .03, trackScale = .5, colors = c("#ff4a74", "grey20", "#71e300"), colorMax = 12, trackOverhang = 50000, order = rev(c("Astro", "Oligo", "OPC", "Micro", "Inh_CGE", "Inh_MGE", "Exc_L1-3_CUX2", "Exc_L2-4_RORB", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4")))
#' @importFrom dplyr filter mutate rowwise
#' @importFrom ggplot2 ggplot geom_tile aes geom_rect geom_segment theme scale_fill_gradientn ylab arrow unit
#' @importFrom gridExtra grid.arrange
#' @importFrom scales squish
#' @importFrom tidyr pivot_longer
histograM <- function(obj,
                      genes = NULL,
                      matrix,
                      colors = NULL,
                      trackOverhang = 5000,
                      arrowOverhang = 3000,
                      ncol = length(genes),
                      legend = TRUE,
                      removeNA = TRUE,
                      trackScale = 1.5,
                      arrowScale = 0.025,
                      colorMax = 100,
                      order = NULL,
                      baseline = 0) {

  if (!is.null(colors)) {
    pal <- colors
  } else {
    pal <- c("#FF0082", "#dbdbdb", "#cccccc", "#999999")
  }

  if (is.null(obj@ref)) {
    stop("Please make sure a genome annotation file has been added to the obj@ref slot with makeRef().")
  }

  if (is.null(obj@genomeMatrices[[matrix]])) {
    stop("Specified matrix does not exist. Use names(obj@genomeMatrices) to check available matrices.")
  }

  p <- vector("list", length(genes)) # empty plot list
  for (i in 1:length(genes)) {

    ref <- obj@ref |> dplyr::filter(gene_name == genes[i])
    aggregated <- obj@genomeMatrices[[matrix]]

    toplot <- aggregated[c((aggregated$chr == ref$seqid[ref$type == "gene"] &
                              aggregated$start > (ref$start[ref$type == "gene"] - trackOverhang) &
                              aggregated$end < (ref$end[ref$type == "gene"] + trackOverhang))), ]
    ngroups <- ncol(toplot) - 3
    trackHeight <- ngroups * trackScale
    toplot <- tidyr::pivot_longer(toplot, cols = c(4:ncol(toplot)), names_to = "group", values_to = "pct_m") |> dplyr::rowwise() |> dplyr::mutate(middle = mean(c(start, end), na.rm = TRUE))

    if (removeNA) {
      toplot <- toplot |> dplyr::filter(!is.na(pct_m))
    }

    if (baseline == "mean") {
      glob_m <- data.frame(group = colnames(aggregated[, 4:ncol(aggregated)]),
                           glob_m = colMeans(aggregated[, 4:ncol(aggregated)], na.rm = T))
      toplot <- dplyr::left_join(toplot, glob_m, by = "group")
    }

    if (!is.null(order)) {
      toplot$group <- factor(toplot$group, levels = order)
    }

    p[[i]] <- ggplot2::ggplot() +

      # plotting gene tracks
      ggplot2::geom_rect(fill = "pink", data = ref |> dplyr::filter(type == "gene") |>
                           dplyr::mutate(promoter_start = ifelse(strand == "+", (start - 1500), (end + 1500)),
                                         promoter_end = ifelse(strand == "+", (promoter_start+3000), (promoter_start-3000))),
                         ggplot2::aes(xmin = promoter_start, xmax = promoter_end, ymin = -trackHeight*2, ymax = -trackHeight)) +
      ggplot2::geom_rect(fill = "black", data = ref |> dplyr::filter(type == "exon"),
                         ggplot2::aes(xmin = start, xmax = end, ymin = -trackHeight*2, ymax = -trackHeight)) +
      ggplot2::geom_segment(data = ref, aes(x = ifelse(strand == "+", (min(start) - arrowOverhang), (max(end)) + arrowOverhang),
                                            y = -(trackHeight*1.5),
                                            xend = ifelse(strand == "+", (max(end) + arrowOverhang), (min(start)) - arrowOverhang),
                                            yend = -(trackHeight*1.5)), arrow = ggplot2::arrow(length = unit(trackHeight*arrowScale, "cm"))) + xlab(genes[i]) +

      # plotting histogram
      {if (baseline == 0)ggplot2::geom_col(data = toplot, ggplot2::aes(x = middle, y = pct_m, fill = pct_m), width = mean(toplot$end - toplot$start))} +
      {if (baseline == "mean")ggplot2::geom_rect(data = toplot, ggplot2::aes(xmin = start, xmax = end, ymin = glob_m, ymax = pct_m, fill = pct_m))} + ggplot2::facet_grid(vars(group)) +
      {if (!legend)ggplot2::theme(legend.position = "none")} +
      ggplot2::scale_fill_gradientn(colors = pal, limits = c(0,colorMax), oob = scales::squish) + theme(axis.title.y = element_blank()) +
      ggplot2::theme(panel.background = element_blank(), axis.ticks = element_blank(), panel.grid.major.y = element_line(color = "#dbdbdb", linetype = "dashed"))
  }
  gridExtra::grid.arrange(grobs = p, ncol = ncol)
}

############################################################################################################################
#' @title heatMap
#' @description Plot methylation levels aggregated across smoothed genomic windows over a given gene body or genomic region
#'
#' @param obj Object containing the matrix to plot
#' @param genes Gene list to plot
#' @param matrix Matrix contained in the genomeMatrices slot. Note: unlike other functions, designed for makeSlidingWindows output
#' @param colors Optional list of colors to include
#' @param trackOverhang Number of base pairs to extend beyond the gene
#' @param regions Optional. Instead of genes, regions can be plotted, e.g. "chr1_1000_2000".
#' Any genes or exons that fall within the region will be plotted below.
#' @param nrow When multiple genes or regions are plotted, number of rows the output is divided across.
#' @param legend Boolean indicating whether to plot the color legend
#' @param removeNA Boolean indicating whether to remove NA values as opposed to plotting them grey
#' @param trackScale Enables adjustment of gene track height
#' @param arrowScale Enables adjustment of gene directional arrow size
#' @param colorMax Set upper bound on color scale. Everything exceeding this threshold will be plotted at the max value
#' @param order Arrange the order of groups
#' @param arrowOverhang Number of base pairs the track arrow should extend beyond the gene
#' @return A ggplot geom_tile object with colors indicating % methylation over tiled windows and the gene of interest beneath
#' @export
#' @examples heatMap(brain, genes = c("SLC17A7"), matrix = "cg_type_tracks")
#' @examples heatMap(brain, regions = c("chr2_199230000_199520000"), matrix = "ch_type_tracks", arrowScale = 0.3, trackScale = 0.07, colors = c("black", "red", "yellow"), colorMax = 20, order = c("Astro", "Oligo", "OPC", "Micro", "Inh_CGE", "Inh_MGE", "Exc_L1-3_CUX2", "Exc_L2-4_RORB", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4"))
#' @importFrom dplyr filter mutate rowwise cur_group_id
#' @importFrom ggplot2 ggplot geom_tile aes geom_rect geom_segment theme scale_fill_gradientn ylab arrow unit
#' @importFrom gridExtra grid.arrange
#' @importFrom tidyr pivot_longer
heatMap <- function(obj,
                    genes = NULL,
                    regions = NULL,
                    matrix,
                    colors = NULL,
                    trackOverhang = 5000,
                    arrowOverhang = 3000,
                    legend = TRUE,
                    removeNA = TRUE,
                    trackScale = .07,
                    arrowScale = 0.25,
                    colorMax = 100,
                    order = NULL,
                    nrow = max(length(genes), length(regions))) {

  if (!is.null(colors)) {
    pal <- colors
  } else {
    pal <- c("#FF0082", "#dbdbdb", "#cccccc", "#999999")
  }

  if (is.null(obj@ref)) {
    stop("Please make sure a genome annotation file has been added to the obj@ref slot with makeRef().")
  }

  if (!is.null(genes)) {
    # make empty plot list
    p <- vector("list", length(genes)) # empty plot list
    for (i in 1:length(genes)) {

      ref <- obj@ref |> dplyr::filter(gene_name == genes[i])
      aggregated <- obj@genomeMatrices[[matrix]]

      toplot <- aggregated[c((aggregated$chr == ref$seqid[ref$type == "gene"] & aggregated$start > (ref$start[ref$type == "gene"] - trackOverhang) & aggregated$end < (ref$end[ref$type == "gene"] + trackOverhang))), ]
      ngroups <- ncol(toplot) - 3
      trackHeight <- ngroups * trackScale
      toplot <- tidyr::pivot_longer(toplot, cols = c(4:ncol(toplot)), names_to = "group", values_to = "pct_m") |> dplyr::rowwise() |> dplyr::mutate(middle = mean(c(start, end), na.rm = TRUE))
      if (removeNA) {
        toplot <- toplot |> dplyr::filter(!is.na(pct_m))
      }
      if (!is.null(order)) {
        toplot$group <- factor(toplot$group, levels = order)
      }

      p[[i]] <- ggplot2::ggplot() +
        ggplot2::geom_tile(data = toplot, ggplot2::aes(x = middle, y = group, fill = pct_m), width = mean(toplot$end - toplot$start)) +
        ggplot2::geom_rect(fill = "pink", data = ref |> dplyr::filter(type == "gene") |> dplyr::mutate(promoter_start = ifelse(strand == "+", (start - 1500), (end + 1500)), promoter_end = ifelse(strand == "+", (promoter_start+3000), (promoter_start-3000))),
                           ggplot2::aes(xmin = promoter_start, xmax = promoter_end, ymin = -trackHeight, ymax = 0)) +
        ggplot2::geom_rect(fill = "black", data = ref |> dplyr::filter(type == "exon"), ggplot2::aes(xmin = start, xmax = end, ymin = -trackHeight, ymax = 0)) +
        ggplot2::geom_segment(data = ref, aes(x = ifelse(strand == "+", (min(start) - arrowOverhang), (max(end)) + arrowOverhang), y = -(trackHeight/2), xend = ifelse(strand == "+", (max(end) + arrowOverhang), (min(start)) - arrowOverhang), yend = - (trackHeight/2)), arrow = arrow(length = unit(trackHeight*arrowScale, "cm"))) + xlab(genes[i]) +
        ggplot2::geom_rect(fill = "black", data = ref |> dplyr::filter(type == "gene"), ggplot2::aes(xmin = start, xmax = end, ymin = (-(trackHeight/2) -(trackHeight/20)), ymax = (-(trackHeight/2) + (trackHeight/20)))) +
        {if (!legend)ggplot2::theme(legend.position = "none")} +
        ggplot2::scale_fill_gradientn(colors = pal, limits = c(0,colorMax), oob = scales::squish) +
        ggplot2::theme(panel.background = element_blank(), axis.ticks = element_blank(),  panel.grid.major.y = element_line(color = "#dbdbdb", linetype = "dashed")) + ggplot2::ylab("Group ID")
    }
  }

  if (!is.null(regions)) {
    p <- vector("list", length(regions))
    for (i in 1:length(regions)) {

      # get genome location
      split_parts <- data.table::tstrsplit(regions[[i]], "_", fixed = TRUE)
      chrom = split_parts[[1]]
      min = as.numeric(split_parts[[2]])
      max = as.numeric(split_parts[[3]])

      aggregated <- obj@genomeMatrices[[matrix]]
      toplot <- aggregated[(chr == chrom & start >= min & end <= max)]
      ngroups <- ncol(toplot) - 3

      ref <- obj@ref |> dplyr::filter(seqid == chrom & start >= min & end <= max & type %in% c("gene", "exon")) |>
        dplyr::distinct(gene_name, start, end, type, strand) |> dplyr::group_by(gene_name) |> dplyr::mutate(label = dplyr::cur_group_id()) |>
        dplyr::mutate(trackHeight = ngroups * label * trackScale,
                      ymin = -(trackHeight - (ngroups * 0.03)),
                      ymax = -(trackHeight + (ngroups * 0.03)))

      trackHeight <- mean(ref$trackHeight)

      toplot <- tidyr::pivot_longer(toplot, cols = c(4:ncol(toplot)), names_to = "group", values_to = "pct_m") |> dplyr::rowwise() |> dplyr::mutate(middle = mean(c(start, end), na.rm = TRUE))
      if (removeNA) {
        toplot <- toplot |> dplyr::filter(!is.na(pct_m))
      }
      if (!is.null(order)) {
        toplot$group <- factor(toplot$group, levels = order)
      }

      p[[i]] <- ggplot2::ggplot() +
        ggplot2::geom_tile(data = toplot, ggplot2::aes(x = middle, y = group, fill = pct_m), width = mean(toplot$end - toplot$start)) +
        {if (!legend)ggplot2::theme(legend.position = "none")} +
        ggplot2::scale_fill_gradientn(colors = pal,limits = c(0,colorMax), oob = scales::squish) +
        ggplot2::theme(panel.background = element_blank(), axis.ticks = element_blank(),  panel.grid.major.y = element_line(color = "#dbdbdb")) + ggplot2::ylab("Group ID") +
        ggplot2::xlab(chrom) +
        # plotting genes beneath
        ggplot2::geom_rect(fill = "pink", data = ref |> dplyr::filter(type == "gene") |> dplyr::mutate(promoter_start = ifelse(strand == "+", (start - 1500), (end + 1500)), promoter_end = ifelse(strand == "+", (promoter_start+3000), (promoter_start-3000))),
                           ggplot2::aes(xmin = promoter_start, xmax = promoter_end, ymin = ymin, ymax = ymax)) +
        ggplot2::geom_rect(fill = "black", data = ref |> dplyr::filter(type == "exon"),
                           ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax)) +
        ggplot2::geom_segment(data = ref |> dplyr::filter(type == "gene") |> dplyr::mutate(x = ifelse(strand == "+", start - arrowOverhang, end + arrowOverhang),
                                                                                           xend = ifelse(strand == "+", end + arrowOverhang, start - arrowOverhang)),
                              aes(x = x, xend = xend, y = -trackHeight, yend = -trackHeight),
                              arrow = ggplot2::arrow(length = unit(trackHeight*arrowScale, "cm"))) +
        geom_text(data = ref |> dplyr::group_by(gene_name) |> dplyr::arrange(desc(end)) |> dplyr::slice_head(n = 1),
                  aes(x = end + 100, y = -trackHeight, label = gene_name), hjust = 0)
    }
  }
  gridExtra::grid.arrange(grobs = p, nrow = nrow)
}

######################################################################
#' @title makePalette
#'
#' @description Generates a color palette of any length
#' @param option Which palette to choose out of options 1-30. Use testPalette to see the colors in each.
#' @param n How many colors are needed in the palette.
#' @return Returns a character vector of n hex codes.
#' @examples pal <- makePalette(option = 1, n = 20)
#' @export
#' @importFrom grDevices colorRampPalette
makePalette <- function(
    option,
    n) {
  if (option == 1) {pal <- c("#004A4A", "#419e9e", "#75C8D2", "#aae3e3", "#cfe3dc", "#d4c8b2", "mistyrose1", "#ffdea3", "#FFD554", "#fcba2b", "#FFA976",  "#FF8B73", "#F05252", "#bd083e")}
  if (option == 2) {pal <- c("#d84c4b", "#15485d", "#439f9e", "#f2e5dd", "#f08e29", "#df431a")}
  if (option == 3) {pal <- c("#150808","#ec5151","#ffdd6d","#a8dadc")}
  if (option == 4) {pal <- c("#293241","#3d5a80","#98c1d9","#e0fbfc","#ee6c4d")}
  if (option == 5) {pal <- c("#00050b", "#0b704e", "#6f9d80", "#e1e9d1", "#fb8b01", "#f34509")}
  if (option == 6) {pal <- c("#3D7577", "#4A4B1F", "#BDBEBF", "#F69E09", "#E32D12", "#540101")}
  if (option == 7) {pal <- c("#B5DCA5","#F9AB60","#E7576E", "#630661", "#220D50")}
  if (option == 8) {pal <- c("#0D353F","#72CDAE","#E6DAC6","#F5562A","#AB2E44")}
  if (option == 9) {pal <- c("#611c35","#a63446","#f44e3f","#ffa630","#f3d9dc","#d1c8e1","#2e5077","#373f51","#4da1a9", "#B4DDE1")}
  if (option == 10) {pal <- c("#fd5145","#ff7165","#ffbaa4","#87d0bf","#157d88", "#043E44")}
  if (option == 11) {pal <- c("#FBD8B0", "#DCF2C4", "#74DFD5", "#134077","#DF4275")}
  if (option == 12) {pal <- c("#c05761","#734f5a","#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51","#941c2f")}
  if (option == 13) {pal <- c("#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51")}
  if (option == 14) {pal <- c("#870022", "#db5375","#e86c5f","#f58549","#f2a65a","#eec170", "#ccd0b5","#bbd7d7", "#288989")}
  if (option == 15) {pal <- c("#f7b0be","#ed8e83","#ef3c23","#f15a42","#fac92c","#cfe5cc","#2178ae","#1b4793")}
  if (option == 16) {pal <- c("#50514f","#f25f5c","#ffe066","#247ba0","#70c1b3","#c0e8f9")}
  if (option == 17) {pal <- c("#C5C9BC", "#114F5A", "#D08222", "#402742", "#B69BB2")}
  if (option == 18) {pal <- c("#f7ede2","#f6bd60","#f5cac3","#84a59d","#f28482")}
  if (option == 19) {pal <- c("#218380","#73d2de","#ffbc42","#d81159","#8f2d56")}
  if (option == 20) {pal <- c("#383536", "#361F1C", "#A23226", "#D7412E", "#EDBBBC","#979296" )}
  if (option == 21) {pal <- c("#c46c9f", "#f46b72", "#f0898c", "#feb46b", "#fe9e6c", "#acd2c7", "#bed4bd", "#d3cd79")}
  if (option == 22) {pal <- c("#B7CFCF", "#D3C75C", "#8BA2CD", "#E2987B", "#37737D")}
  if (option == 23) {pal <- c("#50514f","#f25f5c","#ffe066","#247ba0","#70c1b3","#c0e8f9")}
  if (option == 24) {pal <- c("#7EAA9F", "#892B69", "#6E0B1B", "#C61B24", "#EE5C47", "#fcb9ac")}
  if (option == 25) {pal <- c("#7C2C47", "#B4141F", "#DF6C26", "#E9B60D", "#BABC36")}
  if (option == 26) {pal <- c("#1c5253", "#306b34","#c3eb78", "#f3ffc6","#b6174b")}
  if (option == 27) {pal <- c("#b6e2dd","#c8ddbb","#e9e5af","#fbdf9d","#fbc99d","#fbb39d","#fba09d")}
  if (option == 28) {pal <- c("#247ba0","#a15856","#f25f5c","#f9a061","#ffe066","#92ae83","#70c1b3","#4a9eaa","#50514f")}
  if (option == 29) {pal <- c("#579393", "#A95862", "#FB856C", "#FEB780", "#7D8EA8")}
  if (option == 30) {pal <- c("#C0CFE0", "#C5C463", "#F2AF3F","#FFA090", "#EF5356")}

  if (n < length(pal)) {
    colors <- sample(pal, size = n)
  } else {
    colors <- grDevices::colorRampPalette(pal)(n)
  }
  return(colors)
}

############################################################################################################################
#' @title testPalette
#' @description Visualize the colors in each palette option.
#'
#' @param output Either "swatch", which shows the colors in each palette in rows of tiles,
#' or "dimFeature", which shows your amethyst object colored by cluster_id in each palette option.
#' @param n If output = "swatch", how many colors should be shown with each palette option.
#' @param obj If output = "dimFeature", name of the amethyst object to test.
#' @return If output = "swatch", returns a ggplot object with each row as a palette.
#' If output = "dimFeature", returns the amethyst object in each palette option.
#' @export
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_tile theme_void scale_fill_identity facet_wrap ggtitle
#' @examples testPalette(output = "swatch", n = 30)
#' @examples testPalette(output = "dimFeature", obj = obj)
testPalette <- function(output,
                        n = NULL,
                        obj = NULL) {
  if (output == "swatch") {
    colors <- list()
    for (i in 1:30) {
      colors[[i]] <- makePalette(option = i, n = n)
    }
    colors <- as.data.frame(do.call(cbind, colors))
    colnames(colors) <- paste0(1:30)
    colors <- tidyr::gather(colors, key="key", value="color")

    p1 <- ggplot2::ggplot(colors, ggplot2::aes(x = color, y = key, fill = color)) +
      ggplot2::geom_tile() + ggplot2::theme_void() +
      ggplot2::scale_fill_identity() +
      ggplot2::facet_wrap(. ~ key, scales = "free", ncol = 1, strip.position = "left")
    return(p1)
  } else if (output == "dimFeature") {
    p <- vector("list", 30) # empty plot list
    for (i in 1:30) {
      colors <- sample(makePalette(option = i, n = length(unique(obj@metadata$cluster_id))))
      # get average methylation across a gene
      p[[i]] <- dimFeature(obj = obj, colorBy = cluster_id, pointSize = 0.1, colors = colors, reduction = "umap") +
        theme(legend.position = "none") + ggplot2::ggtitle(paste0("pal ", i))
    }
    gridExtra::grid.arrange(grobs = p, nrow = 5)
  }
}
