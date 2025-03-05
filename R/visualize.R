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
    {if (!is.null(splitBy)) ggplot2::facet_wrap(vars(.data[[splitBy]]), nrow = nrow, scales = "free_y")} +
    ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

}

############################################################################################################################
#' @title histograM
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
#' @param extraTracks Optional; include additional tracks (like regulatory elements) below gene tracks.
#' Input should be a named list of coordinate files containing chr, start, and end columns.
#' @param arrowOverhang Number of base pairs the track arrow should extend beyond the gene
#' @param remove Option to remove genes containing specific patterns from the plot, e.g. "^ENS|^LINC".
#' Helpful when plotting very large regions.
#' @param baseline Option for histogram to begin at 0 or at the "mean" methylation values for that group.
#' @param orientation Defines whether groups are plotted in "rows" or "cols".
#'
#' @return A ggplot object with color and bar height indicating % methylation over windows.
#' Genes and any imported tracks (optional) are plotted beneath.
#' @export
#' @examples histograM(brain, matrix = "ch_type_tracks", genes = "SATB2")
#' @examples histograM(brain, matrix = "ch_type_tracks", regions = c("chr2_170813213_170861151", "chr2_199269505_199471266"), orientation = "cols", trackOverhang = 1000000, remove = "ENS|LINC")
#' @importFrom dplyr filter mutate rowwise cur_group_id
#' @importFrom ggplot2 ggplot geom_tile aes geom_rect geom_segment theme scale_fill_gradientn ylab arrow unit annotate
#' @importFrom gridExtra grid.arrange
#' @importFrom tidyr pivot_longer
histograM <- function(obj,
                      matrix,
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
                      extraTracks = NULL,
                      remove = NULL,
                      baseline = "mean",
                      nrow = max(length(genes), length(regions)),
                      orientation = "rows") {

  # check correct matrix exists
  if (is.null(obj@genomeMatrices[[matrix]])) {
    stop("Specified matrix does not exist. Use names(obj@genomeMatrices) to check available matrices.")
  }

  # check to make sure reference is added
  if (is.null(obj@ref)) {
    stop("Please make sure a genome annotation file has been added to the obj@ref slot with makeRef().")
  }

  if (!"data.table" %in% class(obj@ref)) {
    obj@ref <- data.table(obj@ref)
  }

  if (!is.null(remove)) {
    obj@ref <- obj@ref[!grepl(x = obj@ref$gene_name, pattern = remove), ]
  }

  # check tracks are in proper format
  if (!is.null(extraTracks)) {
    for (i in 1:length(extraTracks)) {
      if (!(all(c("chr", "start", "end") %in% names(extraTracks[[i]])))) {
        stop("Please check each additional track contains chr, start, and end columns.")
      }
      if (!(is.numeric(extraTracks[[i]]$start) && is.numeric(extraTracks[[i]]$end))) {
        stop("Please check all start and end columns contain numeric values.")
      }
    }
  }

  # convert genes to genome ranges
  if (!is.null(regions)) {
    ranges <- regions
  }

  if (!is.null(genes)) {
    ref <- obj@ref[gene_name %in% genes & type %in% c("gene", "exon")]
    ranges <- ref[type == "gene", .SD[which.max(end - start)], by = gene_name][order(match(gene_name, genes))]$location
  }

  # for each genome range...
  p <- vector("list", length(ranges))
  for (i in 1:length(ranges)) {

    # get genome coordinates for this iteration
    split_parts <- data.table::tstrsplit(ranges[[i]], "_", fixed = TRUE)
    chrom = split_parts[[1]]
    min = as.numeric(split_parts[[2]]) - trackOverhang
    max = as.numeric(split_parts[[3]]) + trackOverhang

    # get methylation values from genome matrix
    values <- obj@genomeMatrices[[matrix]][(chr == chrom & start >= min & end <= max)]
    ngroups <- ncol(values) - 3

    # calculate track heights
    if (!is.null(regions)) {
      ref <- obj@ref |> dplyr::filter(seqid == chrom & start >= min & end <= max & type %in% c("gene", "exon")) |>
        dplyr::distinct(gene_name, start, end, type, strand) |> dplyr::group_by(gene_name) |> dplyr::mutate(label = dplyr::cur_group_id()) |>
        dplyr::mutate(trackHeight = ngroups * label * trackScale,
                      ymin = -(trackHeight - (ngroups * 0.03)),
                      ymax = -(trackHeight + (ngroups * 0.03)))
    }

    if (!is.null(genes)) {
      ref <- obj@ref |> dplyr::filter(gene_name == genes[i] & type %in% c("gene", "exon")) |>
        dplyr::distinct(gene_name, start, end, type, strand) |>
        dplyr::mutate(trackHeight = ngroups * trackScale,
                      ymin = -(trackHeight - (ngroups * 0.03)),
                      ymax = -(trackHeight + (ngroups * 0.03)))
    }

    # isolate gene coordinates to make plotting more legible
    genenames <- ref |> dplyr::group_by(gene_name) |> dplyr::arrange(desc(end)) |> dplyr::slice_head(n = 1) |> dplyr::mutate(text_start = ifelse(strand == "+", (end + 2500), (end + 1000)))
    genebody <- ref |> dplyr::filter(type == "gene") |> dplyr::mutate(x = ifelse(strand == "+", start, end),
                                                                      xend = ifelse(strand == "+", end + arrowOverhang, start - arrowOverhang))
    promoters <- ref |> dplyr::filter(type == "gene") |> dplyr::mutate(promoter_start = ifelse(strand == "+", (start - 1500), (end - 1500)),
                                                                       promoter_end = ifelse(strand == "+", (start + 1500), (end + 1500)))
    exons <- ref |> dplyr::filter(type == "exon")

    if (!is.null(arrowScale)) {
      arrowHeight <- arrowScale
    } else {
      arrowHeight <- ( 5 / (ngroups*nrow(genebody)))
    }

    values <- tidyr::pivot_longer(values, cols = c(4:ncol(values)), names_to = "group", values_to = "pct_m") |>
      dplyr::rowwise() |> dplyr::mutate(middle = mean(c(start, end), na.rm = TRUE))

    if (baseline == "mean") {
      glob_m <- data.frame(group = colnames(obj@genomeMatrices[[matrix]][, 4:ncol(obj@genomeMatrices[[matrix]])]),
                           glob_m = colMeans(obj@genomeMatrices[[matrix]][, 4:ncol(obj@genomeMatrices[[matrix]])], na.rm = T))
      values <- dplyr::left_join(values, glob_m, by = "group")
    }

    if (removeNA) {
      values <- values |> dplyr::filter(!is.na(pct_m))
    }
    if (!is.null(order)) {
      values$group <- factor(values$group, levels = order)
    }

    # define other aesthetics that change based on data
    trackHeight <- mean(ref$trackHeight, na.rm = T)
    width = mean(values$end - values$start, na.rm = T)

    if (!is.null(colors)) {
      pal <- colors
    } else if (is.null(colors)) {
      pal <- c("#0073bf", "grey60", "#ffa600")
    }

    # define scale
    if (!is.null(colorMax)) {
      colorMax <- colorMax
    } else if (is.null(colorMax)) {
      if (mean(values$pct_m, na.rm = T) > 10) { # assume it's CG methylation; always use 100
        colorMax <- 100
      } else { # assume it's CH methylation
        colorMax <- as.numeric(quantile(values$pct_m,probs=c(.999), na.rm = T))
      }
    }

    p[[i]] <- ggplot2::ggplot() +

      # plotting histogram
      {if (baseline == 0) ggplot2::geom_col(data = values, ggplot2::aes(x = middle, y = pct_m, fill = pct_m), width = mean(toplot$end - toplot$start))} +
      {if (baseline == "mean") ggplot2::geom_rect(data = values, ggplot2::aes(xmin = start, xmax = end, ymin = glob_m, ymax = pct_m, fill = pct_m))} +
      {if (orientation == "rows") ggplot2::facet_grid(rows = vars(group))} +
      {if (orientation == "cols") ggplot2::facet_grid(cols = vars(group))} +
      {if (!legend) ggplot2::theme(legend.position = "none")} +
      ggplot2::scale_fill_gradientn(colors = pal, limits = c(0,colorMax), oob = scales::squish) + theme(axis.title.y = element_blank()) +
      ggplot2::theme(panel.background = element_blank(), axis.ticks = element_blank(), panel.grid.major.y = element_line(color = "#dbdbdb", linetype = "dashed")) +

      # plotting genes beneath
      ggplot2::geom_rect(data = promoters, fill = "pink", ggplot2::aes(xmin = promoter_start, xmax = promoter_end, ymin = ymin, ymax = ymax)) +
      ggplot2::geom_rect(data = exons, fill = "black", ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax)) +
      ggplot2::geom_segment(data = genebody, aes(x = x, xend = xend, y = - trackHeight, yend = - trackHeight),
                            arrow = ggplot2::arrow(length = unit(trackHeight * arrowHeight , "cm"))) +
      {if (!is.null(regions)) ggplot2::geom_text(data = genenames, aes(x = text_start, y = -trackHeight, label = gene_name), hjust = 0) } +

      # dynamically add extra tracks
      {
        if (!is.null(extraTracks)) {
          lapply(seq_along(extraTracks), function(j) {
            track <- extraTracks[[j]] |>
              dplyr::filter(chr == chrom & start >= min & end <= max) |>
              dplyr::mutate(
                trackHeight = ngroups * (n_groups(ref) + j) * trackScale,
                ymin = -(trackHeight - (ngroups * 0.03)),
                ymax = -(trackHeight + (ngroups * 0.03)),
                track_name = names(extraTracks[j])
              )
            set.seed(222)
            ggplot2::geom_rect(data = track, fill = makePalette(option = 13, n = length(extraTracks))[j],
                               ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax))
          })
        }
      } +

      # dynamically add names (for some reason couldn't combine)
      {
        if (!is.null(extraTracks) && !is.null(names(extraTracks))) {
          trackNames <- data.frame(
            min = min,
            trackName = names(extraTracks),
            trackHeight = seq_along(extraTracks) %>%
              sapply(function(j) (ngroups * (n_groups(ref) + j) * trackScale) * -1)
          )
          ggplot2::geom_text(data = trackNames, aes(x = min, y = trackHeight, label = trackName))
        }
      } +

      # define additional aesthetics
      theme(legend.justification = "top") +
      {if (!legend)ggplot2::theme(legend.position = "none")} +
      ggplot2::theme(panel.background = element_blank(),
                     axis.ticks = element_blank(),
                     panel.grid.major.y = element_line(color = "#f0f0f0", linetype = "dashed")) +
      ggplot2::ylab("Group ID") +
      ggplot2::scale_fill_gradientn(colors = pal,limits = c(0,colorMax), oob = scales::squish) +
      {if (!is.null(regions)) ggplot2::xlab(paste0(chrom, ": ", min, " - ", max)) } +
      {if (!is.null(genes)) ggplot2::xlab(paste0(genes[i], " (", chrom, ")")) }

  }
  gridExtra::grid.arrange(grobs = p, nrow = nrow)
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
#' Also, can specify "genome". This option is meant to be run with large genomic windows.
#' @param nrow When multiple genes or regions are plotted, number of rows the output is divided across.
#' @param legend Boolean indicating whether to plot the color legend
#' @param removeNA Boolean indicating whether to remove NA values as opposed to plotting them grey
#' @param trackScale Enables adjustment of gene track height
#' @param arrowScale Enables adjustment of gene directional arrow size
#' @param colorMax Set upper bound on color scale. Everything exceeding this threshold will be plotted at the max value
#' @param order Arrange the order of groups
#' @param extraTracks Optional; include additional tracks (like regulatory elements) below gene tracks.
#' Input should be a named list of coordinate files containing chr, start, and end columns.
#' @param arrowOverhang Number of base pairs the track arrow should extend beyond the gene
#' @param remove Option to remove genes containing specific patterns from the plot, e.g. "^ENS|^LINC".
#' Helpful when plotting very large regions.
#'
#' @return A ggplot geom_tile object with colors indicating % methylation over tiled windows and the gene of interest beneath
#' @export
#' @examples heatMap(brain, genes = c("SLC17A7"), matrix = "cg_type_tracks")
#' @examples heatMap(brain, regions = c("chr2_199230000_199520000"), matrix = "ch_type_tracks", arrowScale = 0.3, trackScale = 0.07, colors = c("black", "red", "yellow"), colorMax = 20, order = c("Astro", "Oligo", "OPC", "Micro", "Inh_CGE", "Inh_MGE", "Exc_L1-3_CUX2", "Exc_L2-4_RORB", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4"))
#' @importFrom dplyr filter mutate rowwise cur_group_id
#' @importFrom ggplot2 ggplot geom_tile aes geom_rect geom_segment theme scale_fill_gradientn ylab arrow unit annotate
#' @importFrom gridExtra grid.arrange
#' @importFrom tidyr pivot_longer
heatMap <- function(obj,
                    matrix,
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
                    extraTracks = NULL,
                    remove = NULL,
                    nrow = max(length(genes), length(regions))) {

  # check correct matrix exists
  if (is.null(obj@genomeMatrices[[matrix]])) {
    stop("Specified matrix does not exist. Use names(obj@genomeMatrices) to check available matrices.")
  }

  if (!is.null(regions) && "genome" %in% regions) {
    if (ncol(obj@genomeMatrices[[matrix]]) > 1000) {
      stop ("Please consider using an aggregated matrix (see aggregateMatrix) for genome-wide observation.")
    }
    if (nrow(obj@genomeMatrices[[matrix]]) > 100000) {
      cat("Warning: please consider using windows with a larger resolution.")
    }

    # get methylation values from genome matrix
    values <- obj@genomeMatrices[[matrix]]
    values$location <- rownames(values)
    values <- as.data.table(values)
    values <- values[, c("chr", "start", "end") := data.table::tstrsplit(location, "_", fixed = TRUE)][, location := NULL]
    values <- values[, `:=`(start = as.numeric(start), end = as.numeric(end))][, `:=` (middle = (start + end) / 2)]
    values <- values |>
      tidyr::pivot_longer(cols = -c("chr", "start", "end", "middle"),  names_to = "group", values_to = "pct_m") |>
      dplyr::filter(!is.na(pct_m))

    # sort chr levels
    chr_levels <- unique(gsub("chr", "", values$chr))
    values$chr <- factor(values$chr,
                         levels = paste0("chr", c(sort(as.numeric(chr_levels[grepl("^[0-9]", chr_levels)])), chr_levels[!grepl("^[0-9]", chr_levels)])))

    # sort groups
    if (!is.null(order)) {
      values$group <- factor(values$group, levels = order)
    }

    width = mean(values$end - values$start, na.rm = T)

    # determine appropriate scale range
    if (!is.null(colorMax)) {
      range <- c(0, colorMax)
    } else if (is.null(colorMax)) {
      if (mean(values$pct_m, na.rm = T) > 10) { # assume it's CG methylation; always use 100
        range <- c(0, 100)
      } else if (mean(values$pct_m, na.rm = T ) < 10 && all(values$pct_m >= 0, na.rm = TRUE)) { # assume it's CH methylation
        colorMax <- as.numeric(quantile(values$pct_m,probs=c(.999), na.rm = T))
        range <- c(0, colorMax)
      } else { # assume matrix is methylation score
        range <- c(-1, 1)
      }
    }

    # determine appropriate color palette
    if (!is.null(colors)) {
      pal <- colors
    } else if (is.null(colors)) {
      if (mean(values$pct_m, na.rm = T) > 10) { # assume it's CG methylation; always use 100
        pal <- c("#FF0082", "#dbdbdb", "#cccccc", "#999999")
      } else if (mean(values$pct_m, na.rm = T ) < 10 && all(values$pct_m >= 0, na.rm = TRUE)) { # assume it's CH methylation
        pal <- c("black", "red", "yellow")
      } else { # assume matrix is methylation score
        pal <- c("#FF0082", "#dbdbdb", "#787878")
      }
    }

    # plot
    ggplot(values) +
      ggplot2::geom_vline(xintercept = seq(0, max(values$middle), by = 10e6), color = "grey", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
      ggplot2::geom_tile(ggplot2::aes(x = middle, y = group, fill = pct_m), width = 100000) +
      ggplot2::theme_classic() +
      ggplot2::scale_fill_gradientn(colors = pal,limits = range, oob = scales::squish) +
      facet_grid(rows = vars(chr), space = "free")

  } else { # if regions != genome
    # check to make sure reference is added
    if (is.null(obj@ref)) {
      stop("Please make sure a genome annotation file has been added to the obj@ref slot with makeRef().")
    }

    if (!"data.table" %in% class(obj@ref)) {
      obj@ref <- data.table(obj@ref)
    }

    if (!is.null(remove)) {
      obj@ref <- obj@ref[!grepl(x = obj@ref$gene_name, pattern = remove), ]
    }

    # check tracks are in proper format
    if (!is.null(extraTracks)) {
      for (i in 1:length(extraTracks)) {
        if (!(all(c("chr", "start", "end") %in% names(extraTracks[[i]])))) {
          stop("Please check each additional track contains chr, start, and end columns.")
        }
        if (!(is.numeric(extraTracks[[i]]$start) && is.numeric(extraTracks[[i]]$end))) {
          stop("Please check all start and end columns contain numeric values.")
        }
      }
    }

    # convert genes to genome ranges
    if (!is.null(regions)) {
      ranges <- regions
    }

    if (!is.null(genes)) {
      ref <- obj@ref[gene_name %in% genes & type %in% c("gene", "exon")]
      ranges <- ref[type == "gene", .SD[which.max(end - start)], by = gene_name][order(match(gene_name, genes))]$location
    }

    # for each genome range...
    p <- vector("list", length(ranges))
    for (i in 1:length(ranges)) {

      # get genome coordinates for this iteration
      split_parts <- data.table::tstrsplit(ranges[[i]], "_", fixed = TRUE)
      chrom = split_parts[[1]]
      min = as.numeric(split_parts[[2]]) - trackOverhang
      max = as.numeric(split_parts[[3]]) + trackOverhang

      # get methylation values from genome matrix
      values <- obj@genomeMatrices[[matrix]][(chr == chrom & start >= min & end <= max)]
      ngroups <- ncol(values) - 3

      # calculate track heights
      if (!is.null(regions)) {
        ref <- obj@ref |> dplyr::filter(seqid == chrom & start >= min & end <= max & type %in% c("gene", "exon")) |>
          dplyr::distinct(gene_name, start, end, type, strand) |> dplyr::group_by(gene_name) |> dplyr::mutate(label = dplyr::cur_group_id()) |>
          dplyr::mutate(trackHeight = ngroups * label * trackScale,
                        ymin = -(trackHeight - (ngroups * 0.03)),
                        ymax = -(trackHeight + (ngroups * 0.03)))
      }

      if (!is.null(genes)) {
        ref <- obj@ref |> dplyr::filter(gene_name == genes[i] & type %in% c("gene", "exon")) |>
          dplyr::distinct(gene_name, start, end, type, strand) |>
          dplyr::mutate(trackHeight = ngroups * trackScale,
                        ymin = -(trackHeight - (ngroups * 0.03)),
                        ymax = -(trackHeight + (ngroups * 0.03)))
      }

      # isolate gene coordinates to make plotting more legible
      genenames <- ref |> dplyr::group_by(gene_name) |> dplyr::arrange(desc(end)) |> dplyr::slice_head(n = 1) |> dplyr::mutate(text_start = ifelse(strand == "+", (end + 2500), (end + 1000)))
      genebody <- ref |> dplyr::filter(type == "gene") |> dplyr::mutate(x = ifelse(strand == "+", start, end),
                                                                        xend = ifelse(strand == "+", end + arrowOverhang, start - arrowOverhang))
      promoters <- ref |> dplyr::filter(type == "gene") |> dplyr::mutate(promoter_start = ifelse(strand == "+", (start - 1500), (end - 1500)),
                                                                         promoter_end = ifelse(strand == "+", (start + 1500), (end + 1500)))
      exons <- ref |> dplyr::filter(type == "exon")

      if (!is.null(arrowScale)) {
        arrowHeight <- arrowScale
      } else {
        arrowHeight <- ( 5 / (ngroups*nrow(genebody)))
      }

      values <- tidyr::pivot_longer(values, cols = c(4:ncol(values)), names_to = "group", values_to = "pct_m") |>
        dplyr::rowwise() |> dplyr::mutate(middle = mean(c(start, end), na.rm = TRUE))

      if (removeNA) {
        values <- values |> dplyr::filter(!is.na(pct_m))
      }
      if (!is.null(order)) {
        values$group <- factor(values$group, levels = order)
      }

      # define other aesthetics that change based on data
      trackHeight <- mean(ref$trackHeight, na.rm = T)
      width = mean(values$end - values$start, na.rm = T)

      if (!is.null(colors)) {
        pal <- colors
      } else if (is.null(colors)) {
        if (mean(values$pct_m, na.rm = T) > 10) { # Assume it's CG methylation
          pal <- c("#FF0082", "#dbdbdb", "#cccccc", "#999999")
        } else { # Assume it's CH methylation
          pal <- c("black", "red", "yellow")
        }
      }

      # define scale
      if (!is.null(colorMax)) {
        colorMax <- colorMax
      } else if (is.null(colorMax)) {
        if (mean(values$pct_m, na.rm = T) > 10) { # assume it's CG methylation; always use 100
          colorMax <- 100
        } else { # assume it's CH methylation
          colorMax <- as.numeric(quantile(values$pct_m,probs=c(.999), na.rm = T))
        }
      }

      p[[i]] <- ggplot2::ggplot() +

        # plot methylation values per group
        ggplot2::geom_tile(data = values, ggplot2::aes(x = middle, y = group, fill = pct_m), width = width) +

        # plotting genes beneath
        ggplot2::geom_rect(data = promoters, fill = "pink", ggplot2::aes(xmin = promoter_start, xmax = promoter_end, ymin = ymin, ymax = ymax)) +
        ggplot2::geom_rect(data = exons, fill = "black", ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax)) +
        ggplot2::geom_segment(data = genebody, aes(x = x, xend = xend, y = - trackHeight, yend = - trackHeight),
                              arrow = ggplot2::arrow(length = unit(trackHeight * arrowHeight , "cm"))) +
        {if (!is.null(regions)) ggplot2::geom_text(data = genenames, aes(x = text_start, y = -trackHeight, label = gene_name), hjust = 0) } +

        # dynamically add extra tracks
        {
          if (!is.null(extraTracks)) {
            lapply(seq_along(extraTracks), function(j) {
              track <- extraTracks[[j]] |>
                dplyr::filter(chr == chrom & start >= min & end <= max) |>
                dplyr::mutate(
                  trackHeight = ngroups * (n_groups(ref) + j) * trackScale,
                  ymin = -(trackHeight - (ngroups * 0.03)),
                  ymax = -(trackHeight + (ngroups * 0.03)),
                  track_name = names(extraTracks[j])
                )
              set.seed(222)
              ggplot2::geom_rect(data = track, fill = makePalette(option = 13, n = length(extraTracks))[j],
                                 ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax))
            })
          }
        } +

        # dynamically add names (for some reason couldn't combine)
        {
          if (!is.null(extraTracks) && !is.null(names(extraTracks))) {
            trackNames <- data.frame(
              min = min,
              trackName = names(extraTracks),
              trackHeight = seq_along(extraTracks) %>%
                sapply(function(j) (ngroups * (n_groups(ref) + j) * trackScale) * -1)
            )
            ggplot2::geom_text(data = trackNames, aes(x = min, y = trackHeight, label = trackName))
          }
        } +

        # define additional aesthetics
        theme(legend.justification = "top") +
        {if (!legend)ggplot2::theme(legend.position = "none")} +
        ggplot2::theme(panel.background = element_blank(),
                       axis.ticks = element_blank(),
                       panel.grid.major.y = element_line(color = "#f0f0f0", linetype = "dashed")) +
        ggplot2::ylab("Group ID") +
        ggplot2::scale_fill_gradientn(colors = pal,limits = c(0,colorMax), oob = scales::squish) +
        {if (!is.null(regions)) ggplot2::xlab(paste0(chrom, ": ", min, " - ", max)) } +
        {if (!is.null(genes)) ggplot2::xlab(paste0(genes[i], " (", chrom, ")")) }

    }
    gridExtra::grid.arrange(grobs = p, nrow = nrow)
  }
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
    n,
    sample = T) {
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

  if (n < length(pal) && sample) {
    colors <- sample(pal, size = n)
  } else if (n < length(pal) && !sample) {
    colors <- pal[1:n]
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
