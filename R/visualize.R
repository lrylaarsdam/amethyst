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
#' @title markerHeatmap
#' NEEDS UPDATING
#' @description Plot the methylation levels of sample groups across many genes at once
#' @param obj The amethyst object to plot
#' @param level Hierarchy of reference to plot. i.e., "cell_class" plots "Exc"  "Inh"  "NonN"
#' @param impute Logical indicating whether to impute values using Rmagic
#' @param rownames Logical indicating whether to show row names on the heatmap
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
#' @return Returns a heatmap of average methylation levels of key marker genes
#' @export
#' @examples markerHeatmap(obj = obj, level = "CellClass", impute = FALSE, rownames = FALSE, type = "CH")
#' @importFrom dplyr filter group_by arrange rename mutate select
#' @importFrom janitor make_clean_names
#' @importFrom Rmagic magic
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr pivot_wider
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
    pull(marker)
  markers <- makeWindows(obj, genes = markers, type = {{ type }})
  markers[is.na(markers)] <- 0
  if (impute) {
    markers <- Rmagic::magic(markers , npca = 10)
    markers <- markers[["result"]]
  }
  markers <- pivot_longer(markers %>% rownames_to_column(var = "cell_id"), cols = c(2:(length(colnames(markers)) + 1)), names_to = "gene", values_to = "pct_m")
  markers <- inner_join(markers, keygenestable %>% dplyr::filter(cluster_level == level) %>% dplyr::rename("gene" = "marker"), by = "gene")
  markers <- inner_join(markers, obj@metadata %>% rownames_to_column(var = "cell_id"), by = "cell_id")
  markers <- markers %>% rowwise() %>% dplyr::mutate(group = paste0(janitor::make_clean_names(description), "_", gene)) %>% dplyr::arrange(cluster_id, description)
  toplot <- tidyr::pivot_wider(markers %>% dplyr::select(c(cell_id, pct_m, group)), names_from = "group", values_from = "pct_m") %>% tibble::column_to_rownames(var = "cell_id")

  # plot pheatmap
  annot_col <- markers %>% dplyr::select(group, name) %>% distinct(group, name) %>% column_to_rownames("group")
  annot_row <- markers %>% dplyr::select(cell_id, cluster_id) %>% distinct(cell_id, cluster_id) %>% column_to_rownames("cell_id")

  pheatmap(toplot, annotation_row = annot_row, annotation_col = annot_col, show_rownames = rownames,
           labels_col = (markers %>% dplyr::select(group, name, gene) %>% distinct(group, name, .keep_all = TRUE) %>% pull("gene")),
           cluster_rows = F, cluster_cols = F)

}

############################################################################################################################
#' @title clusterCompare
#' NEEDS UPDATING
#' @description Correlates average percent methylation over 100kb windows aggregated by metadata to a reference matrix to aid in cell type identification.
#' @param matrix Average percent methylation aggregated over 100kb windows by a grouping factor (see makeWindows). Column names should be genomic positions.
#' @param ref The reference structure to correlate to. The default for brain used here is from PMC8494641.
#' @param level The grouping level in the reference to compare to
#' @param type The type of methylation in the reference to compare to (cg or ch). Make sure the input matrix reflects the same type.
#' @param method Correlation coefficient to be computed ("pearson", "kendall", or "spearman")
#' @param n Number of top correlations to plot.
#' @return A ggplot dot plot showing the top 5 correlations of each group to the reference
#' @export
#' @examples clusterCompare(cluster_100kb_cg, level = "major_type", type = "cg")
#' @importFrom dplyr bind_rows group_by arrange filter
#' @importFrom ggplot2 ggplot
#' @importFrom stats cor
#' @importFrom tidyr pivot_longer
clusterCompare <- function(
    matrix,
    ref = celltyperefs,
    level,
    type,
    method = "pearson",
    n = 5,
    aggregate = NULL) {

  matrix <- obj@genomeMatrices[[matrix]]
  withref <- dplyr::bind_rows(matrix, ref[[level]][[type]]) # append the reference matrix to the experimental
  cor <- stats::cor(t(withref), method = {{ method }}, use = "pairwise.complete.obs") # generate correlation matrix
  # determine top correlations
  annotation <- tidyr::pivot_longer(as.data.frame(cor) %>%
                                      rownames_to_column("group1"), cols = c(2:(length(colnames(cor))+1)), names_to = "group2", values_to = "cor") %>%
    dplyr::group_by(group1) %>%
    dplyr::arrange(desc(cor), .by_group = TRUE) %>%
    dplyr::filter((group1 %in% rownames(matrix)) & !(group2 %in% rownames(matrix))) %>%
    top_n(n = {{ n }})

  p <- ggplot2::ggplot(annotation, ggplot2:aes(x = group1, y = group2, size = cor, color = cor)) +
    ggplot2:geom_point() +
    ggplot2:theme_classic() +
    ggplot2:scale_color_viridis_c() +
    ggplot2:labs(x = "cluster", y = paste0("Ref ", type, " methylation by ", level)) +
    ggplot2:theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p
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
#' @title umapFeature
#'
#' @description Plot any feature in the metadata in UMAP space
#' @param obj The amethyst object to plot
#' @param colorBy The metadata feature to plot in UMAP space
#' @param colors Optional color scale
#' @param pointSize Optional adjustment of point size
#' @return Returns a ggplot graph plotting xy coordinates of cells colored according to the specified feature
#' @export
#' @examples umapFeature(obj = obj, colorBy = "cluster_id")
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_classic guides guide_legend element_blank
umapFeature <- function(
    obj,
    colorBy,
    colors = NULL,
    pointSize = 0.5) {
  p <- ggplot2::ggplot(obj@metadata, ggplot2::aes(x = umap_x, y = umap_y, color = {{colorBy}})) +
    ggplot2::geom_point(size = pointSize) +
    {if (!is.null(colors)) ggplot2::scale_color_manual(values = {{colors}})} +
    ggplot2::theme_classic() +
    noAxes() +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)))
  p
}

############################################################################################################################
#' @title umapGeneM
#' @description Plot cumulative % methylation over a gene body on umap
#'
#' @param obj The amethyst object to plot
#' @param genes List of genes to plot % methylation in UMAP space
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
#' @param blend Logical indicating whether to blend methylation levels of two or three genes
#' @param nrow Number of rows to distribute graphs over
#' @param metric If matrix is not specified, what methylation parameter to calculate, i.e. "percent", "score", or "ratio"
#' @param colors Optional color gradient to include
#' @param matrix Pre-computed matrix to use in the genomeMatrices slot, if available. The type parameter is not needed if matrix is definied.
#' @return Returns a plot of methylation over given genomic region in relation to UMAP
#' @export
#' @examples umapGeneM(obj, genes = c("GAD1", "SATB2"), type = "CH", metric = "percent", colors = c("black", "turquoise", "gold", "red"), blend = F)
#' @examples umapGeneM(both, genes = c("GAD1", "SATB2"), matrix = "gene_ch", blend = T, nrow = 1)
#' @importFrom dplyr inner_join mutate
#' @importFrom ggplot2 ggplot aes geom_point theme_classic guides scale_color_viridis_c scale_color_gradientn labs scale_color_identity geom_tile scale_fill_identity theme
#' @importFrom grDevices rgb
#' @importFrom gridExtra grid.arrange
#' @importFrom plotly plot_ly subplot
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom scales squish
umapGeneM <-  function(
    obj,
    genes,
    type = NULL,
    metric = "percent",
    matrix = NULL,
    blend = FALSE,
    colors = NULL,
    squish = NULL,
    nrow = round(sqrt(length(genes)))) {

  if (is.null(matrix) && is.null(type)) {
    stop("If matrix is not provided, type and metric must both be specified. Default metric is % methylation.")
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
      if (!is.null(matrix)) {
        genem <- obj@genomeMatrices[[matrix]]
        genem <- genem[genes[i],] |> tibble::rownames_to_column(var = "gene")
        genem <- tidyr::pivot_longer(genem, cols = c(2:ncol(genem)), names_to = "cell_id", values_to = "pctm")
        genem <- genem |> dplyr::filter(!is.na(pctm))
      }
      if (is.null(matrix)) {
        genem <- makeWindows({{obj}}, genes = genes[i], type = {{type}}, metric = {{metric}}) |> tibble::rownames_to_column(var = "gene")
        genem <- tidyr::pivot_longer(genem, cols = c(2:ncol(genem)), names_to = "cell_id", values_to = "pctm")
        genem <- genem |> dplyr::filter(!is.na(pctm))
      }

      # merge with metadata
      plot <- dplyr::inner_join(genem, obj@metadata |> tibble::rownames_to_column(var = "cell_id"), by = "cell_id")
      # plot
      p[[i]] <- ggplot2::ggplot(plot, ggplot2::aes(x = umap_x, y = umap_y, color = pctm)) +
        ggplot2::geom_point(size = 0.1) + ggplot2::theme_classic() +
        ggplot2::scale_color_gradientn(colors = pal) +
        {if (!is.null(squish)) ggplot2::scale_color_gradientn(colors = pal, limits = c(0, squish), oob = scales::squish)} +
        ggplot2::labs(title = paste0(genes[i])) + noAxes()
    }
    gridExtra::grid.arrange(grobs = p, nrow = nrow)
  }

  if (blend) {
    if (!is.null(matrix)) {
      genem <- obj@genomeMatrices[[matrix]]
      genem <- genem[genes,]
      genem <- as.data.frame(t(genem))
    }
    if (is.null(matrix)) {
      genem <- makeWindows(obj = {{obj}}, genes = {{genes}}, type = {{type}}, metric = {{metric}})
      genem <- as.data.frame(t(genem))
    }
    genem <- genem * 0.01
    if (length(genes) == 2) {
      names(genem) <- c("gene1", "gene2")
      genem <- genem |> dplyr::mutate(mix = grDevices::rgb(red = gene1, green = gene2, blue = 0, maxColorValue = max(genem[1:2])))
      plot <- merge(genem, obj@metadata, by = 0)
      legend <- expand.grid(red = seq(0, max(genem[1:2]), by = 0.02), green = seq(0, max(genem[1:2]), by = 0.02))
      legend <- within(legend, mix <- grDevices::rgb(green = green, red = red, blue = 0, maxColorValue = max(genem[1:2])))
      p1 <- ggplot2::ggplot(plot, ggplot2::aes(x = umap_x, y = umap_y, color = mix)) +
        ggplot2::geom_point(size = 0.5) +
        ggplot2::theme_classic() +
        ggplot2::guides(colour = guide_legend(override.aes = list(size=3))) +
        ggplot2::scale_color_identity() +
        ggplot2::labs(title = paste0(type, " methylation across ", genes[1], " and ", genes[2])) + noAxes()
      p2 <- ggplot2::ggplot(legend, ggplot2::aes(x=red, y=green)) +
        ggplot2::geom_tile(aes(fill=mix), color="white") +
        ggplot2::scale_fill_identity() +
        ggplot2::theme_classic() +
        ggplot2::labs(x = paste(genes[1]), y = paste(genes[2]), title = "Legend")
      p2 <- p2 + coord_fixed(ratio = 1/1)
      plot_grid(p1, p2, rel_widths = c(3,1))
    }
    else if (length(genes) == 3) {
      names(genem) <- c("gene1", "gene2", "gene3")
      genem <- genem |> dplyr::mutate(mix = grDevices::rgb(red = gene1, green = gene2, blue = gene3, maxColorValue = max(genem[1:3])))
      plot <- merge(genem, obj@metadata, by = 0)
      p1 <- ggplot2::ggplot(plot, ggplot2::aes(x = umap_x, y = umap_y, color = mix)) +
        ggplot2::geom_point(size = 0.5) +
        ggplot2::theme_classic() +
        ggplot2::guides(colour = guide_legend(override.aes = list(size=3))) +
        ggplot2::scale_color_identity() +
        ggplot2::labs(title = paste0(type, " methylation across ", genes[1], " (red), ", genes[2], " (green), and ", genes[3], " (blue)")) + noAxes() +
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
#' @title vlnGeneM
#' @description Generates a violin plot of percent methylation over a gene body
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
#' @examples vlnGeneM(obj, genes = c("SATB2", "GAD1"), matrix = "test_gene_ch", groupBy = "cluster_id")
#' @importFrom dplyr inner_join
#' @importFrom ggplot2 ggplot aes geom_violin geom_jitter theme_classic scale_fill_manual scale_color_manual facet_wrap
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
vlnGeneM <-  function(
    obj,
    genes,
    type = NULL,
    metric = "percent",
    matrix = NULL,
    colors = NULL,
    groupBy,
    nrow = round(sqrt(length(genes)))) {

  if (is.null(matrix) && is.null(type)) {
    stop("If matrix is not provided, 'type' must be specified.")
  }

  if (!is.null(matrix)) {
    genem <- obj@genomeMatrices[[matrix]]
    genem <- genem[genes,] |> tibble::rownames_to_column(var = "gene")
    genem <- tidyr::pivot_longer(genem, cols = c(2:ncol(genem)), names_to = "cell_id", values_to = "pctm")
  }
  if (is.null(matrix)) {
    genem <- makeWindows({{obj}}, genes = genes, type = {{type}}, metric = {{metric}}) |> tibble::rownames_to_column(var = "gene")
    genem <- tidyr::pivot_longer(genem, cols = c(2:ncol(genem)), names_to = "cell_id", values_to = "pctm")
  }

  # merge with metadata
  plot <- dplyr::inner_join(genem, obj@metadata |> tibble::rownames_to_column(var = "cell_id"), by = "cell_id")

  ggplot2::ggplot(plot, ggplot2::aes(x = .data[[groupBy]], y = pctm, fill = .data[[groupBy]], color = .data[[groupBy]])) +
    ggplot2::geom_jitter(ggplot2::aes(color = .data[[groupBy]]), size = 0.1, alpha = 0.2) +
    ggplot2::geom_violin(alpha = 0.3) +
    ggplot2::theme_classic() +
    {if (!is.null(colors)) scale_fill_manual(values = colors) } +
    {if (!is.null(colors)) scale_color_manual(values = colors) } +
    ggplot2::facet_wrap(vars(gene), nrow = nrow)
}

############################################################################################################################
#' @title dotGeneM
#' @description Generates a dot plot of percent methylation over a gene body
#'
#' @param obj The amethyst object to plot
#' @param genes  List of genes to plot methylation on the x axis
#' @param type What type of methylation to retrieve; i.e. gene body mCH, gene body mCG, or promoter mCG.
#' @param matrix Optional name of a pre-calculated matrix of values contained in the genomeMatrices slot to use
#' @param metric If matrix is not specified, what methylation parameter to calculate, i.e. "percent", "score", or "ratio"
#' @param groupBy A categorical variable contained in the metadata to group cells by on the y axis
#'
#' @return Returns a ggplot object displaying average methylation values for each gene by the grouping variable
#' @export
#' @examples dotGeneM(obj = obj, genes = c("SATB2", "GAD1"), type = "gene_ch", groupBy = "cluster_id")
#' @importFrom dplyr inner_join summarise group_by
#' @importFrom ggplot2 ggplot aes geom_point theme_classic scale_color_viridis_c theme
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
dotGeneM <-  function(
    obj,
    genes,
    groupBy,
    matrix = NULL,
    type = NULL,
    metric = NULL) {

  if (!is.null(matrix)) {
    genem <- obj@genomeMatrices[[matrix]]
    genem <- genem[genes,] |> tibble::rownames_to_column(var = "gene")
  }
  if (is.null(matrix)) {
    genem <- makeWindows({{obj}}, gene = {{genes}}, type = {{type}}, metric = {{metric}}) |> tibble::rownames_to_column(var = "gene")
  }
  genem <- tidyr::pivot_longer(genem, cols = c(2:ncol(genem)), names_to = "cell_id", values_to = "pctm")
  genem <- dplyr::inner_join(genem, obj@metadata |> tibble::rownames_to_column(var = "cell_id"), by = "cell_id")
  genem <- genem |> dplyr::group_by(.data[[groupBy]], gene) |> dplyr::summarise(pctm = mean(pctm, na.rm = TRUE))

  ggplot2::ggplot(genem, ggplot2::aes(x = gene, y = .data[[groupBy]])) +
    ggplot2::geom_point(aes(size = pctm, color = pctm)) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_viridis_c(name = paste0("pctm", type)) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

}

############################################################################################################################
#' @title histGeneM
#' @description Plot methylation levels over a gene body in histogram format
#' @param obj Amethyst object to plot
#' @param type Type of methylation to plot - either "CG", "CG", or "both"
#' @param groupBy A categorical variable contained in the metadata to facet the plot by
#' @param genes Which genes to plot methylation over
#' @param bins Number of bins to divide the histogram by
#' @param colors Optional color palette to plot
#' @param smooth Boolean indicating whether to plot a smoothed average line of histogram values
#' @param track Boolean indicating whether to plot a gene track below the histogram
#' @param promoter Boolean indicating whether to plot promoter CG levels (appears as a red block)
#' @param track.color Option to change the color in which the gene track is plotted
#' @param alpha Opacity of the histogram
#' @param nrow Number of rows to divide the output across
#' @param legend Boolean indicating whether to plot the groupBy variable legend
#' @return Returns a ggplot histogram showing methylation levels over a gene body
#' @export
#' @examples histGeneM(obj = obj, type = "both", groupBy = cluster_id, genes = "SATB2")
#' @importFrom data.table rbindlist
#' @importFrom dplyr pull filter group_by mutate relocate summarise inner_join arrange select distinct
#' @importFrom ggplot2 ggplot
#' @importFrom plyr round_any
#' @importFrom gridExtra grid.arrange
#' @importFrom tibble rownames_to_column
histGeneM <- function(
    obj,
    type, # "CG" or "CH" methylation, or "both"
    groupBy = NULL, # metadata to divide by
    genes,
    bins = 100, # number of columns to plot across gene body
    colors = makeAmethystPalette(),
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
    stop("No reference found. Make sure obj@ref slot contains the proper annotation file.")
  }
  for (i in 1:length(genes)) {

    # determine size of bins
    binwidth <- (obj@ref |> dplyr::filter(type == "gene" & gene_type == "protein_coding")|>
                   dplyr::mutate(length = end - start) |>
                   dplyr::group_by(gene_name) |>
                   dplyr::mutate(replicate = n()) |>
                   dplyr::group_by(gene_name) |>
                   arrange(desc(length)) |>
                   dplyr::filter(row_number()==1 & gene_name == genes[i]) |>
                   dplyr::pull(length))/bins

    # get average methylation across a gene
    if (type == "CG" | type == "CH") {
      genem <- getGeneM({{obj}}, gene = genes[i], type = {{type}})
    } else if (type == "both") {
      both <- as.list(c("CH", "CG"))
      genem <- map(.x = both, .f = function(x) {getGeneM({{obj}}, gene = genes[i], type = x)})
      names(genem) <- c("CH", "CG")
      genem <- data.table::rbindlist(genem, idcol = "type") |> dplyr::relocate(type, .after = last_col())
    }
    genem <- genem |> dplyr::mutate(start = round_any(pos, binwidth, floor), end = plyr::round_any(pos, binwidth, ceiling), bin = paste0(chr, "_", plyr::round_any(pos, binwidth, floor), "_", plyr::round_any(pos, binwidth, ceiling)))
    if (type == "CG" | type == "CH") {
      genem <- genem |> dplyr::group_by(cell_id, bin) |>
        dplyr::summarise(n = n(), start = mean(start), end = mean(end), middle = mean(c(start, end)), pct_m = (sum(c != 0)*100/(sum(c != 0) + sum(t != 0))), .groups = "keep")
    } else if (type == "both") {
      genem <- genem |> dplyr::group_by(cell_id, bin, type) |>
        dplyr::summarise(n = n(), start = mean(start), end = mean(end), middle = mean(c(start, end)), pct_m = (sum(c != 0)*100/(sum(c != 0) + sum(t != 0))), .groups = "keep")
    }
    # merge with metadata for grouping information
    genem <- dplyr::inner_join(genem, obj@metadata |> tibble::rownames_to_column(var = "cell_id"), by = "cell_id")
    if (type == "CG" | type == "CH") {
      summary <- genem |> dplyr::group_by(start, middle, end, bin, {{groupBy}}) |> dplyr::summarise(mean = mean(pct_m))
      summary$type <- as.character(type)
    } else if (type == "both") {
      summary <- genem |> dplyr::group_by(start, middle, end, bin, {{groupBy}}, type) |> dplyr::summarise(mean = mean(pct_m))
    }
    if (promoter) {
      promoters <- getGeneM({{obj}}, gene = genes[i], type = "promoters")
      promoters <- dplyr::inner_join(promoters, obj@metadata |> tibble::rownames_to_column(var = "cell_id"), by = "cell_id")
      promoters <- promoters |> dplyr::group_by(cell_id, {{groupBy}}) |> dplyr::summarise(n = n(), pct_m = (sum(c != 0)*100/(sum(c != 0) + sum(t != 0))), type = "promoter", .groups = "keep")
      promoters$middle <- obj@ref |> dplyr::filter(type == "gene" & gene_type == "protein_coding") |>
        dplyr::mutate(length = end - start) |> dplyr::group_by(gene_name) |>
        dplyr::mutate(replicate = n()) |>
        group_by(gene_name) |>
        arrange(desc(length)) |>
        dplyr::filter(row_number()==1) |>
        dplyr::mutate(tss = ifelse(strand == "+", start, end)) |>
        dplyr::mutate(promoter_start = tss - 1500, promoter_end = tss + 1500) |>
        dplyr::distinct(seqid, start, end, .keep_all = T) |>
        dplyr::filter(seqid != "chrM") |> dplyr::select("seqid", "promoter_start", "promoter_end", "gene_name", "location") |>
        dplyr::arrange(seqid, promoter_start) |>
        dplyr::filter(gene_name %in% genes[i]) |> dplyr::mutate(middle = mean(c(promoter_start, promoter_end))) |>
        dplyr::pull(middle)
      promoters <- promoters |> dplyr::mutate(start = (middle - 1500), end = (middle + 1500), bin = "promoter") |>
        dplyr::group_by(start, middle, end, {{groupBy}}, type) |>
        dplyr::summarise(mean = mean(pct_m))
      summary <- do.call(rbind, list(summary, promoters))
    }
    if (track) {
      exon <- obj@ref |> dplyr::filter(gene_name == genes[i] & type == "exon") |> dplyr::distinct(seqid, start, end)
    }
    # plot
    p[[i]] <- ggplot2::ggplot() +
      {if(track & promoter)geom_rect(data = promoters, aes(xmin = start, xmax = end, ymin = (0-(max(summary$mean)*.06)), ymax = (0-(max(summary$mean)*.03))), fill = "red", alpha = .4)} +
      {if(track)geom_rect(data = exon, aes(xmin = start, xmax = end, ymin = (0-(max(summary$mean)*.06)), ymax = (0-(max(summary$mean)*.03))), color = track.color, fill = track.color, alpha = 1)} + # if exons is true, plot exons track
      {if(track)geom_rect(data = exon, aes(xmin = min(start), xmax = max(end), ymin = (0-(max(summary$mean)*.045)), ymax = (0-(max(summary$mean)*.045))), color = track.color, alpha = 1)} + # also plot light line for intron track
      {if(type == "CG" | type == "CH")geom_rect(data = summary %>% dplyr::filter(type != "promoter"), aes(xmin = start, xmax = end, ymin = 0, ymax = mean, fill = {{groupBy}}), alpha = 0.4)} +
      {if(type == "both")geom_rect(data = summary %>% dplyr::filter(type != "promoter"), aes(xmin = start, xmax = end, ymin = 0, ymax = mean, fill = {{groupBy}}, alpha = type))} +
      {if(promoter)geom_rect(data = summary %>% dplyr::filter(type == "promoter"), aes(xmin = start, xmax = end, ymin = 0, ymax = mean), fill = "red", alpha = 0.2)} +
      {if(smooth & type == "CG")geom_smooth(data = genem, aes(x = middle, y= pct_m, color = {{groupBy}}))} +      # smoothed lines of mean % CH of exons and introns
      {if(smooth & type == "CH")geom_smooth(data = genem, aes(x = middle, y= pct_m, color = {{groupBy}}))} +
      {if(smooth & type == "both")geom_smooth(data = genem, aes(x = middle, y= pct_m, color = {{groupBy}}, alpha = type))} +
      facet_wrap(vars({{groupBy}}), nrow = 1) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + # wrap by cell class; make graph look nice
      scale_fill_manual(values = colors) + scale_color_manual(values = colors) + # colors
      {if(!legend)theme(legend.position = "none")} + scale_alpha_manual(values = c(.2, .6)) +
      {if(type != "both")ggtitle(paste0({{type}}, " methylation along ", genes[i], " body"))} +
      {if(type == "both")ggtitle(paste0("CG+CH methylation along ", genes[i], " body"))} +
      labs(x = paste0("Location along ", unique(separate(data = genem, col = bin, into = "chr", sep = "_", extra = "drop")$chr)))
  }
  gridExtra::grid.arrange(grobs = p, nrow = nrow)
}

############################################################################################################################
#' @title tileGeneM
#' @description Plot CG methylation levels over a gene body in 500bp smoothed windows
#'
#' @param obj Object containing the matrix to plot
#' @param genes Gene list to plot
#' @param matrix Matrix contained in the genomeMatrices slot. Note: unlike other functions, designed for makeSlidingWindows output
#' @param colors Optional list of colors to include
#' @param trackOverhang Number of base pairs to extend beyond the gene
#' @param arrowOverhang Number of base pairs the track arrow should extend beyond the gene
#' @return A ggplot geom_tile object with colors indicating % methylation over 500 bp windows and the gene of interest beneath
#' @export
#' @examples tileGeneM(obj, gene = "SYT7", matrix = "cluster_cg_500_slidingwindows")
#' @importFrom dplyr filter mutate
#' @importFrom ggplot2 ggplot geom_tile aes geom_rect geom_segment theme scale_fill_gradientn ylab
#' @importFrom gridExtra grid.arrange
#' @importFrom tidyr pivot_longer
tileGeneM <- function(obj,
                      genes,
                      matrix,
                      colors = NULL,
                      trackOverhang = 5000,
                      arrowOverhang = 3000,
                      nrow = length(genes),
                      legend = TRUE,
                      width = 500) {
  # make empty plot list
  p <- vector("list", length(genes)) # empty plot list
  for (i in 1:length(genes)) {

    ref <- obj@ref |> dplyr::filter(gene_name == genes[i])
    aggregated <- obj@genomeMatrices[[matrix]]

    if (!is.null(colors)) {
      colors <- colors
    } else {
      colors <- c("#FF0082", "#dbdbdb", "#cccccc", "#999999")
    }

    toplot <- aggregated[c((aggregated$chr == ref$seqid[ref$type == "gene"] & aggregated$start > (ref$start[ref$type == "gene"] - trackOverhang) & aggregated$end < (ref$end[ref$type == "gene"] + trackOverhang))), ]
    ngroups <- ncol(toplot) - 3
    trackHeight <- ngroups * .07
    toplot <- tidyr::pivot_longer(toplot, cols = c(4:ncol(toplot)), names_to = "cluster_id", values_to = "pct_mCG") |> rowwise() |> dplyr::mutate(middle = mean(c(start, end), na.rm = TRUE))

    p[[i]] <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = toplot, ggplot2::aes(x = middle, y = cluster_id, fill = pct_mCG), width = width) +
      ggplot2::geom_rect(fill = "pink", data = ref |> dplyr::filter(type == "gene") |> dplyr::mutate(promoter_start = ifelse(strand == "+", (start - 1500), (end + 1500)), promoter_end = ifelse(strand == "+", (promoter_start+3000), (promoter_start-3000))),
                         ggplot2::aes(xmin = promoter_start, xmax = promoter_end, ymin = -trackHeight, ymax = 0)) +
      ggplot2::geom_rect(fill = "black", data = ref |> dplyr::filter(type == "exon"), ggplot2::aes(xmin = start, xmax = end, ymin = -trackHeight, ymax = 0)) +
      ggplot2::geom_segment(data = ref, aes(x = ifelse(strand == "+", (min(start) - arrowOverhang), (max(end)) + arrowOverhang), y = -(trackHeight/2), xend = ifelse(strand == "+", (max(end) + arrowOverhang), (min(start)) - arrowOverhang), yend = - (trackHeight/2)), arrow = arrow(length = unit(trackHeight/4, "cm"))) + xlab(genes[i]) +
      ggplot2::geom_rect(fill = "black", data = ref |> dplyr::filter(type == "gene"), ggplot2::aes(xmin = start, xmax = end, ymin = (-(trackHeight/2) -(trackHeight/20)), ymax = (-(trackHeight/2) + (trackHeight/20)))) +
      {if (!legend)ggplot2::theme(legend.position = "none")} +
      ggplot2::scale_fill_gradientn(colors = colors) +
      ggplot2::theme(panel.background = element_blank(), axis.ticks = element_blank()) + ggplot2::ylab("Cluster ID")
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
  if (option == 2) {pal <- c("#233A32", "#5b8e7d","#8cb369","#f4e285","#f4a259","#bc4b51")}
  if (option == 3) {pal <- c("#150808","#ec5151","#ffdd6d","#a8dadc")}
  if (option == 4) {pal <- c("#293241","#3d5a80","#98c1d9","#e0fbfc","#ee6c4d")}
  if (option == 5) {pal <- c("#150808","#ec5151","#ffdd6d","#a8dadc")}
  if (option == 6) {pal <- c("#0D3B66","#FAF0CA","#F4D35E","#EE964B","#F95738")}
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
  if (option == 17) {pal <- c("#227c9d","#17c3b2","#ffcb77","#fef9ef","#fe6d73")}
  if (option == 18) {pal <- c("#f7ede2","#f6bd60","#f5cac3","#84a59d","#f28482")}
  if (option == 19) {pal <- c("#218380","#73d2de","#ffbc42","#d81159","#8f2d56")}
  if (option == 20) {pal <- c("#562c2c","#f2542d","#f5dfbb","#0e9594","#127475")}
  if (option == 21) {pal <- c("#73b7b8","#52a1a3","#76c8b1","#50b99b","#f6cb52","#f3b816","#d23f0f","#f05a29","#dc244b","#af1d3c")}
  if (option == 22) {pal <- c("#F05252","#FFA976","#66C2A5","#aae3e3","#004A4A","#FF8B73","#A6D854","#419e9e","#FFD92F","#bd083e", "#fcba2b", "#8DA0CB","#75C8D2","#E78AC3","#B2DF8A","#1F78B4","#E5C494","#B3B3B3","#FC8D62","#FB9A99")}
  if (option == 23) {pal <- c("#50514f","#f25f5c","#ffe066","#247ba0","#70c1b3","#c0e8f9")}
  if (option == 24) {pal <- c("#023047","#219ebc","#8ecae6","#ffb703","#fb8500")}
  if (option == 25) {pal <- c("#885053","#fe5f55","#777da7","#94c9a9","#c6ecae")}
  if (option == 26) {pal <- c("#1c5253", "#306b34","#c3eb78", "#f3ffc6","#b6174b")}
  if (option == 27) {pal <- c("#b6e2dd","#c8ddbb","#e9e5af","#fbdf9d","#fbc99d","#fbb39d","#fba09d")}
  if (option == 28) {pal <- c("#247ba0","#a15856","#f25f5c","#f9a061","#ffe066","#92ae83","#70c1b3","#4a9eaa","#50514f")}
  if (option == 29) {pal <- c("#84a59d","#f6bd60","#f7ede2","#f5cac3","#f28482")}
  if (option == 30) {pal <- c("#dad2d8","#143642","#0f8b8d","#a8201a","#ec9a29")}

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
#' or "umap", which shows your amethyst object colored by cluster_id in each palette option.
#' @param n If output = "swatch", how many colors should be shown with each palette option.
#' @param obj If output = "umap", name of the amethyst object to test.
#' @return If output = "swatch", returns a ggplot object with each row as a palette.
#' If output = "umap", returns the amethyst object in each palette option.
#' @export
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_tile theme_void scale_fill_identity facet_wrap
#' @examples testPalette(output = "swatch", n = 30)
#' @examples testPalette(output = "umap", obj = obj)
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
  } else if (output == "umap") {
    p <- vector("list", 30) # empty plot list
    for (i in 1:30) {
      colors <- sample(makePalette(option = i, n = length(unique(obj@metadata$cluster_id))))
      # get average methylation across a gene
      p[[i]] <- umapFeature(obj = obj, colorBy = cluster_id, pointSize = 0.1, colors = colors) +
        theme(legend.position = "none") + ggtitle(paste0("pal ", i))
    }
    gridExtra::grid.arrange(grobs = p, nrow = 5)
  }
}
