#' @importFrom cowplot plot_grid
#' @importFrom data.table as.data.table %between% %inrange% fread frollmean frollsum frollapply rbindlist setDT tstrsplit
#' @import dplyr
#' @importFrom furrr future_map
#' @import future
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom irlba irlba
#' @importFrom janitor make_clean_names clean_names
#' @import methods
#' @import patchwork
#' @importFrom pheatmap pheatmap
#' @importFrom plotly layout plot_ly subplot
#' @importFrom plyr round_any
#' @importFrom purrr map
#' @importFrom Rmagic magic
#' @importFrom Rphenograph Rphenograph
#' @importFrom rhdf5 h5read
#' @importFrom rtracklayer readGFF
#' @importFrom stats cor
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom umap umap

pal <- c("#004A4A", "#F05252", "#419e9e", "#fcba2b", "#bd083e", "#FB9A99", "#75C8D2",  "#FF8B73", "#B2DF8A", "#1F78B4", "#E31A1C",  "#aae3e3",  "#FFA976")
options(scipen = 999)
options(future.globals.maxSize = (4000*1024^2)) # 4 GB
