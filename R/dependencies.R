packages <- c("rhdf5", "rtracklayer", "ggnewscale", "Rphenograph", "irlba",
              "Rcpp", "umap", "ggplot2", "pheatmap", "gridExtra", "grid",
              "patchwork", "lattice", "cowplot", "plotly", "janitor",
              "data.table", "readr", "tidyverse", "tibble", "plyr",
              "tidyr", "Rmagic", "future", "furrr", "purrr", "dplyr")

# Load packages with tryCatch
for (pkg in packages) {
  tryCatch(library(pkg, character.only = TRUE),
           error = function(e) {
             message(paste("Error loading package", pkg, ":", conditionMessage(e)))
           }
  )
}

pal <- c("#004A4A", "#F05252", "#419e9e", "#fcba2b", "#bd083e", "#FB9A99", "#75C8D2",  "#FF8B73", "#B2DF8A", "#1F78B4", "#E31A1C",  "#aae3e3",  "#FFA976")
options(scipen = 999)
options(future.globals.maxSize = (20000*1024^2)) # 20 GB; adjust if local
