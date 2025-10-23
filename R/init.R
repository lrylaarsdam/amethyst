# Authors: Stephen Coleman, PhD; Lauren Rylaarsdam, PhD
# 2024-2025
############################################################################################################################
##  amethyst -- R package for single cell methylation data analysis and visualization
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Attaching amethyst version ",
                        packageDescription("amethyst")$Version, ".")
  options(scipen = 999)
  options(future.globals.maxSize = (4000*1024^2)) # 4 GB
  cat("Welcome to Amethyst! If you use Amethyst in your work, please cite our associated manuscript:
      Rylaarsdam, L.E. et al. Single-cell DNA methylation analysis tool Amethyst resolves distinct non-CG methylation patterns in human astrocytes and oligodendrocytes.
      Commun. Biol. 8, 1468 (2025). https://www.nature.com/articles/s42003-025-08859-2")
}
