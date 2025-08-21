# Authors: Stephen Coleman, PhD; Lauren Rylaarsdam, PhD
# 2024-2025
############################################################################################################################
##  amethyst -- R package for single cell methylation data analysis and visualization
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Attaching amethyst version ",
                        packageDescription("amethyst")$Version, ".")
  options(scipen = 999)
  options(future.globals.maxSize = (4000*1024^2)) # 4 GB
}
