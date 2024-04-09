##  amethyst -- R package for visualisation of single cell methylation data
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Attaching amethyst version ",
                        packageDescription("amethyst")$Version, ".")
  options(scipen = 999)
  options(future.globals.maxSize = (4000*1024^2)) # 4 GB
}
