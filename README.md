# Amethyst: A METHyl-c Single-cell analysis Toolkit

<!-- badges: start -->
<!-- badges: end -->

<p align="center">
  <img src="https://github.com/lrylaarsdam/amethyst/blob/main/images/amethyst.png?raw=true" alt="Amethyst" width="250" />
</p>

Single-cell sequencing technologies have revolutionized biomedical research by enabling deconvolution of cell type-specific properties in highly heterogeneous tissue. While robust tools have been developed to handle bioinformatic challenges posed by single-cell RNA and ATAC data, options for emergent modalities such methylation are much more limited, impeding the utility of results. Here we present Amethyst, the first comprehensive R package for atlas-scale single-cell methylation sequencing data analysis. Amethyst takes base-level methylation calls and facilitates batch integration, doublet detection, dimensionality reduction, clustering, cell type annotation, differentially methylated region calling, and interpretation of results all in one streamlined platform. Versatile visualization functions mediate rapid interaction with the data in a local environment. Efforts like Amethyst will increase accessibility to single-cell methylation data interpretation, accelerating progress in understanding principles of this critical epigenetic modification across diverse contexts.  

## Installation

Installation of Amethyst can be done using devtools:

```{r}
library(devtools)
devtools::install_github("lrylaarsdam/amethyst")
```

You will likely need to install one or more dependencies:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install()

library(BiocManager)
BiocManager::install(c("caret", "devtools", "data.table", "dplyr", "furrr", "future", "future.apply",
  "ggplot2", "grDevices", "gridExtra", "igraph", "irlba", "janitor", "Matrix", "methods", 
  "plotly", "plyr", "purrr", "randomForest", "rhdf5", "rtracklayer", "scales", "stats", "stringr", 
  "tibble", "tidyr", "umap", "utils"))

devtools::install_github("JinmiaoChenLab/Rphenograph")
devtools::install_github("KrishnaswamyLab/MAGIC/Rmagic")
```

## Getting Started

Amethyst begins with base-level methylation calls per cell wrapped into .h5 files. If you need to generate this file from your sequencing data, scripts for initial processing of reads are available at the Adey Lab [Premethyst](https://github.com/adeylab/premethyst) repo. Please see [vignettes](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/pbmc_vignette/pbmc_vignette.html) for example Premethyst outputs and subsequent analysis steps.

If you are using the Scale Biosciences pipeline, we have written a helper function to generate an Amethyst object from the output. createScaleObject automatically populates the metadata and h5path slots for you. In its most basic form, all that is needed is the directory path:

```{r}
obj <- createScaleObject(directory = "path/to/scalebio/output/folder")
```

You may also want to load the automatically generated matrices. Below is an example of how to load the "CG.score" matrix. However, please make sure your computational resources are equipped to handle the entire matrix size first.

```{r}
obj <- createScaleObject(directory = "path/to/scalebio/output/folder", genomeMatrices = c("CG.score"))
```

## Vignettes

To become familiar with the Amethyst workflow, we recommend beginning with the [pbmc vignette](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/pbmc_vignette/pbmc_vignette.html), which is focused on CG methylation and applicable to any tissue. 

Certain tissues - such as the brain and stem cells - also contain high levels of non-CG methylation and necessitate a very different workflow. After completing the pbmc vignette, we recommend going over the [brain vignette](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/brain_vignette/brain_vignette.html) for CH-specific analysis.

Additional vignettes for doublet detection, batch integration, subsetting, and other utilities will be coming shortly.

## Issues

Amethyst is still a work in progress. Please let us know if any [issues](https://github.com/lrylaarsdam/amethyst/issues) come up. 

## License

Amethyst is distributed under the MIT License. Please see LICENSE.txt for further information. 

