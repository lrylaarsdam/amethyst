# Amethyst
*A comprehensive toolkit for single-cell methylation sequencing data analysis*

<!-- badges: start -->
<!-- badges: end -->


<p align="center">
  <img src="https://github.com/lrylaarsdam/amethyst/blob/main/images/amethyst.png?raw=true" alt="Amethyst" width="250" />
</p>

Single-cell sequencing technologies have revolutionized biomedical research by enabling deconvolution of cell type-specific properties in highly heterogeneous tissue. While robust tools have been developed to handle bioinformatic challenges posed by single-cell RNA and ATAC data, options for emergent modalities such methylation are much more limited, impeding the utility of results. Here we present Amethyst, the first comprehensive R package for atlas-scale single-cell methylation sequencing data analysis. Amethyst takes base-level methylation calls and facilitates batch integration, doublet detection, dimensionality reduction, clustering, cell type annotation, differentially methylated region calling, and interpretation of results all in one streamlined platform. Versatile visualization functions mediate rapid interaction with the data in a local environment. Efforts like Amethyst will increase accessibility to single-cell methylation data interpretation, accelerating progress in understanding principles of this critical epigenetic modification across diverse contexts. To learn more, see our [preprint](https://www.biorxiv.org/content/10.1101/2024.08.13.607670v2.full.pdf+html)!  


## Dev branch notes

We are currently implementing widespread improvements to Amethyst. The biggest change is base-resolution methylation information is now expected to be stored in the h5 file under context/barcode/1. This structure will allow aggregated metrics to be stored in the same file: for example, methylation values over 100kb windows could be stored in context/barcode/100000. This change will 1) eliminate the need to store massive genomeMatrices within the amethyst object, and 2) allow the most computationally-heavy steps to be calculated outside of R. Users can still do all steps in R if desired, but this method will help accommodate efficient processing of very large datasets.

In addition: 
  - new slots are being added to store methylation tracks and results
  - multiple clustering and dimensionality reductions can be stored within the same object
  - chromosome white list can be accommodated (to do)
  - experiments with overlapping barcode sets can be combined (to do)

We appreciate your patience as these changes are being incorporated.

## Installation

You will likely need to install one or more dependencies first:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install()

library(BiocManager)
BiocManager::install(c("caret", "devtools", "data.table", "dplyr", "furrr", "future", "future.apply",
  "ggplot2", "grDevices", "gridExtra", "igraph", "irlba", "janitor", "Matrix", "methods", "pheatmap",
  "plotly", "plyr", "purrr", "randomForest", "rhdf5", "rtracklayer", "Rtsne", "scales", "stats", "stringr", 
  "tibble", "tidyr", "umap", "utils"))

devtools::install_github("JinmiaoChenLab/Rphenograph")
devtools::install_github("KrishnaswamyLab/MAGIC/Rmagic")
```

Installation of Amethyst can then be done using devtools:

```{r}
library(devtools)
devtools::install_github("lrylaarsdam/amethyst")
```

If that doesn't work, try installing with remotes:

```{r}
library("remotes")
remotes::install_github("lrylaarsdam/amethyst")
```

## Getting started

Amethyst begins with base-level methylation calls per cell wrapped into h5 files. The structure of the h5 file is illustrated in the diagram below. If desired, aggregate methylation levels over features can be calculated with [facet](https://pypi.org/project/amethyst-facet/) and stored in this h5 file as well. 

<p align="center">
  <img src="https://github.com/lrylaarsdam/amethyst/blob/dev/images/h5structure.png?raw=true" alt="Amethyst" width="650" />
</p>

If you need to generate this file from your sequencing data, scripts for initial processing of reads are available at the Adey Lab [Premethyst](https://github.com/adeylab/premethyst) repo. Please see [vignettes](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/pbmc_vignette/pbmc_vignette.html) for example Premethyst outputs and subsequent analysis steps.

If you are using the Scale Biosciences pipeline, we have written a helper function to load the output into an Amethyst object. createScaleObject automatically populates the metadata and h5path slots for you. In its most basic form, all that is needed is the directory path:

```{r}
obj <- createScaleObject(directory = "~/Downloads/ScaleMethyl.out/samples")
```

You may also wish to load any pre-generated matrices, which would allow one to skip past the makeWindows step in the vignette. Below is an example of how to load the "CG.score" matrix. Double-check your computational resources are capable of handling the entire matrix size first.

```{r}
obj <- createScaleObject(directory = "~/Downloads/ScaleMethyl.out/samples", genomeMatrices = list("CG.score", "CH"))
```


## Vignettes

To become familiar with the Amethyst workflow, we recommend beginning with the [pbmc vignette](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/pbmc_vignette/pbmc_vignette.html), which is focused on CG methylation and applicable to any tissue. 

Certain tissues - such as the brain and stem cells - also contain high levels of non-CG methylation and necessitate a very different workflow. After completing the pbmc vignette, we recommend going over the [brain vignette](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/brain_vignette/brain_vignette.html) for CH-specific analysis.

In addition to these general workflow examples, we have specific vignettes for:
- [Doublet detection](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/doublet_detection/doublet_detection.html)
- [Batch integration](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/batch_correction/batch_correction.html)
- [Additional utilities: subsetting, merging, imputation](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/additional_utilities/additional_utilities.html)


## Issues

Amethyst is still a work in progress. Please let us know if any [issues](https://github.com/lrylaarsdam/amethyst/issues) come up. 


## News

Please see our [news](https://github.com/lrylaarsdam/amethyst/blob/HEAD/NEWS.md) for a change log and version update information. Widespread improvements are currently being implemented on the development branch.


## License

Amethyst is distributed under the MIT License. Please see LICENSE.txt for further information. 


