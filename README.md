# Amethyst

*A comprehensive toolkit for single-cell methylation sequencing data analysis*

<!-- badges: start -->
<!-- badges: end -->

<p align="center">
<img src="https://github.com/lrylaarsdam/amethyst/blob/main/images/amethyst.png?raw=true" alt="Amethyst" width="250"/>
</p>

## Welcome to Amethyst!

Single-cell sequencing technologies have revolutionized biomedical research by enabling deconvolution of cell type-specific properties in highly heterogeneous tissue. While robust tools have been developed to handle bioinformatic challenges posed by single-cell RNA and ATAC data, options for emergent modalities such methylation are much more limited, impeding the utility of results. Here we present Amethyst, the first comprehensive R package for atlas-scale single-cell methylation sequencing data analysis. Amethyst takes base-level methylation calls and facilitates batch integration, doublet detection, dimensionality reduction, clustering, cell type annotation, differentially methylated region calling, and interpretation of results all in one streamlined platform. See our [manuscript](https://www.nature.com/articles/s42003-025-08859-2) to learn more! 


## Vignettes

To become familiar with the Amethyst workflow, we recommend beginning with the [PBMC vignette](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/pbmc_vignette/pbmc_vignette.html), which is focused on CG methylation analysis and applicable to any tissue.

Certain tissues - such as the brain and stem cells - also contain high levels of non-CG methylation and necessitate a very different analysis approach. After completing the PBMC vignette, we recommend going over the [brain vignette](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/brain_vignette/brain_vignette.html) for mCH-specific analysis.

In addition to these general workflow examples, we have specific vignettes for: 
  - [Combining datasets with overlapping barcodes](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/combining_overlapping_barcodes/combining_overlapping_barcodes.html) 
  - [Batch integration](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/batch_correction/batch_correction.html) 
  - [Alternative clustering approaches](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/clustering_alternatives/clustering_alternatives.html)
  - [Doublet detection](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/doublet_detection/doublet_detection.html) 
  - [Additional utilities: subsetting, merging, imputation](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/additional_utilities/additional_utilities.html)


## Installation

You will likely need to install one or more dependencies first:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

library("BiocManager")
BiocManager::install(c("caret", "devtools", "data.table", "dplyr", "furrr", "future", "future.apply",
  "ggplot2", "grDevices", "gridExtra", "igraph", "irlba", "janitor", "Matrix", "methods", "pheatmap",
  "plotly", "plyr", "purrr", "randomForest", "rhdf5", "rtracklayer", "Rtsne", "scales", "stats", "stringr", 
  "tibble", "tidyr", "umap", "utils"))

devtools::install_github("JinmiaoChenLab/Rphenograph")
devtools::install_github("KrishnaswamyLab/MAGIC/Rmagic")
devtools::install_github("TomKellyGenetics/leiden")
```

Installation of Amethyst can then be done using devtools or remotes:

```{r}
library("devtools")
devtools::install_github("lrylaarsdam/amethyst")

library("remotes")
remotes::install_github("lrylaarsdam/amethyst")
```


## Input data structures

Amethyst begins with base-level methylation calls per cell wrapped into h5 files. The structure of the h5 file is illustrated in the diagram below. If desired, aggregate methylation levels over features can be calculated with [Facet](https://pypi.org/project/amethyst-facet/) and stored in the h5 file as well.

<p align="center">
<img src="https://github.com/lrylaarsdam/amethyst/blob/main/images/h5structure.png?raw=true" alt="Amethyst" width="550"/>
</p>

Scripts for initial processing of sequencing data to produce this input format are available at the Adey Lab [Premethyst](https://github.com/adeylab/premethyst) repo. Please see [vignettes](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/pbmc_vignette/pbmc_vignette.html) for example Premethyst outputs and subsequent analysis steps.

### Converting the h5 object from v0.0.0.9000 to v1.0+

We recently made many improvements to Amethyst that required upgrades to the input data structures. If you have h5 files in the old v0.0.0.9000 format, [Facet](https://pypi.org/project/amethyst-facet/) has a helper function to convert the files so base-level observations are stored under context/barcode/1.

```{bash}
facet convert new_format.h5 old_format.h5
```

If you are using the Scale Biosciences pipeline, we have written a helper function to load the output into an Amethyst object. *createScaleObject* automatically populates the metadata and path slots for you. In its most basic form, all that is needed is the directory path:

```{r}
obj <- createScaleObject(directory = "~/Downloads/ScaleMethyl.out/samples")
```

You may also wish to load any pre-generated matrices, which would allow one to skip past the *makeWindows* step in the vignette. Below is an example of how to load the "CG.score" matrix. Double-check your computational resources are capable of handling the entire matrix size first.

```{r}
obj <- createScaleObject(directory = "~/Downloads/ScaleMethyl.out/samples", genomeMatrices = list("CG.score", "CH"))
```

If using neither Premethyst nor ScaleMethyl, any pre-processing platform can be used. The important thing is that the data is in the expected structure. The diagram below illustrates each slot and its respective data contents.

<p align="center">
<img src="https://github.com/lrylaarsdam/amethyst/blob/main/images/objectstructure.png?raw=true" alt="Amethyst" width="750"/>
</p>

### Converting the Amethyst object structure from v0.0.0.9000 to v1.0+

You might notice in the diagram above that we implemented some minor changes to the object structure in v1.0. To make this transition as smooth as possible, we have provided a helper function *convertObject* to convert format v0.0.0.9000 to v1.0.

```{r}
new_obj <- convertObject(obj = old_obj)
```


## Issues

Amethyst is still a work in progress. Please let us know if any [issues](https://github.com/lrylaarsdam/amethyst/issues) come up.


## News

Our [manuscript](https://www.nature.com/articles/s42003-025-08859-2) is now officially published in *Communications Biology!* Please cite if you use Amethyst for your analysis.

We have recently implemented widespread improvements from v0.0.0.9000. Changes include:
  - Base-resolution methylation information is now expected to be stored in the h5 file under context/barcode/1 (see **Getting started**)
  - Aggregated values can also be calculated and stored in the h5 file under context/barcode/name using the efficient helper package [Facet](https://pypi.org/project/amethyst-facet/)
  - New slots added to store methylation tracks and results (see **Getting started**)
  - Multiple clustering and dimensionality reductions can be stored within the same object 
  - Chromosome white lists can be accommodated 
  - Experiments with overlapping barcode sets can be combined (see [vignette](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/combining_overlapping_barcodes/combining_overlapping_barcodes.html))

Please see our [news](https://github.com/lrylaarsdam/amethyst/blob/HEAD/NEWS.md) for a more comprehensive change log and version update information. Detailed object structure explanations and conversion instructions are included above. 


## License

Amethyst is distributed under the MIT License. Please see LICENSE.txt for further information.
