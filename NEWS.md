# v1.0.0 (dev)

# v0.0.0.9000
### 11-14-2024
  - **Added** vignettes
    - [Doublet detection](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/doublet_detection/doublet_detection.html)
    - [Batch integration](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/batch_correction/batch_correction.html)
    - [Additional utilities: subsetting, merging, imputation](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/additional_utilities/additional_utilities.html)

### 11-01-2024
  - **Added** flexibility with visualization parameters
	  - Affected functions: *histograM, heatMap, dotM*
    - Examples:
        - histograM baseline can either be 0 or mean methylation (credit: Ryan Mulqueen, PhD)
        - heatMap color scale and max value can be adjusted
        - dotM can be further faceted by a variable in the metadata
  - **Changed** *findClusterMarkers* to logFC = log2(mean_1 / mean_2) (credit: Joe Verity-Legg)
	  - Affected functions: *findClusterMarkers*

### 08-27-2024
  - **Added** compatibility for more genome builds 
	  - Affected functions: *makeRef, fetchMarkers, makeWindows, calcSmoothedWindows*
