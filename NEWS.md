# v1.0.3
### 08-21-2025
  - **Fixed**
    - *loadWindows* issue thanks to @shengzha
    - *filterDMR* and *collapseDMR*; test names now added as "member_id" column instead of order
    - Fixed *heatMap/histograM* frameshifting issue when duplicate genes are used as input
    
  - **Added**
    - More threshold parameters to *aggregateMatrix* (minCells, minValues)
    - *addPrefix* command suggested by @hug-cr
    - Scripts corresponding to upcoming manuscript release in *Communications Biology*
    - Integrated error message when obj@h5paths isn't correctly structured
    - Color palette names! With help from Marissa Co

# v1.0.2
### 07-24-2025
  - **Fixed**
    - *loadSmoothedWindows* issue thanks to @shengzha
    - *createScaleObject* issue (hopefully) thanks @yosefellenbogen1 @edogiuili
    - *testDMR* bug; huge thanks to @stcolema. We would recommend re-running testDMR result from previous versions.
    - *filterDMR* and *collapseDMR*; test names now added as "member_id" column instead of order
    
# v1.0.1
### 07-07-2025
  - **Fixed**
    - prefix integration issues affecting *indexChr*, *makeWindows*, and *calcSmoothedWindows*
    - *createObject* locked error (thanks to @yosefellenbogen1)
    - minor updates throughout the vignettes (thanks to Jack Henry Kotnik)
    
  - **Added**
    - started tracking version releases

# v1.0.0
### 06-02-2025 - v1.0.0 release!
  - **Added**
    - Enabled processing of datasets with overlapping barcodes
      - Necessitated change in obj@h5path structure
      - Necessitated change in functions: *makeWindows, indexChr, calcSmoothedWindows*
      - [Vignette for combining datasets with overlapping barcodes](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/combining_overlapping_barcodes/combining_overlapping_barcodes.html) 
    - [Vignette for alternative clustering approaches](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/clustering_alternatives/clustering_alternatives.html)
    - facet.R script containing functions to help load Facet results
      - Affected functions: *loadWindows, loadSmoothedWindows*
    - Enabled chromosome whitelist functionality
      - Affected functions: *makeWindows, indexChr, calcSmoothedWindows*
    - *getGeneCoords* helper function
    - *heatMapGenome* function

  - **Changed**
    - Updated the [PBMC vignette](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/pbmc_vignette/pbmc_vignette.html) to reflect v1.0.0 changes
    - Updated the [brain vignette](http://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/brain_vignette/brain_vignette.html) to reflect v1.0.0 changes
    - Changed h5 path input structure
      - base information in /context/barcode/1
      - aggregated information in /context/barcode/name
    - Changed amethyst object structure
      - **h5paths** slot now has "barcode" and "prefix" (optional) columns
      - added **tracks** and **results** slots
        - Affected functions: *createObject, convertObject*
      - changed **obj@ref** slot to data.table
      - moved UMAP and TSNE results to **reductions** slot
        - Affected functions: *runUmap, runTsne, dimM, dimFeature*
      - enabled multiple clustering results to be stored in metadata
        - Affected functions: *runCluster*
      - provided *convertObject* function to automatically initiate these changes
      - *heatMap* and *histograM* can plot tracks underneath
      - *dotM* can accept two matrices for size aesthetics instead of one

  - **Deprecated**
    - Moved *clusterCompare, findRanges, indexGenes,* and *getGeneM* to deprecated
    
  - **Fixed**
    - *makeWindows* promoters can now also calculate "score" and "ratio"
    - Specified DESCRIPTION imports
    - Added gene bar across *heatMap* and *histograM* regions where result spans the whole window
    - Fixed minor *heatMap* and *histograM* issues

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
        - histograM baseline can either be 0 or mean methylation (thanks to Ryan Mulqueen)
        - heatMap color scale and max value can be adjusted
        - dotM can be further faceted by a variable in the metadata
  - **Changed** *findClusterMarkers* to logFC = log2(mean_1 / mean_2) (thanks to Joe Verity-Legg)
	  - Affected functions: *findClusterMarkers*

### 08-27-2024
  - **Added** compatibility for more genome builds 
	  - Affected functions: *makeRef, fetchMarkers, makeWindows, calcSmoothedWindows*
