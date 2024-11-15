# Amethyst Version: 0.0.0.9000

### Added
- Nov 14: Vignettes added
  - Doublet detection
  - Batch integration
  - Additional utilities: subsetting, merging, imputation

### Changed
- Nov 1 2024: Added flexibility with visualization parameters
	- Affected functions: histograM, heatMap, dotM
  	- Examples:
		- histograM baseline can either be 0 or mean methylation (credit: Ryan Mulqueen, PhD)
		- heatMap color scale and max value can be adjusted
		- dotM can be further faceted by a variable in the metadata


### Deprecated

### Removed

### Fixed
- Aug 27 2024: Added compatibility for more genome builds
	- Affected functions: makeRef, fetchMarkers, makeWindows, calcSmoothedWindows
- Nov 1 2024: switched to logFC = log2(mean_1 / mean_2) (credit: Joe Verity-Legg)
	- Affected functions: findClusterMarkers

### Security

