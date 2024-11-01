# Amethyst Version: 0.0.0.9000

### Added

### Changed
- Nov 1 2024: Added flexibility with visualization parameters
  Examples: histograM baseline can either be 0 or mean methylation (credit: Ryan Mulqueen, PhD)
	    heatMap color scale and max value can be adjusted
	    dotM can be further faceted by a variable in the metadata
  Affected functions: histograM, heatMap, dotM

### Deprecated

### Removed

### Fixed
- Aug 27 2024: Added compatibility for more genome builds. 
  Affected functions: makeRef, fetchMarkers, makeWindows, calcSmoothedWindows
- Nov 1 2024: findClusterMarkers switched to logFC = log2(mean_1 / mean_2) (credit: Joe Verity-Legg)

### Security
