# monocle

# monocle 2.9.2

* **bugs**:
  * Replaced deprecated igraph APIs used by trajectory ordering and plotting, including `dfs(father = TRUE)`, `shortest.paths()`, `graph.empty()`, and `get.edgelist()`.
  * Fixed `estimateDispersions()` argument forwarding so dispersion fitting warnings are suppressed by default unless `verbose = TRUE`.
  * Replaced remaining `ggplot2::aes_string()` usages in plotting code with tidy-evaluation mappings.
  * Updated tests to use current message expectations, explicit test data loading, and namespaced `parallel::detectCores()`.

# monocle 2.9.0

* **bugs**:
  * Fixed roxygen2 documentation generation errors: added UTF-8 encoding declaration, corrected function imports from BiocGenerics to parallel package, and fixed special character encoding issues.
  * Fixed `class()` string comparison issues by replacing with `inherits()` for better type checking in multiple functions.
  * Fixed missing imports: added `quantile` from stats package and `read.delim` from utils package.
  * Fixed Rd cross-reference for `densityClust` package.
  * `exportCDS()`: Added check for deprecated `scater::newSCESet` function with clear error message.

* **func**:
  * Updated import statements: changed `clusterApply`, `clusterCall`, `parRapply`, and `parCapply` from BiocGenerics to parallel package.

# monocle 2.6.1

* **bugs**:
  * `clusterCells()`: Rolled back to previous densityPeak clustering algorithm as default. The knn-based density peak clustering was not general for all datasets.
  * `importCDS()`, `exportCDS()`: Fixed various bugs.

* **func**:
  * `clusterCells()`: Added new Louvain clustering algorithm for dealing with large datasets (> 50k cells).

# monocle 2.6.0

* **release**:
  * Official release of the features in 2.5.0 through BioC.

# monocle 2.5.0

* **func**:
  * Added `plot_complex_cell_trajectory()`: New utility for visualizing complex developmental trajectory.
  * Added `plot_multiple_branches_pseudotime()`: Multi-way kinetic curve visualization.
  * Added `plot_multiple_branches_heatmap()`: Multi-way heatmap visualization.
  * `clusterCells()`: Now supports clustering for 100k+ cells using kNN based density peak clustering algorithm. Requires densityClust package version 0.3.
  * Added `importCDS()`: Function for importing scater or Seurat objects into monocle.
  * Added `exportCDS()`: Function for exporting monocle CellDataSet to scater or Seurat formats.

# monocle 2.4.0

* **func**:
  * Changed default `expressionFamily` from Tobit to `negbinomial.size()`. Users with TPM or FPKM data are urged to convert to relative transcript counts with `relative2abs()`.
  * `clusterCells()`: Revamped functionality based on t-SNE and densityPeak.
  * Added new procedure for selected ordering genes called "dpFeature". See vignette for details.

# monocle 2.1.1

* **bugs**:
  * `reduceDimension()`: Fixed problem where function would return different results on repeated runs given the same inputs. The issue was in DDRTree's kmeans and irlba implementations. Now uses deterministically initialized eigenvectors and deterministically selected rows.
  * `classifyCells()`: Fixed problem related to joining factors and levels that generated annoying warnings.
  * `differentialGeneTest()`: Fixed check for valid sizeFactors prior to testing. Without this check, the function would report FAIL on all genes because of factors not having enough levels.

# monocle 2.1.0

* **func**:
  * `relative2abs()`: Re-designed Census algorithm for converting relative expression values (e.g. TPMs) into absolute transcript counts. The new version is much more accurate. The interface has changed, and output values will be quite different from the old version. Now reports estimates of mRNA counts in the lysate instead of cDNA counts.
  * Added `plot_pseudotime_heatmap()`: New heatmap function that replaces the old `plot_genes_heatmap()`.
  * Added extensive new documentation.

* **bugs**:
  * `orderCells()`: Fixed issue where ordering cells with DDRTree would compress cells at the tips of trajectories.
  * Fixed pseudo counts application to work correctly regardless of underlying distribution used to model expression.
  * Fixed variance stabilization to be applied correctly when a CellDataSet object's expressionFamily is `negbinomial.size()`.
  * `calculateMarkerSpecificity()`: Changed gene_id field from factor to character, fixing indexing errors and nonsensical semi-supervised clustering and ordering results.
  * `BEAM()`, `plot_branched_heatmap()`: Fixed sparse matrix issue by using `cBind` instead of `cbind`.

# monocle 1.99.0

* **release**:
  * First public release of the Monocle 2 series. For a summary of new features and changes, please see: http://cole-trapnell-lab.github.io/monocle-release/features/

# monocle 1.1.5

* **data**:
  * Changed data layout in HSMMSingleCell package due to Bioconductor build issues related to VGAM updates. This caused some changes in the vignette.

# monocle 1.1.1

* **bugs**:
  * `responseMatrix()`: Fixed bug that occurs when you don't have any genes that fail VGAM fitting.
