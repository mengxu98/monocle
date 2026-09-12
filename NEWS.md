# monocle

# monocle 2.9.3

* **func**:
  * `orderCells()` and `reduceDimension()` DDRTree hot paths now run natively in
    C++ (ported from scop): cell-to-MST projection plus the projected-cell
    minimum spanning tree (`Prim`, replacing the complete-graph `stats::dist()`
    + `igraph::mst()` construction), the tree ordering traversal behind
    `extract_ddrtree_ordering()`, and state-based root selection in
    `select_root_cell()`. Outputs match the R implementation (pseudotime,
    states, and root cell are identical); when several projected cells share
    identical coordinates the spanning tree may pick a different but equally
    minimal arrangement among those tied cells, which only perturbs the
    pseudotime of the tied cells themselves. The fast paths are used by
    default and fall back to the original R code on error; set
    `options(monocle.fast_ordering = FALSE)` to disable them.
  * Ordering a 3,000-cell dataset drops from ~9.3 s to ~1.2 s; the R path
    scales quadratically with cell number while the C++ path stays flat.
  * DDRTree ordering no longer materializes the dense cell-by-cell distance
    matrix for large datasets. `project2MST()` reuses the projected-cell MST
    edge weights returned by the C++ primitive and stores `cellPairwiseDistances`
    as a sparse `Matrix` once a dataset exceeds
    `getOption("monocle.max_dense_pairwise_cells", 10000)` cells. A dense
    `n_cells x n_cells` numeric matrix costs 8 * n_cells^2 bytes (~0.8 GB at 10k
    cells, ~7 GB at 30k cells), which is what made large datasets impossible to
    order; below the limit the stored matrix is unchanged.
    `extract_ddrtree_ordering()` now indexes the distances instead of coercing
    them to a dense matrix, so it works with either representation.

* **bugs**:
  * `is_negbinomial_cds()` returns a single logical again. `negbinomial()` carries
    `vfamily = c("negbinomial", "VGAMcategorical")`, so the old `%in%` result had
    length two and callers such as `differentialGeneTest()`, `genSmoothCurves()`
    and `genSmoothCurveResiduals()` failed with
    `'length = 2' in coercion to 'logical(1)'`.
  * Genes with zero counts in every cell are reported as `status = "OK"` with
    `pval = 1` (and zero fitted curves/residuals) instead of `status = "FAIL"`
    with `NA`s, matching the VGAM path. Dropping them from the BH adjustment had
    rescaled every q-value in the dataset (~30% smaller when 30% of the genes
    were all-zero), so the fast paths are drop-in replacements again.
  * `estimateDispersions()`, `reduceDimension()`, and related row-variance
    calculations no longer convert sparse expression matrices to dense
    `TsparseMatrix` objects via `(x - rowMeans(x))^2`. That path overflows R's
    32-bit index limit when `nrow * ncol` exceeds `2^31 - 1`. Variance is now
    `E[X^2] - mean(X)^2`, which stays sparse.

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
