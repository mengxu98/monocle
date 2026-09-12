fast_ns_fun <- function(name) {
  getFromNamespace(name, "monocle")
}

force_vgam_path <- function(code) {
  ns <- asNamespace("monocle")
  orig <- getFromNamespace("is_negbinomial_cds", "monocle")
  unlockBinding("is_negbinomial_cds", ns)
  assign("is_negbinomial_cds", function(cds) FALSE, envir = ns)
  lockBinding("is_negbinomial_cds", ns)
  on.exit({
    unlockBinding("is_negbinomial_cds", ns)
    assign("is_negbinomial_cds", orig, envir = ns)
    lockBinding("is_negbinomial_cds", ns)
  }, add = TRUE)
  force(code)
}

make_synthetic_nb_cds <- function(n_genes = 40, n_cells = 80, seed = 42) {
  set.seed(seed)
  sf <- stats::runif(n_cells, 0.7, 1.3)
  pt <- seq(0, 10, length.out = n_cells)
  branch <- factor(rep(c("A", "B"), length.out = n_cells))
  mu <- outer(stats::runif(n_genes, 0.5, 8), exp(0.15 * (pt - mean(pt))))
  half <- max(1L, floor(n_genes / 2))
  mu[seq_len(half), branch == "B"] <- mu[seq_len(half), branch == "B"] * 2.5
  theta <- 5
  exprs_mat <- matrix(0, n_genes, n_cells)
  for (g in seq_len(n_genes)) {
    exprs_mat[g, ] <- stats::rnbinom(n_cells, size = theta, mu = mu[g, ] * sf)
  }
  rownames(exprs_mat) <- paste0("g", seq_len(n_genes))
  colnames(exprs_mat) <- paste0("c", seq_len(n_cells))
  pd <- data.frame(
    Pseudotime = pt,
    Branch = branch,
    Size_Factor = sf,
    row.names = colnames(exprs_mat)
  )
  fd <- data.frame(
    gene_short_name = rownames(exprs_mat),
    row.names = rownames(exprs_mat)
  )
  cds <- monocle::newCellDataSet(
    exprs_mat,
    phenoData = Biobase::AnnotatedDataFrame(pd),
    featureData = Biobase::AnnotatedDataFrame(fd),
    expressionFamily = VGAM::negbinomial.size()
  )
  cds <- BiocGenerics::estimateSizeFactors(cds)
  cds <- BiocGenerics::estimateDispersions(cds)
  cds
}

jaccard_sets <- function(a, b) {
  u <- sum(a | b)
  if (u == 0) {
    return(1)
  }
  sum(a & b) / u
}