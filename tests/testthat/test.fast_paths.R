context("C++ fast paths")

test_that("find_closest_point_cpp uses all dimensions", {
  find_closest_point_cpp <- fast_ns_fun("find_closest_point_cpp")
  set.seed(7)
  Z <- cbind(c(0, 0, 0), c(0, 0, 5))
  Y <- cbind(c(0, 0, 0), c(0, 0, 10), c(8, 0, 5))
  idx <- as.integer(find_closest_point_cpp(Z, Y, 1L))
  dmat <- as.matrix(dist(t(cbind(Z, Y))))[1:2, 3:5]
  idx_r <- apply(dmat, 1, which.min)
  expect_equal(idx, as.integer(idx_r))
})

test_that("find_closest_point_cpp matches full-dimension R nearest neighbor", {
  find_closest_point_cpp <- fast_ns_fun("find_closest_point_cpp")
  set.seed(1)
  D <- 4
  N <- 20
  K <- 7
  Z <- matrix(rnorm(D * N), D)
  Y <- matrix(rnorm(D * K), D)
  idx <- as.integer(find_closest_point_cpp(Z, Y, 1L))
  idx_r <- integer(N)
  for (i in seq_len(N)) {
    idx_r[i] <- which.min(colSums((Y - Z[, i])^2))
  }
  expect_equal(idx, idx_r)
})

test_that("find_closest_point_cpp is stable across thread counts", {
  find_closest_point_cpp <- fast_ns_fun("find_closest_point_cpp")
  set.seed(3)
  Z <- matrix(rnorm(5 * 50), 5)
  Y <- matrix(rnorm(5 * 11), 5)
  expect_equal(
    as.integer(find_closest_point_cpp(Z, Y, 1L)),
    as.integer(find_closest_point_cpp(Z, Y, 4L))
  )
})

test_that("nb_formula_for_lm rewrites sm.ns to ns", {
  nb_formula_for_lm <- fast_ns_fun("nb_formula_for_lm")
  expect_equal(
    nb_formula_for_lm("~sm.ns(Pseudotime, df=3)*Branch"),
    "~splines::ns(Pseudotime, df=3)*Branch"
  )
})

test_that("nb_design_for_newdata keeps spline knots from training data", {
  nb_design_for_newdata <- fast_ns_fun("nb_design_for_newdata")
  train <- data.frame(
    Pseudotime = seq(0, 10, length.out = 21),
    Branch = factor(rep(c("A", "B"), length.out = 21))
  )
  newdata <- data.frame(
    Pseudotime = c(1, 5, 9),
    Branch = factor(c("A", "B", "A"), levels = levels(train$Branch))
  )
  des <- nb_design_for_newdata("~sm.ns(Pseudotime, df=3)*Branch", train, newdata)
  expect_equal(ncol(des$X_fit), ncol(des$X_new))
  expect_equal(colnames(des$X_fit), colnames(des$X_new))
  expect_equal(nrow(des$X_new), 3)
})

test_that("differentialGeneTest C++ matches VGAM ranks and calls on synthetic NB", {
  cds <- make_synthetic_nb_cds(n_genes = 40, n_cells = 80, seed = 42)

  res_cpp <- monocle::differentialGeneTest(
    cds,
    fullModelFormulaStr = "~Branch",
    reducedModelFormulaStr = "~1",
    cores = 1
  )
  res_vgam <- force_vgam_path(
    monocle::differentialGeneTest(
      cds,
      fullModelFormulaStr = "~Branch",
      reducedModelFormulaStr = "~1",
      cores = 1
    )
  )

  expect_true(all(res_cpp$status %in% c("OK", "FAIL")))
  expect_true(all(res_vgam$status %in% c("OK", "FAIL")))

  both_ok <- res_cpp$status == "OK" & res_vgam$status == "OK"
  expect_gt(mean(both_ok), 0.85)

  pf <- res_cpp$pval[both_ok]
  pv <- res_vgam$pval[both_ok]
  expect_gt(suppressWarnings(cor(pf, pv, method = "spearman")), 0.99)

  log_diff <- abs(log10(pmax(pf, 1e-300)) - log10(pmax(pv, 1e-300)))
  expect_lt(median(log_diff), 0.05)
  expect_lt(max(log_diff), 0.5)

  for (alpha in c(0.05, 0.01)) {
    expect_gt(jaccard_sets(pf < alpha, pv < alpha), 0.95)
  }
})

test_that("differentialGeneTest C++ p-values match across cores", {
  cds <- make_synthetic_nb_cds(n_genes = 30, n_cells = 60, seed = 7)
  r1 <- monocle::differentialGeneTest(
    cds,
    fullModelFormulaStr = "~Branch",
    reducedModelFormulaStr = "~1",
    cores = 1
  )
  r4 <- monocle::differentialGeneTest(
    cds,
    fullModelFormulaStr = "~Branch",
    reducedModelFormulaStr = "~1",
    cores = 4
  )
  expect_equal(r1$status, r4$status)
  expect_equal(r1$pval, r4$pval, tolerance = 1e-12)
})

test_that("genSmoothCurves C++ matches VGAM response curves on synthetic NB", {
  cds <- make_synthetic_nb_cds(n_genes = 20, n_cells = 60, seed = 11)
  new_data <- data.frame(
    Pseudotime = seq(0, 10, length.out = 15),
    Branch = factor("A", levels = levels(Biobase::pData(cds)$Branch))
  )

  M_cpp <- monocle::genSmoothCurves(
    cds,
    new_data = new_data,
    trend_formula = "~sm.ns(Pseudotime, df=3)",
    relative_expr = TRUE,
    response_type = "response",
    cores = 1
  )
  M_vgam <- force_vgam_path(
    monocle::genSmoothCurves(
      cds,
      new_data = new_data,
      trend_formula = "~sm.ns(Pseudotime, df=3)",
      relative_expr = TRUE,
      response_type = "response",
      cores = 1
    )
  )

  expect_equal(dim(M_cpp), dim(M_vgam))
  ok <- which(is.finite(rowMeans(M_cpp)) & is.finite(rowMeans(M_vgam)))
  expect_gt(length(ok), 0.8 * nrow(M_cpp))

  rel <- abs(M_cpp[ok, , drop = FALSE] - M_vgam[ok, , drop = FALSE]) /
    pmax(abs(M_vgam[ok, , drop = FALSE]), 1e-6)
  expect_lt(median(rel), 0.05)
  expect_lt(as.numeric(stats::quantile(rel, 0.95)), 0.15)

  expect_gt(
    suppressWarnings(cor(
      rowMeans(M_cpp[ok, , drop = FALSE]),
      rowMeans(M_vgam[ok, , drop = FALSE]),
      method = "spearman"
    )),
    0.99
  )
})

test_that("zero-count genes are reported as OK with p = 1, as the VGAM path does", {
  cds <- make_synthetic_nb_cds(n_genes = 40, n_cells = 60, seed = 5)
  counts <- Biobase::exprs(cds)
  zero_genes <- rownames(cds)[1:2]
  counts[zero_genes, ] <- 0
  Biobase::exprs(cds) <- counts

  res <- monocle::differentialGeneTest(
    cds,
    fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)",
    reducedModelFormulaStr = "~1",
    cores = 1
  )

  expect_identical(as.character(res[zero_genes, "status"]), c("OK", "OK"))
  expect_equal(as.numeric(res[zero_genes, "pval"]), c(1, 1))
  # q-values are adjusted over every gene, so the zero-count genes must not be
  # dropped from the BH denominator.
  expect_equal(res$qval, stats::p.adjust(res$pval, method = "BH"))
})

test_that("zero-count genes yield zero curves and residuals instead of NA", {
  cds <- make_synthetic_nb_cds(n_genes = 40, n_cells = 60, seed = 6)
  counts <- Biobase::exprs(cds)
  zero_genes <- rownames(cds)[1:2]
  counts[zero_genes, ] <- 0
  Biobase::exprs(cds) <- counts

  new_data <- data.frame(Pseudotime = seq(0, 10, length.out = 5))
  row.names(new_data) <- paste0("pt", seq_len(nrow(new_data)))

  curves <- monocle::genSmoothCurves(
    cds,
    new_data = new_data,
    trend_formula = "~sm.ns(Pseudotime, df=3)",
    cores = 1
  )
  expect_false(any(is.na(curves[zero_genes, , drop = FALSE])))
  expect_equal(unname(curves[zero_genes, , drop = FALSE]), matrix(0, 2, nrow(new_data)))

  residuals <- fast_ns_fun("genSmoothCurveResiduals")(
    cds,
    trend_formula = "~sm.ns(Pseudotime, df=3)",
    residual_type = "response",
    cores = 1
  )
  residuals <- as.matrix(residuals)
  expect_false(any(is.na(residuals[zero_genes, , drop = FALSE])))
  expect_equal(unname(residuals[zero_genes, , drop = FALSE]), matrix(0, 2, ncol(cds)))
})

test_that("negbinomial() CellDataSets are detected without a length-2 condition", {
  is_negbinomial_cds <- fast_ns_fun("is_negbinomial_cds")
  counts <- matrix(
    stats::rpois(6 * 8, 4),
    nrow = 6,
    dimnames = list(paste0("g", 1:6), paste0("c", 1:8))
  )
  pd <- data.frame(Pseudotime = seq_len(8), row.names = colnames(counts))
  fd <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  cds <- monocle::newCellDataSet(
    counts,
    phenoData = Biobase::AnnotatedDataFrame(pd),
    featureData = Biobase::AnnotatedDataFrame(fd),
    expressionFamily = VGAM::negbinomial()
  )

  # vfamily for negbinomial() is c("negbinomial", "VGAMcategorical"), so the
  # helper has to collapse it: a length-2 result breaks if () and && callers.
  detected <- is_negbinomial_cds(cds)
  expect_length(detected, 1L)
  expect_identical(detected, TRUE)
  expect_true(is_negbinomial_cds(cds) && TRUE)

  size_cds <- make_synthetic_nb_cds(n_genes = 6, n_cells = 20, seed = 2)
  expect_identical(is_negbinomial_cds(size_cds), TRUE)
})

test_that("negbinomial() CellDataSets run through dispersions and DE", {
  set.seed(8)
  n_genes <- 40
  n_cells <- 60
  size_factor <- stats::runif(n_cells, 0.7, 1.3)
  pt <- seq(0, 10, length.out = n_cells)
  mu <- outer(stats::runif(n_genes, 0.5, 8), exp(0.15 * (pt - mean(pt))))
  counts <- matrix(0, n_genes, n_cells)
  for (g in seq_len(n_genes)) {
    counts[g, ] <- stats::rnbinom(n_cells, size = 5, mu = mu[g, ] * size_factor)
  }
  rownames(counts) <- paste0("g", seq_len(n_genes))
  colnames(counts) <- paste0("c", seq_len(n_cells))
  pd <- data.frame(Pseudotime = pt, row.names = colnames(counts))
  fd <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  cds <- monocle::newCellDataSet(
    counts,
    phenoData = Biobase::AnnotatedDataFrame(pd),
    featureData = Biobase::AnnotatedDataFrame(fd),
    expressionFamily = VGAM::negbinomial()
  )

  # estimateDispersions() and the size factor guard used to compare vfamily with
  # ==, which is length two for this family and errored out. Outlier removal is
  # off because the refit on the trimmed table is fragile for synthetic data.
  cds <- BiocGenerics::estimateSizeFactors(cds)
  cds <- BiocGenerics::estimateDispersions(cds, remove_outliers = FALSE)

  res <- monocle::differentialGeneTest(
    cds,
    fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)",
    reducedModelFormulaStr = "~1",
    cores = 1
  )
  expect_identical(unique(as.character(res$status)), "OK")
  expect_identical(unique(as.character(res$family)), "negbinomial")

  new_data <- data.frame(Pseudotime = seq(0, 1, length.out = 5))
  row.names(new_data) <- paste0("pt", seq_len(nrow(new_data)))
  curves <- monocle::genSmoothCurves(
    cds,
    new_data = new_data,
    trend_formula = "~sm.ns(Pseudotime, df=3)",
    cores = 1
  )
  expect_true(all(is.finite(curves)))
})