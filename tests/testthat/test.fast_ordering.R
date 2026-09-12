make_branching_counts <- function(n_genes = 60, n_cells = 80, seed = 42) {
  set.seed(seed)
  n_per <- ceiling(n_cells / 3)
  branch <- rep(c("root", "a", "b"), length.out = n_cells)
  t <- rep(seq(0, 1, length.out = n_per), length.out = n_cells)
  mu <- matrix(1, n_genes, n_cells)
  mu[1:20, branch == "root"] <- 5 + t[branch == "root"] * 2
  mu[21:40, branch == "a"] <- 3 + t[branch == "a"] * 8
  mu[41:60, branch == "b"] <- 3 + t[branch == "b"] * 8
  counts <- matrix(
    stats::rpois(n_genes * n_cells, lambda = as.vector(mu)),
    nrow = n_genes,
    dimnames = list(paste0("g", seq_len(n_genes)), paste0("c", seq_len(n_cells)))
  )
  counts
}

make_branching_cds <- function(n_genes = 60, n_cells = 80, seed = 42) {
  counts <- make_branching_counts(n_genes, n_cells, seed)
  pd <- data.frame(x = seq_len(n_cells), row.names = colnames(counts))
  fd <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  cds <- monocle::newCellDataSet(
    Matrix::Matrix(counts, sparse = TRUE),
    phenoData = Biobase::AnnotatedDataFrame(pd),
    featureData = Biobase::AnnotatedDataFrame(fd),
    expressionFamily = VGAM::uninormal()
  )
  monocle::reduceDimension(cds, max_components = 2, norm_method = "none", pseudo_expr = 0)
}

compare_orderings <- function(cds_fast, cds_r, tolerance = 1e-8) {
  expect_equal(
    as.numeric(pData(cds_fast)$Pseudotime),
    as.numeric(pData(cds_r)$Pseudotime),
    tolerance = tolerance
  )
  expect_identical(
    as.character(pData(cds_fast)$State),
    as.character(pData(cds_r)$State)
  )
  expect_identical(
    cds_fast@auxOrderingData[["DDRTree"]]$root_cell,
    cds_r@auxOrderingData[["DDRTree"]]$root_cell
  )
}

test_that("ddrtree_order_from_edges_cpp matches igraph::dfs-based ordering", {
  set.seed(11)
  for (trial in seq_len(20L)) {
    n <- sample(2:60, 1)
    edges <- matrix(0L, n - 1, 2)
    for (v in 2:n) edges[v - 1, ] <- c(sample(v - 1, 1), v)
    edges <- edges[sample(n - 1), , drop = FALSE]
    weights <- runif(n - 1, 0.01, 5)
    vnames <- paste0("v", sample(n))
    root <- sample(n, 1)

    g <- igraph::graph_from_edgelist(edges, directed = FALSE)
    igraph::V(g)$name <- vnames
    dp <- matrix(0, n, n, dimnames = list(vnames, vnames))
    dp[cbind(vnames[edges[, 1]], vnames[edges[, 2]])] <- weights
    dp[cbind(vnames[edges[, 2]], vnames[edges[, 1]])] <- weights

    traversal <- igraph::dfs(
      g, root = vnames[root], mode = "all",
      unreachable = FALSE, parent = TRUE
    )
    mst_parent <- as.numeric(traversal$parent)
    states <- rep(1, n)
    pseudotimes <- rep(0, n)
    parents <- rep(NA_character_, n)
    names(states) <- vnames
    names(pseudotimes) <- vnames
    names(parents) <- vnames
    curr_state <- 1
    for (i in seq_along(traversal$order)) {
      curr_name <- igraph::V(g)[traversal$order[i]]$name
      if (!is.na(mst_parent[traversal$order[i]])) {
        par_name <- igraph::V(g)[mst_parent[traversal$order[i]]]$name
        if (igraph::degree(g, v = par_name) > 2) curr_state <- curr_state + 1
        pseudotimes[curr_name] <- pseudotimes[par_name] + dp[curr_name, par_name]
        parents[curr_name] <- par_name
      }
      states[curr_name] <- curr_state
    }

    out <- fast_ns_fun("ddrtree_order_from_edges_cpp")(n, edges, weights, root)
    expect_identical(as.character(factor(states)), as.character(factor(out$state)))
    expect_equal(unname(pseudotimes), as.numeric(out$pseudotime), tolerance = 0)
    cpp_parents <- rep(NA_character_, n)
    hit <- out$parent > 0L
    cpp_parents[hit] <- vnames[out$parent[hit]]
    expect_identical(unname(parents), cpp_parents)
  }
})

test_that("orderCells fast ordering matches the R path on small data", {
  cds_fast <- make_branching_cds(n_cells = 80)
  cds_r <- make_branching_cds(n_cells = 80)

  options(monocle.fast_ordering = FALSE)
  cds_r <- orderCells(cds_r)
  options(monocle.fast_ordering = TRUE)
  cds_fast <- orderCells(cds_fast)
  compare_orderings(cds_fast, cds_r)

  states <- sort(unique(as.character(pData(cds_fast)$State)))
  states <- states[!is.na(states)]
  skip_if(length(states) < 2L)

  options(monocle.fast_ordering = FALSE)
  cds_r <- orderCells(cds_r, root_state = states[length(states)])
  options(monocle.fast_ordering = TRUE)
  cds_fast <- orderCells(cds_fast, root_state = states[length(states)])
  compare_orderings(cds_fast, cds_r)
})

test_that("orderCells fast ordering matches the R path with DDRTree centers", {
  cds_fast <- make_branching_cds(n_cells = 150)
  cds_r <- make_branching_cds(n_cells = 150)

  options(monocle.fast_ordering = FALSE)
  cds_r <- orderCells(cds_r)
  options(monocle.fast_ordering = TRUE)
  cds_fast <- orderCells(cds_fast)
  compare_orderings(cds_fast, cds_r)
})

test_that("fast projection leaves a valid spanning projected tree", {
  cds <- make_branching_cds(n_cells = 80)
  cds <- orderCells(cds)

  tree <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree
  expect_true(igraph::ecount(tree) == ncol(cds) - 1L)
  expect_true(igraph::vcount(tree) == ncol(cds))
  expect_identical(sort(igraph::V(tree)$name), sort(colnames(cds)))
  expect_true(igraph::is_connected(tree))
  expect_true(all(is.finite(igraph::E(tree)$weight)))
  expect_true(all(igraph::E(tree)$weight > 0))
  expect_identical(
    row.names(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex),
    colnames(cds)
  )
  expect_true(all(is.finite(pData(cds)$Pseudotime)))
})

test_that("large cell counts keep only sparse MST distances and order identically", {
  project2MST <- fast_ns_fun("project2MST")
  project_point_to_line_segment <- fast_ns_fun("project_point_to_line_segment")
  cell_mst_adjacency <- fast_ns_fun("cell_mst_adjacency")

  cds_dense <- make_branching_cds(n_cells = 80)
  options(monocle.max_dense_pairwise_cells = 0L)
  cds_sparse <- make_branching_cds(n_cells = 80)
  on.exit(options(monocle.max_dense_pairwise_cells = NULL), add = TRUE)

  cds_sparse <- orderCells(cds_sparse)
  options(monocle.max_dense_pairwise_cells = NULL)
  cds_dense <- orderCells(cds_dense)
  compare_orderings(cds_sparse, cds_dense)

  # force the sparse representation and check it still describes the tree
  options(monocle.max_dense_pairwise_cells = 0L)
  cds_proj <- project2MST(cds_dense, project_point_to_line_segment)
  options(monocle.max_dense_pairwise_cells = NULL)

  dp <- cellPairwiseDistances(cds_proj)
  n_cells <- length(colnames(cds_dense))
  expect_s4_class(dp, "Matrix")
  expect_identical(unname(dim(dp)), c(n_cells, n_cells))
  expect_identical(dimnames(dp), list(colnames(cds_dense), colnames(cds_dense)))

  tree <- cds_proj@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree
  edges <- igraph::as_edgelist(tree, names = TRUE)
  idx <- cbind(match(edges[, 1], colnames(cds_dense)), match(edges[, 2], colnames(cds_dense)))
  expect_equal(as.numeric(dp[idx]), igraph::E(tree)$weight)
  expect_true(all(is.finite(dp[idx])))
  expect_true(all(dp[idx] > 0))

  # and the same tree, expressed densely, carries the same weights
  adj <- cell_mst_adjacency(
    matrix(match(edges, colnames(cds_dense)), ncol = 2L),
    igraph::E(tree)$weight,
    n_cells,
    colnames(cds_dense)
  )
  expect_equal(as.matrix(adj), as.matrix(dp))
})
