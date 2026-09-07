run_pca <- function(data, ...) {
  res <- stats::prcomp(data, center = F, scale = F)
  res$x
}

scale_pseudotime <- function(cds, verbose = FALSE) {
  Parent <- NA
  pd <- pData(cds)
  pd$Cell_name <- row.names(pd)
  range_df <- plyr::ddply(pd, .(State), function(x) {
    min_max <- range(x$Pseudotime)
    min_cell <- subset(x, Pseudotime %in% min_max[1])
    max_cell <- subset(x, Pseudotime %in% min_max[2])
    min_cell$fate_type <- "Start"
    max_cell$fate_type <- "End"
    rbind(min_cell, max_cell)
  })

  adj_list <- data.frame(Source = subset(range_df, length(Parent) > 0 & !is.na(Parent))[, "Parent"], Target = subset(range_df, length(Parent) > 0 & !is.na(Parent))[, "Cell_name"])
  adj_list$Source <- pd[as.character(adj_list$Source), "State"]
  adj_list$Target <- pd[as.character(adj_list$Target), "State"]

  uniq_cell_list <- unique(c(as.character(adj_list$Source), as.character(adj_list$Target)))
  adj_mat <- matrix(rep(0, length(uniq_cell_list)^2), nrow = length(uniq_cell_list), ncol = length(uniq_cell_list), dimnames = list(uniq_cell_list, uniq_cell_list))
  adj_mat[thisutils::as_matrix(adj_list)] <- 1
  net <- igraph::graph_from_adjacency_matrix(thisutils::as_matrix(adj_mat), mode = "directed", weighted = NULL, diag = FALSE)

  net_leaves <- which(igraph::degree(net, v = V(net), mode = "out") == 0, useNames = T)

  pd$scale_pseudotime <- NA

  for (i in net_leaves) {
    path_vertex <- as.vector(
      igraph::all_shortest_paths(
        net,
        from = 1, to = i, mode = "out"
      )$res[[1]]
    )
    pd_subset <- subset(pd, State %in% path_vertex & is.na(scale_pseudotime))

    min_cell_name <- row.names(subset(pd_subset, Pseudotime == min(Pseudotime)))

    if (!is.na(pd[min_cell_name, "Parent"])) {
      parent_min_cell <- as.character(pd[min_cell_name, "Parent"])
      subset_min_pseudo <- pd[parent_min_cell, "Pseudotime"]
      scale_pseudotime_ini <- pd[parent_min_cell, "scale_pseudotime"]
      scaling_factor <- (100 - pd[parent_min_cell, "scale_pseudotime"]) / c(max(pd_subset$Pseudotime) - subset_min_pseudo)
    } else {
      subset_min_pseudo <- min(pd[, "Pseudotime"])
      scale_pseudotime_ini <- 0
      scaling_factor <- 100 / c(max(pd_subset$Pseudotime) - min(pd_subset$Pseudotime))
    }

    pseudotime_scaled <- (pd_subset$Pseudotime - subset_min_pseudo) * scaling_factor + scale_pseudotime_ini

    log_message("{i}\t{range(pseudotime_scaled)[1]}\t{range(pseudotime_scaled)[2]}", verbose = verbose)

    pd[row.names(pd_subset), "ori_pseudotime"] <- pd[row.names(pd_subset), "Pseudotime"]
    pd[row.names(pd_subset), "Pseudotime"] <- pseudotime_scaled
  }
  scale_pseudotime <- (pd_subset$Pseudotime - subset_min_pseudo) * scaling_factor + scale_pseudotime_ini
  log_message("{i}\t{range(scale_pseudotime)[1]}\t{range(scale_pseudotime)[2]}")
  pd[row.names(pd_subset), "scale_pseudotime"] <- scale_pseudotime

  pData(cds) <- pd

  return(cds)
}

# get_next_node_id <- function() {
#   next_node <<- next_node + 1
#   return(next_node)
# }
make_get_next_node_id <- function(start = 0L) {
  next_node <- start
  function() {
    next_node <<- next_node + 1L
    next_node
  }
}
get_next_node_id <- make_get_next_node_id()

#' Recursively builds and returns a PQ tree for the MST
#' @param mst The minimum spanning tree, as an igraph object.
#' @param use_weights Whether to use edge weights when finding the diameter path of the tree.
#' @param root_node The name of the root node to use for starting the path finding.
pq_helper <- function(mst, use_weights = TRUE, root_node = NULL) {
  new_subtree <- igraph::make_empty_graph()

  root_node_id <- paste("Q_", get_next_node_id(), sep = "")

  new_subtree <- new_subtree + igraph::vertex(root_node_id, type = "Q", color = "black")

  if (is.null(root_node) == FALSE) {
    sp <- igraph::all_shortest_paths(mst, from = V(mst)[root_node])
    sp_lengths <- sapply(sp$res, length)
    target_node_idx <- which(sp_lengths == max(sp_lengths))[1]
    diam <- V(mst)[unlist(sp$res[target_node_idx])]
  } else {
    if (use_weights) {
      diam <- V(mst)[igraph::get_diameter(mst)]
    } else {
      diam <- V(mst)[igraph::get_diameter(mst, weights = NA)]
    }
  }

  V(new_subtree)[root_node_id]$diam_path_len <- length(diam)

  diam_decisiveness <- igraph::degree(mst, v = diam) > 2
  ind_nodes <- diam_decisiveness[diam_decisiveness == TRUE]

  first_diam_path_node_idx <- head(as.vector(diam), n = 1)
  last_diam_path_node_idx <- tail(as.vector(diam), n = 1)
  if (sum(ind_nodes) == 0 ||
    (igraph::degree(mst, first_diam_path_node_idx) == 1 &&
      igraph::degree(mst, last_diam_path_node_idx) == 1)) {
    ind_backbone <- diam
  } else {
    last_bb_point <- names(tail(ind_nodes, n = 1))[[1]]
    first_bb_point <- names(head(ind_nodes, n = 1))[[1]]
    diam_path_names <- V(mst)[as.vector(diam)]$name
    last_bb_point_idx <- which(diam_path_names == last_bb_point)[1]
    first_bb_point_idx <- which(diam_path_names == first_bb_point)[1]
    ind_backbone_idxs <- as.vector(diam)[first_bb_point_idx:last_bb_point_idx]
    ind_backbone <- V(mst)[ind_backbone_idxs]
  }

  mst_no_backbone <- mst - ind_backbone

  for (backbone_n in ind_backbone) {
    if (igraph::degree(mst, v = backbone_n) > 2) {
      new_p_id <- paste("P_", get_next_node_id(), sep = "")
      new_subtree <- new_subtree + igraph::vertex(new_p_id, type = "P", color = "grey")
      new_subtree <- new_subtree + igraph::vertex(V(mst)[backbone_n]$name, type = "leaf", color = "white")
      new_subtree <- new_subtree + igraph::edge(new_p_id, V(mst)[backbone_n]$name)
      new_subtree <- new_subtree + igraph::edge(root_node_id, new_p_id)

      nb <- igraph::make_ego_graph(mst, 1, nodes = backbone_n)[[1]]

      for (n_i in V(nb))
      {
        n <- V(nb)[n_i]$name
        if (n %in% V(mst_no_backbone)$name) {
          sc <- igraph::subcomponent(mst_no_backbone, n)

          sg <- igraph::induced_subgraph(mst_no_backbone, sc, impl = "copy_and_delete")

          if (igraph::ecount(sg) > 0) {
            sub_pq <- pq_helper(sg, use_weights)

            for (v in V(sub_pq$subtree))
            {
              new_subtree <- new_subtree + igraph::vertex(V(sub_pq$subtree)[v]$name, type = V(sub_pq$subtree)[v]$type, color = V(sub_pq$subtree)[v]$color, diam_path_len = V(sub_pq$subtree)[v]$diam_path_len)
            }

            edge_list <- igraph::as_edgelist(sub_pq$subtree)
            for (i in 1:nrow(edge_list))
            {
              new_subtree <- new_subtree + igraph::edge(V(sub_pq$subtree)[edge_list[i, 1]]$name, V(sub_pq$subtree)[edge_list[i, 2]]$name)
            }

            new_subtree <- new_subtree + igraph::edge(new_p_id, V(sub_pq$subtree)[sub_pq$root]$name)
          } else {
            new_subtree <- new_subtree + igraph::vertex(n, type = "leaf", color = "white")
            new_subtree <- new_subtree + igraph::edge(new_p_id, n)
          }
        }
      }
    } else {
      new_subtree <- new_subtree + igraph::vertex(V(mst)[backbone_n]$name, type = "leaf", color = "white")
      new_subtree <- new_subtree + igraph::edge(root_node_id, V(mst)[backbone_n]$name)
    }
  }

  return(list(root = root_node_id, subtree = new_subtree))
}

make_canonical <- function(pq_tree) {
  type <- NA
  canonical_pq <- pq_tree

  V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode = "out") == 2]$color <- "black"
  V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode = "out") == 2]$type <- "Q"

  single_child_p <- V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode = "out") == 1]
  V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode = "out") == 1]$color <- "blue"

  for (p_node in single_child_p)
  {
    child_of_p_node <- igraph::neighbors(canonical_pq, p_node, mode = "out")
    parent_of_p_node <- igraph::neighbors(canonical_pq, p_node, mode = "in")

    for (child_of_p in child_of_p_node)
    {
      canonical_pq[parent_of_p_node, child_of_p] <- TRUE
    }
  }

  canonical_pq <- igraph::delete_vertices(
    canonical_pq, V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, mode = "out") == 1]
  )
  return(canonical_pq)
}

count_leaf_descendents <- function(pq_tree, curr_node, children_counts) {
  if (V(pq_tree)[curr_node]$type == "leaf") {
    children_counts[curr_node] <- 0
    return(children_counts)
  } else {
    children_count <- 0
    children <- igraph::neighbors(pq_tree, curr_node, mode = "out")
    for (i in seq_along(children))
    {
      child <- names(children)[i]
      children_counts <- count_leaf_descendents(pq_tree, child, children_counts)
      if (V(pq_tree)[child]$type == "leaf") {
        children_count <- children_count + 1
      } else {
        children_count <- children_count + children_counts[child]
      }
    }
    children_counts[curr_node] <- children_count
    return(children_counts)
  }
}

#' Return an ordering for a P node in the PQ tree
#' @param q_level_list A list of Q nodes in the PQ tree
#' @param dist_matrix A symmetric matrix of pairwise distances between cells
order_p_node <- function(q_level_list, dist_matrix) {
  q_order_res <- combinat::permn(q_level_list, fun = order_q_node, dist_matrix)
  all_perms <- lapply(q_order_res, function(x) {
    x$ql
  })
  all_perms_weights <- unlist(lapply(q_order_res, function(x) {
    x$wt
  }))

  opt_perm_idx <- head((which(all_perms_weights == min(all_perms_weights))), 1)
  opt_perm <- all_perms[[opt_perm_idx]]

  stopifnot(length(opt_perm) == length(q_level_list))

  return(opt_perm)
}

order_q_node <- function(q_level_list, dist_matrix) {
  new_subtree <- igraph::make_empty_graph()

  if (length(q_level_list) == 1) {
    return(list(ql = q_level_list, wt = 0))
  }
  for (i in 1:length(q_level_list))
  {
    new_subtree <- new_subtree + igraph::vertex(paste(i, "F"), type = "forward")
    new_subtree <- new_subtree + igraph::vertex(paste(i, "R"), type = "reverse")
  }

  for (i in (1:(length(q_level_list) - 1)))
  {
    cost <- dist_matrix[q_level_list[[i]][length(q_level_list[[i]])], q_level_list[[i + 1]][1]]
    new_subtree <- new_subtree + igraph::edge(paste(i, "F"), paste(i + 1, "F"), weight = cost)

    cost <- dist_matrix[q_level_list[[i]][length(q_level_list[[i]])], q_level_list[[i + 1]][length(q_level_list[[i + 1]])]]
    new_subtree <- new_subtree + igraph::edge(paste(i, "F"), paste(i + 1, "R"), weight = cost)

    cost <- dist_matrix[q_level_list[[i]][1], q_level_list[[i + 1]][1]]
    new_subtree <- new_subtree + igraph::edge(paste(i, "R"), paste(i + 1, "F"), weight = cost)

    cost <- dist_matrix[q_level_list[[i]][1], q_level_list[[i + 1]][length(q_level_list[[i + 1]])]]
    new_subtree <- new_subtree + igraph::edge(paste(i, "R"), paste(i + 1, "R"), weight = cost)
  }

  first_fwd <- V(new_subtree)[paste(1, "F")]
  first_rev <- V(new_subtree)[paste(1, "R")]
  last_fwd <- V(new_subtree)[paste(length(q_level_list), "F")]
  last_rev <- V(new_subtree)[paste(length(q_level_list), "R")]

  FF_path <- unlist(
    igraph::shortest_paths(
      new_subtree,
      from = as.vector(first_fwd),
      to = as.vector(last_fwd),
      mode = "out",
      output = "vpath"
    )$vpath
  )
  FR_path <- unlist(
    igraph::shortest_paths(
      new_subtree,
      from = as.vector(first_fwd),
      to = as.vector(last_rev),
      mode = "out",
      output = "vpath"
    )$vpath
  )
  RF_path <- unlist(
    igraph::shortest_paths(
      new_subtree,
      from = as.vector(first_rev),
      to = as.vector(last_fwd),
      mode = "out",
      output = "vpath"
    )$vpath
  )
  RR_path <- unlist(
    igraph::shortest_paths(
      new_subtree,
      from = as.vector(first_rev),
      to = as.vector(last_rev),
      mode = "out",
      output = "vpath"
    )$vpath
  )

  FF_weight <- sum(E(new_subtree, path = FF_path)$weight)
  FR_weight <- sum(E(new_subtree, path = FR_path)$weight)
  RF_weight <- sum(E(new_subtree, path = RF_path)$weight)
  RR_weight <- sum(E(new_subtree, path = RR_path)$weight)

  paths <- list(FF_path, FR_path, RF_path, RR_path)
  path_weights <- c(FF_weight, FR_weight, RF_weight, RR_weight)
  opt_path_idx <- head((which(path_weights == min(path_weights))), 1)
  opt_path <- paths[[opt_path_idx]]

  stopifnot(length(opt_path) == length(q_level_list))

  directions <- V(new_subtree)[opt_path]$type
  q_levels <- list()
  for (i in 1:length(directions))
  {
    if (directions[[i]] == "forward") {
      q_levels[[length(q_levels) + 1]] <- q_level_list[[i]]
    } else {
      q_levels[[length(q_levels) + 1]] <- rev(q_level_list[[i]])
    }
  }

  return(list(ql = q_levels, wt = min(path_weights)))
}

measure_diameter_path <- function(pq_tree, curr_node, path_lengths) {
  if (V(pq_tree)[curr_node]$type != "Q") {
    path_lengths[curr_node] <- 0
    return(path_lengths)
  } else {
    children_count <- 0
    children <- igraph::neighbors(pq_tree, curr_node, mode = "out")
    for (i in seq_along(children))
    {
      child <- names(children)[i]
      children_counts <- count_leaf_descendents(pq_tree, child, children_counts)
      if (V(pq_tree)[child]$type == "leaf") {
        children_count <- children_count + 1
      } else {
        children_count <- children_count + children_counts[child]
      }
    }

    path_lengths[curr_node] <- children_count
    return(children_counts)
  }
}

assign_cell_lineage <- function(pq_tree, curr_node, assigned_state, node_states) {
  if (V(pq_tree)[curr_node]$type == "leaf") {
    node_states[V(pq_tree)[curr_node]$name] <- assigned_state
    return(node_states)
  } else {
    children <- igraph::neighbors(pq_tree, curr_node, mode = "out")
    for (i in seq_along(children))
    {
      child <- names(children)[i]
      node_states <- assign_cell_lineage(pq_tree, child, assigned_state, node_states)
    }
    return(node_states)
  }
}

extract_good_ordering <- function(pq_tree, curr_node, dist_matrix) {
  if (V(pq_tree)[curr_node]$type == "leaf") {
    return(V(pq_tree)[curr_node]$name)
  } else if (V(pq_tree)[curr_node]$type == "P") {
    p_level <- list()
    children <- igraph::neighbors(pq_tree, curr_node, mode = "out")
    for (i in seq_along(children))
    {
      child <- names(children)[i]
      p_level[[length(p_level) + 1]] <- extract_good_ordering(pq_tree, child, dist_matrix)
    }
    p_level <- order_p_node(p_level, dist_matrix)
    p_level <- unlist(p_level)
    return(p_level)
  } else if (V(pq_tree)[curr_node]$type == "Q") {
    q_level <- list()
    children <- igraph::neighbors(pq_tree, curr_node, mode = "out")
    for (i in seq_along(children))
    {
      child <- names(children)[i]
      q_level[[length(q_level) + 1]] <- extract_good_ordering(pq_tree, child, dist_matrix)
    }
    q_level <- order_q_node(q_level, dist_matrix)
    q_level <- q_level$ql
    q_level <- unlist(q_level)
    return(q_level)
  }
}

#' Extract a linear ordering of cells from a PQ tree
#'
#' @param orig_pq_tree The PQ object to use for ordering
#' @param curr_node The node in the PQ tree to use as the start of ordering
#' @param dist_matrix A symmetric matrix containing pairwise distances between cells
#' @param num_branches The number of outcomes allowed in the trajectory.
#' @param reverse_main_path Whether to reverse the direction of the trajectory
extract_good_branched_ordering <- function(orig_pq_tree, curr_node, dist_matrix, num_branches, reverse_main_path = FALSE) {
  requireNamespace("plyr")
  nei <- NULL
  type <- NA
  pseudo_time <- NA

  pq_tree <- orig_pq_tree

  branch_node_counts <- V(pq_tree)[type == "Q"]$diam_path_len
  names(branch_node_counts) <- V(pq_tree)[type == "Q"]$name
  if (length(names(branch_node_counts)) < num_branches) {
    stop("Number of branches attempted is larger than the branches constructed from pq_tree algorithm")
  }

  branch_node_counts <- sort(branch_node_counts, decreasing = TRUE)

  cell_states <- rep(NA, length(as.vector(V(pq_tree)[type == "leaf"])))
  names(cell_states) <- V(pq_tree)[type == "leaf"]$name

  cell_states <- assign_cell_lineage(pq_tree, curr_node, 1, cell_states)

  branch_point_roots <- list()

  branch_tree <- igraph::make_empty_graph()

  for (i in 1:num_branches)
  {
    branch_point_roots[[length(branch_point_roots) + 1]] <- names(branch_node_counts)[i]
    branch_id <- names(branch_node_counts)[i]
    branch_tree <- branch_tree + igraph::vertex(branch_id)
    parents <- igraph::neighbors(pq_tree, names(branch_node_counts)[i], mode = "in")
    if (length(parents) > 0 && V(pq_tree)[parents]$type == "P") {
      p_node_parent <- names(parents)[1]
      parent_branch_id <- names(igraph::neighbors(pq_tree, p_node_parent, mode = "in"))[1]
      branch_tree <- branch_tree + igraph::edge(parent_branch_id, branch_id)
    }
    pq_tree[parents, names(branch_node_counts)[i]] <- FALSE
  }

  branch_pseudotimes <- list()

  for (i in 1:length(branch_point_roots))
  {
    branch_ordering <- extract_good_ordering(pq_tree, branch_point_roots[[i]], dist_matrix)
    branch_ordering_time <- weight_of_ordering(branch_ordering, dist_matrix)
    names(branch_ordering_time) <- branch_ordering
    branch_pseudotimes[[length(branch_pseudotimes) + 1]] <- branch_ordering_time
    names(branch_pseudotimes)[length(branch_pseudotimes)] <- branch_point_roots[[i]]
  }

  cell_ordering_tree <- igraph::make_empty_graph()
  curr_branch <- "Q_1"

  extract_branched_ordering_helper <- function(branch_tree, curr_branch, cell_ordering_tree, branch_pseudotimes, dist_matrix, reverse_ordering = FALSE) {
    nei <- NULL

    curr_branch_pseudotimes <- branch_pseudotimes[[curr_branch]]
    curr_branch_root_cell <- NA
    for (i in 1:length(curr_branch_pseudotimes))
    {
      cell_ordering_tree <- cell_ordering_tree + igraph::vertex(names(curr_branch_pseudotimes)[i])
      if (i > 1) {
        if (reverse_ordering == FALSE) {
          cell_ordering_tree <- cell_ordering_tree + igraph::edge(names(curr_branch_pseudotimes)[i - 1], names(curr_branch_pseudotimes)[i])
        } else {
          cell_ordering_tree <- cell_ordering_tree + igraph::edge(names(curr_branch_pseudotimes)[i], names(curr_branch_pseudotimes)[i - 1])
        }
      }
    }

    if (reverse_ordering == FALSE) {
      curr_branch_root_cell <- names(curr_branch_pseudotimes)[1]
    } else {
      curr_branch_root_cell <- names(curr_branch_pseudotimes)[length(curr_branch_pseudotimes)]
    }

    children <- igraph::neighbors(branch_tree, curr_branch, mode = "out")
    for (j in seq_along(children))
    {
      child <- names(children)[j]
      child_cell_ordering_subtree <- igraph::make_empty_graph()

      child_head <- names(branch_pseudotimes[[child]])[1]
      child_tail <- names(branch_pseudotimes[[child]])[length(branch_pseudotimes[[child]])]

      curr_branch_cell_names <- names(branch_pseudotimes[[curr_branch]])
      head_dist_to_curr <- dist_matrix[child_head, curr_branch_cell_names]
      closest_to_head <- names(head_dist_to_curr)[which(head_dist_to_curr == min(head_dist_to_curr))]

      head_dist_to_anchored_branch <- NA
      branch_index_for_head <- NA

      head_dist_to_anchored_branch <- dist_matrix[closest_to_head, child_head]

      tail_dist_to_curr <- dist_matrix[child_tail, curr_branch_cell_names]
      closest_to_tail <- names(tail_dist_to_curr)[which(tail_dist_to_curr == min(tail_dist_to_curr))]

      tail_dist_to_anchored_branch <- NA
      branch_index_for_tail <- NA

      tail_dist_to_anchored_branch <- dist_matrix[closest_to_tail, child_tail]

      if (tail_dist_to_anchored_branch < head_dist_to_anchored_branch) {
        reverse_child <- TRUE
      } else {
        reverse_child <- FALSE
      }

      res <- extract_branched_ordering_helper(branch_tree, child, child_cell_ordering_subtree, branch_pseudotimes, dist_matrix, reverse_child)
      child_cell_ordering_subtree <- res$subtree
      child_subtree_root <- res$root

      for (v in V(child_cell_ordering_subtree))
      {
        cell_ordering_tree <- cell_ordering_tree + igraph::vertex(V(child_cell_ordering_subtree)[v]$name)
      }

      edge_list <- igraph::as_edgelist(child_cell_ordering_subtree)
      for (i in 1:nrow(edge_list))
      {
        cell_ordering_tree <- cell_ordering_tree + igraph::edge(V(cell_ordering_tree)[edge_list[i, 1]]$name, V(cell_ordering_tree)[edge_list[i, 2]]$name)
      }

      if (tail_dist_to_anchored_branch < head_dist_to_anchored_branch) {
        cell_ordering_tree <- cell_ordering_tree + igraph::edge(closest_to_tail, child_subtree_root)
      } else {
        cell_ordering_tree <- cell_ordering_tree + igraph::edge(closest_to_head, child_subtree_root)
      }
    }

    return(list(subtree = cell_ordering_tree, root = curr_branch_root_cell, last_cell_state = 1, last_cell_pseudotime = 0.0))
  }

  res <- extract_branched_ordering_helper(branch_tree, curr_branch, cell_ordering_tree, branch_pseudotimes, dist_matrix, reverse_main_path)
  cell_ordering_tree <- res$subtree

  curr_state <- 1

  assign_cell_state_helper <- function(ordering_tree_res, curr_cell) {
    nei <- NULL

    cell_tree <- ordering_tree_res$subtree
    V(cell_tree)[curr_cell]$cell_state <- curr_state

    children <- igraph::neighbors(cell_tree, curr_cell, mode = "out")
    ordering_tree_res$subtree <- cell_tree

    if (length(children) == 1) {
      ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, names(children)[1])
    } else {
      for (i in seq_along(children)) {
        curr_state <<- curr_state + 1
        ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, names(children)[i])
      }
    }
    return(ordering_tree_res)
  }

  res <- assign_cell_state_helper(res, res$root)

  assign_pseudotime_helper <- function(ordering_tree_res, dist_matrix, last_pseudotime, curr_cell) {
    nei <- NULL

    cell_tree <- ordering_tree_res$subtree
    curr_cell_pseudotime <- last_pseudotime
    V(cell_tree)[curr_cell]$pseudotime <- curr_cell_pseudotime
    parent_nodes <- igraph::neighbors(cell_tree, curr_cell, mode = "in")
    V(cell_tree)[curr_cell]$parent <- if (length(parent_nodes) > 0) names(parent_nodes)[1] else NA

    ordering_tree_res$subtree <- cell_tree
    children <- igraph::neighbors(cell_tree, curr_cell, mode = "out")

    for (i in seq_along(children)) {
      next_node <- names(children)[i]
      delta_pseudotime <- dist_matrix[curr_cell, next_node]
      ordering_tree_res <- assign_pseudotime_helper(ordering_tree_res, dist_matrix, last_pseudotime + delta_pseudotime, next_node)
    }

    return(ordering_tree_res)
  }

  res <- assign_pseudotime_helper(res, dist_matrix, 0.0, res$root)

  cell_names <- V(res$subtree)$name
  cell_states <- V(res$subtree)$cell_state
  cell_pseudotime <- V(res$subtree)$pseudotime
  cell_parents <- V(res$subtree)$parent
  ordering_df <- data.frame(
    sample_name = cell_names,
    cell_state = factor(cell_states),
    pseudo_time = cell_pseudotime,
    parent = cell_parents
  )

  ordering_df <- plyr::arrange(ordering_df, pseudo_time)
  return(list("ordering_df" = ordering_df, "cell_ordering_tree" = cell_ordering_tree))
}

reverse_ordering <- function(pseudo_time_ordering) {
  pt <- pseudo_time_ordering$pseudo_time
  names(pt) <- pseudo_time_ordering$sample_name
  rev_pt <- -((pt - max(pt)))
  rev_df <- pseudo_time_ordering
  rev_df$pseudo_time <- rev_pt
  return(rev_df)
}

weight_of_ordering <- function(ordering, dist_matrix) {
  time_delta <- c(0)
  curr_weight <- 0
  ep <- 0.01
  for (i in 2:length(ordering))
  {
    d <- dist_matrix[ordering[[i]], ordering[[i - 1]]]
    curr_weight <- curr_weight + d + ep
    time_delta <- c(time_delta, curr_weight)
  }

  return(time_delta)
}

#' Marks genes for clustering
#' @description The function marks genes that will be used for clustering in subsequent calls to clusterCells.
#' The list of selected genes can be altered at any time.
#'
#' @param cds the CellDataSet upon which to perform this operation
#' @param ordering_genes a vector of feature ids (from the CellDataSet's featureData) used for ordering cells
#' @return an updated CellDataSet object
#' @export
setOrderingFilter <- function(cds, ordering_genes) {
  fData(cds)$use_for_ordering <- row.names(fData(cds)) %in% ordering_genes
  cds
}

ica_helper <- function(X, n.comp, alg.typ = c("parallel", "deflation"), fun = c("logcosh", "exp"), alpha = 1,
                       row.norm = TRUE, maxit = 200, tol = 1e-4, verbose = FALSE, w.init = NULL, use_irlba = TRUE) {
  dd <- dim(X)
  d <- dd[dd != 1L]
  if (length(d) != 2L) {
    stop("data must be matrix-conformal")
  }
  X <- if (length(d) != length(dd)) {
    matrix(X, d[1L], d[2L])
  } else {
    thisutils::as_matrix(X)
  }
  if (alpha < 1 || alpha > 2) {
    stop("alpha must be in range [1,2]")
  }
  alg.typ <- match.arg(alg.typ)
  fun <- match.arg(fun)
  n <- nrow(X)
  p <- ncol(X)
  if (n.comp > min(n, p)) {
    log_message("'n.comp' is too large: reset to {min(n, p)}")
    n.comp <- min(n, p)
  }
  if (is.null(w.init)) {
    w.init <- matrix(stats::rnorm(n.comp^2), n.comp, n.comp)
  } else {
    if (!is.matrix(w.init) || length(w.init) != (n.comp^2)) {
      stop("w.init is not a matrix or is the wrong size")
    }
  }

  log_message("Centering", verbose = verbose)
  X <- scale(X, scale = FALSE)
  X <- if (row.norm) {
    Matrix::t(scale(X, scale = row.norm))
  } else {
    Matrix::t(X)
  }
  log_message("Whitening", verbose = verbose)
  V <- X %*% Matrix::t(X) / n

  log_message("Finding SVD", verbose = verbose)

  initial_v <- thisutils::as_matrix(
    stats::qnorm(1:(ncol(V) + 1) / (ncol(V) + 1))[1:ncol(V)]
  )
  s <- irlba::irlba(V, n.comp, n.comp, v = initial_v)
  svs <- s$d

  D <- diag(c(1 / sqrt(s$d)))
  K <- D %*% Matrix::t(s$u)
  K <- matrix(K[1:n.comp, ], n.comp, p)
  X1 <- K %*% X

  log_message("Running ICA", verbose = verbose)
  if (alg.typ == "deflation") {
    a <- fastICA::ica.R.def(X1, n.comp,
      tol = tol, fun = fun,
      alpha = alpha, maxit = maxit, verbose = verbose,
      w.init = w.init
    )
  } else if (alg.typ == "parallel") {
    a <- fastICA::ica.R.par(X1, n.comp,
      tol = tol, fun = fun,
      alpha = alpha, maxit = maxit, verbose = verbose,
      w.init = w.init
    )
  }
  w <- a %*% K
  S <- w %*% X
  A <- Matrix::t(w) %*% solve(w %*% Matrix::t(w))
  return(list(X = Matrix::t(X), K = Matrix::t(K), W = Matrix::t(a), A = Matrix::t(A), S = Matrix::t(S), svs = svs))
}

extract_ddrtree_ordering <- function(cds, root_cell, verbose = T) {
  dp <- cellPairwiseDistances(cds)
  dp_mst <- minSpanningTree(cds)

  curr_state <- 1

  res <- list(subtree = dp_mst, root = root_cell)

  states <- rep(1, ncol(dp))
  names(states) <- V(dp_mst)$name

  pseudotimes <- rep(0, ncol(dp))
  names(pseudotimes) <- V(dp_mst)$name

  parents <- rep(NA, ncol(dp))
  names(parents) <- V(dp_mst)$name

  if ("parent" %in% names(formals(igraph::dfs))) {
    mst_traversal <- igraph::dfs(
      dp_mst,
      root = root_cell,
      mode = "all",
      unreachable = FALSE,
      parent = TRUE
    )
    mst_parent <- as.numeric(mst_traversal$parent)
  } else {
    mst_traversal <- igraph::dfs(
      dp_mst,
      root = root_cell,
      mode = "all",
      unreachable = FALSE,
      father = TRUE
    )
    mst_parent <- as.numeric(mst_traversal$father)
  }
  curr_state <- 1

  for (i in 1:length(mst_traversal$order)) {
    curr_node <- mst_traversal$order[i]
    curr_node_name <- V(dp_mst)[curr_node]$name

    if (is.na(mst_parent[curr_node]) == FALSE) {
      parent_node <- mst_parent[curr_node]
      parent_node_name <- V(dp_mst)[parent_node]$name
      parent_node_pseudotime <- pseudotimes[parent_node_name]
      parent_node_state <- states[parent_node_name]
      curr_node_pseudotime <- parent_node_pseudotime + dp[curr_node_name, parent_node_name]
      if (igraph::degree(dp_mst, v = parent_node_name) > 2) {
        curr_state <- curr_state + 1
      }
    } else {
      parent_node <- NA
      parent_node_name <- NA
      curr_node_pseudotime <- 0
    }

    curr_node_state <- curr_state
    pseudotimes[curr_node_name] <- curr_node_pseudotime
    states[curr_node_name] <- curr_node_state
    parents[curr_node_name] <- parent_node_name
  }

  ordering_df <- data.frame(
    sample_name = names(states),
    cell_state = factor(states),
    pseudo_time = as.vector(pseudotimes),
    parent = parents
  )
  row.names(ordering_df) <- ordering_df$sample_name
  return(ordering_df)
}

select_root_cell <- function(cds, root_state = NULL, reverse = FALSE) {
  if (is.null(root_state) == FALSE) {
    if (is.null(pData(cds)$State)) {
      stop("Error: State has not yet been set. Please call orderCells() without specifying root_state, then try this call again.")
    }
    root_cell_candidates <- subset(pData(cds), State == root_state)
    if (nrow(root_cell_candidates) == 0) {
      stop(paste("Error: no cells for State =", root_state))
    }

    dp <- thisutils::as_matrix(
      stats::dist(
        Matrix::t(reducedDimS(cds)[, row.names(root_cell_candidates)])
      )
    )
    gp <- igraph::graph_from_adjacency_matrix(
      dp,
      mode = "undirected", weighted = TRUE
    )
    dp_mst <- igraph::mst(gp)

    tip_leaves <- names(which(igraph::degree(minSpanningTree(cds)) == 1))

    diameter <- igraph::get_diameter(dp_mst)

    if (length(diameter) == 0) {
      stop(paste("Error: no valid root cells for State =", root_state))
    }

    root_cell_candidates <- root_cell_candidates[names(diameter), ]
    if (is.null(cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell) == FALSE &&
      pData(cds)[cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell, ]$State == root_state) {
      root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$Pseudotime == min(root_cell_candidates$Pseudotime))]
    } else {
      root_cell <- row.names(root_cell_candidates)[which(root_cell_candidates$Pseudotime == max(root_cell_candidates$Pseudotime))]
    }
    if (length(root_cell) > 1) {
      root_cell <- root_cell[1]
    }

    if (cds@dim_reduce_type == "DDRTree") {
      graph_point_for_root_cell <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex[root_cell, ]
      root_cell <- V(minSpanningTree(cds))[graph_point_for_root_cell]$name
    }
  } else {
    if (is.null(minSpanningTree(cds))) {
      stop("Error: no spanning tree found for CellDataSet object. Please call reduceDimension before calling orderCells()")
    }
    diameter <- igraph::get_diameter(minSpanningTree(cds))
    if (is.null(reverse) == FALSE && reverse == TRUE) {
      root_cell <- names(diameter[length(diameter)])
    } else {
      root_cell <- names(diameter[1])
    }
  }
  return(root_cell)
}

#' Orders cells according to pseudotime.
#'
#' Learns a "trajectory" describing the biological process the cells are
#' going through, and calculates where each cell falls within that trajectory.
#' Monocle learns trajectories in two steps. The first step is reducing the dimensionality
#' of the data with \code{\link{reduceDimension}()}. The second is this function.
#' function. This function takes as input a CellDataSet and returns it with
#' two new columns: \code{Pseudotime} and \code{State}, which together encode
#' where each cell maps to the trajectory. \code{orderCells()} optionally takes
#' a "root" state, which you can use to specify the start of the trajectory. If
#' you don't provide a root state, one is selected arbitrarily.
#'
#' The \code{reduction_method} argument to \code{\link{reduceDimension}()}
#' determines which algorithm is used by \code{orderCells()} to learn the trajectory.
#' If \code{reduction_method == "ICA"}, this function uses \emph{polygonal reconstruction}
#' to learn the underlying trajectory. If \code{reduction_method == "DDRTree"},
#' the trajectory is specified by the principal graph learned by the
#' \code{\link[DDRTree]{DDRTree}()} function.
#'
#' Whichever algorithm you use, the trajectory will be composed of segments.
#' The cells from a segment will share the same value of \code{State}. One of
#' these segments will be selected as the root of the trajectory arbitrarily.
#' The most distal cell on that segment will be chosen as the "first" cell in the
#' trajectory, and will have a Pseudotime value of zero. \code{orderCells()} will
#' then "walk" along the trajectory, and as it encounters additional cells, it
#' will assign them increasingly large values of Pseudotime.
#'
#' @param cds the CellDataSet upon which to perform this operation
#' @param root_state The state to use as the root of the trajectory.
#' You must already have called orderCells() once to use this argument.
#' @param num_paths the number of end-point cell states to allow in the biological process.
#' @param reverse whether to reverse the beginning and end points of the learned biological process.
#'
#'
#' @return an updated CellDataSet object, in which phenoData contains values for State and Pseudotime for each cell
#' @export
orderCells <- function(
    cds,
    root_state = NULL,
    num_paths = NULL,
    reverse = NULL) {
  if (class(cds)[1] != "CellDataSet") {
    log_message(
      "Error cds is not of type 'CellDataSet'",
      message_type = "error"
    )
  }
  if (is.null(cds@dim_reduce_type)) {
    log_message(
      "Dimensionality not yet reduced. Please call reduceDimension() before calling this function",
      message_type = "error"
    )
  }
  if (any(c(length(cds@reducedDimS) == 0, length(cds@reducedDimK) == 0))) {
    log_message(
      "Dimension reduction didn't prodvide correct results. Please check your reduceDimension() step and ensure correct dimension reduction are performed before calling this function",
      message_type = "error"
    )
  }
  root_cell <- select_root_cell(cds, root_state, reverse)

  cds@auxOrderingData <- new.env(hash = TRUE)
  if (cds@dim_reduce_type == "ICA") {
    if (is.null(num_paths)) {
      num_paths <- 1
    }
    adjusted_s <- Matrix::t(cds@reducedDimS)
    dp <- thisutils::as_matrix(stats::dist(adjusted_s))
    cellPairwiseDistances(cds) <- thisutils::as_matrix(stats::dist(adjusted_s))
    gp <- igraph::graph_from_adjacency_matrix(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- igraph::mst(gp)
    minSpanningTree(cds) <- dp_mst
    next_node <<- 0
    res <- pq_helper(dp_mst, use_weights = FALSE, root_node = root_cell)

    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell

    order_list <- extract_good_branched_ordering(
      res$subtree, res$root,
      cellPairwiseDistances(cds), num_paths, FALSE
    )
    cc_ordering <- order_list$ordering_df
    row.names(cc_ordering) <- cc_ordering$sample_name

    minSpanningTree(cds) <- igraph::as_undirected(order_list$cell_ordering_tree)

    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), ]$pseudo_time
    pData(cds)$State <- cc_ordering[row.names(pData(cds)), ]$cell_state

    mst_branch_nodes <- V(minSpanningTree(cds))[which(igraph::degree(minSpanningTree(cds)) > 2)]$name

    minSpanningTree(cds) <- dp_mst
    cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree <- igraph::as_undirected(order_list$cell_ordering_tree)
  } else if (cds@dim_reduce_type == "DDRTree") {
    if (is.null(num_paths) == FALSE) {
      log_message(
        "num_paths only valid for method 'ICA' in reduceDimension()",
        message_type = "warning"
      )
    }
    cc_ordering <- extract_ddrtree_ordering(cds, root_cell)

    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), ]$pseudo_time

    K_old <- reducedDimK(cds)
    old_dp <- cellPairwiseDistances(cds)
    old_mst <- minSpanningTree(cds)
    old_A <- reducedDimA(cds)
    old_W <- reducedDimW(cds)

    cds <- project2MST(cds, project_point_to_line_segment)
    minSpanningTree(cds) <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_tree

    root_cell_idx <- which(V(old_mst)$name == root_cell, arr.ind = T)
    cells_mapped_to_graph_root <- which(
      cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == root_cell_idx
    )
    if (length(cells_mapped_to_graph_root) == 0) {
      cells_mapped_to_graph_root <- root_cell_idx
    }

    cells_mapped_to_graph_root <- V(minSpanningTree(cds))[cells_mapped_to_graph_root]$name

    tip_leaves <- names(which(igraph::degree(minSpanningTree(cds)) == 1))
    root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% tip_leaves][1]
    if (is.na(root_cell)) {
      root_cell <- select_root_cell(cds, root_state, reverse)
    }

    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell

    cc_ordering_new_pseudotime <- extract_ddrtree_ordering(cds, root_cell)

    pData(cds)$Pseudotime <- cc_ordering_new_pseudotime[row.names(pData(cds)), ]$pseudo_time
    if (is.null(root_state) == TRUE) {
      closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
      pData(cds)$State <- cc_ordering[closest_vertex[, 1], ]$cell_state
    }

    reducedDimK(cds) <- K_old
    cellPairwiseDistances(cds) <- old_dp
    minSpanningTree(cds) <- old_mst
    reducedDimA(cds) <- old_A
    reducedDimW(cds) <- old_W

    mst_branch_nodes <- V(minSpanningTree(cds))[which(igraph::degree(minSpanningTree(cds)) > 2)]$name
  } else if (cds@dim_reduce_type == "SimplePPT") {
    if (is.null(num_paths) == FALSE) {
      log_message(
        "Warning: num_paths only valid for method 'ICA' in reduceDimension()",
        message_type = "warning"
      )
    }
    cc_ordering <- extract_ddrtree_ordering(cds, root_cell)

    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), ]$pseudo_time
    pData(cds)$State <- cc_ordering[row.names(pData(cds)), ]$cell_state

    mst_branch_nodes <- V(minSpanningTree(cds))[which(igraph::degree(minSpanningTree(cds)) > 2)]$name
  }

  cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points <- mst_branch_nodes

  cds
}

normalize_expr_data <- function(
    cds,
    norm_method = c("log", "vstExprs", "none"),
    pseudo_expr = 1,
    relative_expr = TRUE) {
  FM <- exprs(cds)
  use_for_ordering <- NULL
  if (is.null(fData(cds)$use_for_ordering) == FALSE &&
    nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0) {
    FM <- FM[fData(cds)$use_for_ordering, ]
  }

  norm_method <- match.arg(norm_method)
  if (cds@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    if (is.null(pseudo_expr)) {
      if (norm_method == "log") {
        pseudo_expr <- 1
      } else {
        pseudo_expr <- 0
      }
    }

    checkSizeFactors(cds)

    if (norm_method == "vstExprs") {
      if (relative_expr == FALSE) {
        log_message("relative_expr is ignored when using norm_method == 'vstExprs'", message_type = "warning")
      }

      if (is.null(fData(cds)$use_for_ordering) == FALSE &&
        nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0) {
        VST_FM <- vstExprs(cds[fData(cds)$use_for_ordering, ], round_vals = FALSE)
      } else {
        VST_FM <- vstExprs(cds, round_vals = FALSE)
      }

      if (is.null(VST_FM) == FALSE) {
        FM <- VST_FM
      } else {
        stop("Error: set the variance-stabilized value matrix with vstExprs(cds) <- computeVarianceStabilizedValues() before calling this function with use_vst=TRUE")
      }
    } else if (norm_method == "log") {
      if (relative_expr) {
        FM <- Matrix::t(Matrix::t(FM) / BiocGenerics::sizeFactors(cds))
      }

      if (is.null(pseudo_expr)) {
        pseudo_expr <- 1
      }
      FM <- FM + pseudo_expr
      FM <- log2(FM)
    } else if (norm_method == "none") {
      FM <- Matrix::t(Matrix::t(FM) / BiocGenerics::sizeFactors(cds))
      FM <- FM + pseudo_expr
    }
  } else if (cds@expressionFamily@vfamily == "binomialff") {
    if (norm_method == "none") {
      ncounts <- FM > 0
      ncounts[ncounts != 0] <- 1
      FM <- Matrix::t(Matrix::t(ncounts) * log(1 + ncol(ncounts) / rowSums(ncounts)))
    } else {
      stop("Error: the only normalization method supported with binomial data is 'none'")
    }
  } else if (cds@expressionFamily@vfamily == "Tobit") {
    FM <- FM + pseudo_expr
    if (norm_method == "none") {

    } else if (norm_method == "log") {
      FM <- log2(FM)
    } else {
      stop("Error: the only normalization methods supported with Tobit-distributed (e.g. FPKM/TPM) data are 'log' (recommended) or 'none'")
    }
  } else if (cds@expressionFamily@vfamily == "uninormal") {
    if (norm_method == "none") {
      FM <- FM + pseudo_expr
    } else {
      stop("Error: the only normalization method supported with gaussian data is 'none'")
    }
  }
  return(FM)
}

#' Compute a projection of a CellDataSet object into a lower dimensional space
#'
#' @description Monocle aims to learn how cells transition through a biological program of
#' gene expression changes in an experiment. Each cell can be viewed as a point
#' in a high-dimensional space, where each dimension describes the expression of
#' a different gene in the genome. Identifying the program of gene expression
#' changes is equivalent to learning a \emph{trajectory} that the cells follow
#' through this space. However, the more dimensions there are in the analysis,
#' the harder the trajectory is to learn. Fortunately, many genes typically
#' co-vary with one another, and so the dimensionality of the data can be
#' reduced with a wide variety of different algorithms. Monocle provides two
#' different algorithms for dimensionality reduction via \code{reduceDimension}.
#' Both take a CellDataSet object and a number of dimensions allowed for the
#' reduced space. You can also provide a model formula indicating some variables
#' (e.g. batch ID or other technical factors) to "subtract" from the data so it
#' doesn't contribute to the trajectory.
#'
#' @details You can choose two different reduction algorithms: Independent Component
#' Analysis (ICA) and Discriminative Dimensionality Reduction with Trees (DDRTree).
#' The choice impacts numerous downstream analysis steps, including \code{\link{orderCells}}.
#' Choosing ICA will execute the ordering procedure described in Trapnell and Cacchiarelli et al.,
#' which was implemented in Monocle version 1. \code{\link[DDRTree]{DDRTree}} is a more recent manifold
#' learning algorithm developed by Qi Mao and colleages. It is substantially more
#' powerful, accurate, and robust for single-cell trajectory analysis than ICA,
#' and is now the default method.
#'
#' Often, experiments include cells from different batches or treatments. You can
#' reduce the effects of these treatments by transforming the data with a linear
#' model prior to dimensionality reduction. To do so, provide a model formula
#' through \code{residualModelFormulaStr}.
#'
#' Prior to reducing the dimensionality of the data, it usually helps
#' to normalize it so that highly expressed or highly variable genes don't
#' dominate the computation. \code{reduceDimension()} automatically transforms
#' the data in one of several ways depending on the \code{expressionFamily} of
#' the CellDataSet object. If the expressionFamily is \code{negbinomial} or \code{negbinomial.size}, the
#' data are variance-stabilized. If the expressionFamily is \code{Tobit}, the data
#' are adjusted by adding a pseudocount (of 1 by default) and then log-transformed.
#' If you don't want any transformation at all, set norm_method to "none" and
#' pseudo_expr to 0. This maybe useful for single-cell qPCR data, or data you've
#' already transformed yourself in some way.
#'
#' @param cds the CellDataSet upon which to perform this operation
#' @param max_components the dimensionality of the reduced space
#' @param reduction_method A character string specifying the algorithm to use for dimensionality reduction.
#' @param norm_method Determines how to transform expression values prior to reducing dimensionality
#' @param residualModelFormulaStr A model formula specifying the effects to subtract from the data before clustering.
#' @param pseudo_expr amount to increase expression values before dimensionality reduction
#' @param relative_expr When this argument is set to TRUE (default), we intend to convert the expression into a relative expression.
#' @param auto_param_selection when this argument is set to TRUE (default), it will automatically calculate the proper value for the ncenter (number of centroids) parameters which will be passed into DDRTree call.
#' @param scaling When this argument is set to TRUE (default), it will scale each gene before running trajectory reconstruction.
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... additional arguments to pass to the dimensionality reduction function
#' @return an updated CellDataSet object
#'
#' @export
reduceDimension <- function(
    cds,
    max_components = 2,
    reduction_method = c("DDRTree", "ICA", "tSNE", "SimplePPT", "L1-graph", "SGL-tree"),
    norm_method = c("log", "vstExprs", "none"),
    residualModelFormulaStr = NULL,
    pseudo_expr = 1,
    relative_expr = TRUE,
    auto_param_selection = TRUE,
    verbose = FALSE,
    scaling = TRUE,
    ...) {
  extra_arguments <- list(...)
  set.seed(2016)

  FM <- normalize_expr_data(cds, norm_method, pseudo_expr)

  FM <- FM[row_mean_var(FM)$var > 0, ]

  if (is.null(residualModelFormulaStr) == FALSE) {
    log_message("Removing batch effects", verbose = verbose)
    X.model_mat <- sparse.model.matrix(
      stats::as.formula(residualModelFormulaStr),
      data = pData(cds), drop.unused.levels = TRUE
    )

    fit <- limma::lmFit(FM, X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    FM <- thisutils::as_matrix(FM) - beta %*% Matrix::t(X.model_mat[, -1])
  } else {
    X.model_mat <- NULL
  }

  if (scaling) {
    FM <- thisutils::as_matrix(Matrix::t(scale(Matrix::t(FM))))
    FM <- FM[!is.na(row.names(FM)), ]
  } else {
    FM <- thisutils::as_matrix(FM)
  }

  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }

  FM <- FM[apply(FM, 1, function(x) all(is.finite(x))), ]
  if (is.function(reduction_method)) {
    reducedDim <- reduction_method(FM, ...)
    colnames(reducedDim) <- colnames(FM)
    reducedDimW(cds) <- thisutils::as_matrix(reducedDim)
    reducedDimA(cds) <- thisutils::as_matrix(reducedDim)
    reducedDimS(cds) <- thisutils::as_matrix(reducedDim)
    reducedDimK(cds) <- thisutils::as_matrix(reducedDim)
    dp <- thisutils::as_matrix(stats::dist(reducedDim))
    cellPairwiseDistances(cds) <- dp
    gp <- igraph::graph_from_adjacency_matrix(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- igraph::mst(gp)
    minSpanningTree(cds) <- dp_mst
    cds@dim_reduce_type <- "function_passed"
  } else {
    reduction_method <- match.arg(reduction_method)
    if (reduction_method == "tSNE") {
      log_message("Remove noise by PCA ...", verbose = verbose)

      if ("num_dim" %in% names(extra_arguments)) {
        num_dim <- extra_arguments$num_dim
      } else {
        num_dim <- 50
      }

      FM <- (FM)
      irlba_res <- prcomp_irlba(Matrix::t(FM),
        n = min(num_dim, min(dim(FM)) - 1),
        center = TRUE, scale. = TRUE
      )
      irlba_pca_res <- irlba_res$x

      topDim_pca <- irlba_pca_res

      log_message("Reduce dimension by tSNE ...", verbose = verbose)

      tsne_res <- Rtsne::Rtsne(thisutils::as_matrix(topDim_pca), dims = max_components, pca = F, ...)

      tsne_data <- tsne_res$Y[, 1:max_components]
      row.names(tsne_data) <- colnames(tsne_data)

      reducedDimA(cds) <- Matrix::t(tsne_data)

      cds@auxClusteringData[["tSNE"]]$pca_components_used <- num_dim
      cds@auxClusteringData[["tSNE"]]$reduced_dimension <- Matrix::t(topDim_pca)

      cds@dim_reduce_type <- "tSNE"
    } else if (reduction_method == "ICA") {
      log_message("Reducing to independent components", verbose = verbose)
      init_ICA <- ica_helper(Matrix::t(FM), max_components,
        use_irlba = TRUE, ...
      )
      x_pca <- Matrix::t(Matrix::t(FM) %*% init_ICA$K)
      W <- Matrix::t(init_ICA$W)
      weights <- W
      A <- Matrix::t(solve(weights) %*% Matrix::t(init_ICA$K))
      colnames(A) <- colnames(weights)
      rownames(A) <- rownames(FM)
      S <- weights %*% x_pca
      rownames(S) <- colnames(weights)
      colnames(S) <- colnames(FM)
      reducedDimW(cds) <- thisutils::as_matrix(W)
      reducedDimA(cds) <- thisutils::as_matrix(A)
      reducedDimS(cds) <- thisutils::as_matrix(S)
      reducedDimK(cds) <- thisutils::as_matrix(init_ICA$K)
      adjusted_s <- Matrix::t(reducedDimS(cds))
      dp <- thisutils::as_matrix(stats::dist(adjusted_s))
      cellPairwiseDistances(cds) <- dp
      gp <- igraph::graph_from_adjacency_matrix(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- igraph::mst(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "ICA"
    } else if (reduction_method == "DDRTree") {
      log_message("Learning principal graph with DDRTree", verbose = verbose)

      if (auto_param_selection & ncol(cds) >= 100) {
        if ("ncenter" %in% names(extra_arguments)) {
          ncenter <- extra_arguments$ncenter
        } else {
          ncenter <- cal_ncenter(ncol(FM))
        }
        ddr_args <- c(
          list(X = FM, dimensions = max_components, ncenter = ncenter, verbose = verbose),
          extra_arguments[names(extra_arguments) %in% c("initial_method", "maxIter", "sigma", "lambda", "param.gamma", "tol")]
        )
        ddrtree_res <- do.call(DDRTree, ddr_args)
      } else {
        ddrtree_res <- DDRTree(FM, max_components, verbose = verbose, ...)
      }
      if (ncol(ddrtree_res$Y) == ncol(cds)) {
        colnames(ddrtree_res$Y) <- colnames(FM)
      } else {
        colnames(ddrtree_res$Y) <- paste("Y_", 1:ncol(ddrtree_res$Y), sep = "")
      }

      colnames(ddrtree_res$Z) <- colnames(FM)
      reducedDimW(cds) <- ddrtree_res$W
      reducedDimS(cds) <- ddrtree_res$Z
      reducedDimK(cds) <- ddrtree_res$Y
      cds@auxOrderingData[["DDRTree"]]$objective_vals <- ddrtree_res$objective_vals

      adjusted_K <- Matrix::t(reducedDimK(cds))
      dp <- thisutils::as_matrix(stats::dist(adjusted_K))
      cellPairwiseDistances(cds) <- dp
      gp <- igraph::graph_from_adjacency_matrix(dp, mode = "undirected", weighted = TRUE)
      dp_mst <- igraph::mst(gp)
      minSpanningTree(cds) <- dp_mst
      cds@dim_reduce_type <- "DDRTree"
      cds <- findNearestPointOnMST(cds)
    } else {
      stop("Error: unrecognized dimensionality reduction method")
    }
  }
  cds
}

findNearestPointOnMST <- function(cds) {
  dp_mst <- minSpanningTree(cds)
  Z <- reducedDimS(cds)
  Y <- reducedDimK(cds)

  tip_leaves <- names(which(igraph::degree(dp_mst) == 1))

  distances_Z_to_Y <- proxy::dist(Matrix::t(Z), Matrix::t(Y))
  closest_vertex <- apply(distances_Z_to_Y, 1, function(z) {
    which(z == min(z))[1]
  })

  closest_vertex_names <- colnames(Y)[closest_vertex]
  closest_vertex_df <- thisutils::as_matrix(closest_vertex)
  row.names(closest_vertex_df) <- names(closest_vertex)

  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
  cds
}

project2MST <- function(cds, Projection_Method) {
  dp_mst <- minSpanningTree(cds)
  Z <- reducedDimS(cds)
  Y <- reducedDimK(cds)

  cds <- findNearestPointOnMST(cds)
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex

  closest_vertex_names <- colnames(Y)[closest_vertex]
  closest_vertex_df <- thisutils::as_matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)

  tip_leaves <- names(which(igraph::degree(dp_mst) == 1))

  if (!is.function(Projection_Method)) {
    P <- Y[, closest_vertex]
  } else {
    P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
    for (i in 1:length(closest_vertex)) {
      neighbors <- names(igraph::neighbors(dp_mst, closest_vertex_names[i], mode = "all"))
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]

      for (neighbor in neighbors) {
        if (closest_vertex_names[i] %in% tip_leaves) {
          tmp <- projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
        } else {
          tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
        }
        projection <- rbind(projection, tmp)
        distance <- c(distance, stats::dist(rbind(Z_i, tmp)))
      }
      if (!inherits(projection, "matrix")) {
        projection <- thisutils::as_matrix(projection)
      }
      P[, i] <- projection[which(distance == min(distance))[1], ]
    }
  }

  colnames(P) <- colnames(Z)

  dp <- thisutils::as_matrix(stats::dist(Matrix::t(P)))

  min_dist <- min(dp[dp != 0])
  dp <- dp + min_dist
  diag(dp) <- 0

  cellPairwiseDistances(cds) <- dp
  gp <- igraph::graph_from_adjacency_matrix(
    dp,
    mode = "undirected", weighted = TRUE
  )
  dp_mst <- igraph::mst(gp)

  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df

  cds
}

projPointOnLine <- function(point, line) {
  ap <- point - line[, 1]
  ab <- line[, 2] - line[, 1]

  res <- line[, 1] + c((ap %*% ab) / (ab %*% ab)) * ab
  return(res)
}

project_point_to_line_segment <- function(p, df) {
  A <- df[, 1]
  B <- df[, 2]
  AB <- (B - A)
  AB_squared <- sum(AB^2)
  if (AB_squared == 0) {
    q <- A
  } else {
    Ap <- (p - A)
    t <- sum(Ap * AB) / AB_squared

    if (t < 0.0) {
      q <- A
    } else if (t > 1.0) {
      q <- B
    } else {
      q <- A + t * AB
    }
  }
  return(q)
}

traverseTree <- function(g, starting_cell, end_cells) {
  distance <- igraph::distances(g, v = starting_cell, to = end_cells)
  branchPoints <- which(igraph::degree(g) == 3)
  path <- igraph::shortest_paths(g, from = starting_cell, end_cells)

  return(list(shortest_path = path$vpath, distance = distance, branch_points = intersect(branchPoints, unlist(path$vpath))))
}

traverseTreeCDS <- function(cds, starting_cell, end_cells) {
  subset_cell <- c()
  dp_mst <- cds@minSpanningTree

  for (end_cell in end_cells) {
    traverse_res <- traverseTree(dp_mst, starting_cell, end_cell)
    path_cells <- names(traverse_res$shortest_path[[1]])

    subset_cell <- c(subset_cell, path_cells)
  }

  subset_cell <- unique(subset_cell)
  cds_subset <- SubSet_cds(cds, subset_cell)

  root_state <- pData(cds_subset[, starting_cell])[, "State"]
  cds_subset <- orderCells(cds_subset, root_state = as.numeric(root_state))

  return(cds_subset)
}

SubSet_cds <- function(cds, cells) {
  cells <- unique(cells)
  if (ncol(reducedDimK(cds)) != ncol(cds)) {
    stop("SubSet_cds doesn't support cds with ncenter run for now. You can try to subset the data and do the construction of trajectory on the subset cds")
  }

  exprs_mat <- as(thisutils::as_matrix(cds[, cells]), "sparseMatrix")
  cds_subset <- newCellDataSet(exprs_mat,
    phenoData = new("AnnotatedDataFrame", data = pData(cds)[colnames(exprs_mat), ]),
    featureData = new("AnnotatedDataFrame", data = fData(cds)),
    expressionFamily = negbinomial.size(),
    lowerDetectionLimit = 1
  )
  BiocGenerics::sizeFactors(cds_subset) <- BiocGenerics::sizeFactors(cds[, cells])
  cds_subset@dispFitInfo <- cds@dispFitInfo

  cds_subset@reducedDimW <- cds@reducedDimW
  cds_subset@reducedDimS <- cds@reducedDimS[, cells]
  cds_subset@reducedDimK <- cds@reducedDimK[, cells]

  cds_subset@cellPairwiseDistances <- cds@cellPairwiseDistances[cells, cells]

  adjusted_K <- Matrix::t(reducedDimK(cds_subset))
  dp <- thisutils::as_matrix(stats::dist(adjusted_K))
  cellPairwiseDistances(cds_subset) <- dp
  gp <- igraph::graph_from_adjacency_matrix(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- igraph::mst(gp)
  minSpanningTree(cds_subset) <- dp_mst
  cds_subset@dim_reduce_type <- "DDRTree"
  cds_subset <- findNearestPointOnMST(cds_subset)

  cds_subset <- orderCells(cds_subset)
}

reverseEmbeddingCDS <- function(cds) {
  if (nrow(cds@reducedDimW) < 1) {
    stop("You need to first apply reduceDimension function on your cds before the reverse embedding")
  }

  FM <- normalize_expr_data(cds, norm_method = "log")

  FM <- FM[row_mean_var(FM)$var > 0, ]

  reverse_embedding_data <- reducedDimW(cds) %*% reducedDimS(cds)
  row.names(reverse_embedding_data) <- row.names(FM)

  cds_subset <- cds[row.names(FM), ]

  reverse_embedding_data <- Matrix::t(apply(reverse_embedding_data, 1, function(x) x + abs(min(x))))

  raw_data <- thisutils::as_matrix(exprs(cds)[row.names(FM), ])
  reverse_embedding_data <- reverse_embedding_data * (apply(raw_data, 1, function(x) stats::quantile(x, 0.99))) / apply(reverse_embedding_data, 1, max)

  Biobase::exprs(cds_subset) <- reverse_embedding_data

  return(cds_subset)
}

cal_ncenter <- function(ncells, ncells_limit = 100) {
  round(2 * ncells_limit * log(ncells) / (log(ncells) + log(ncells_limit)))
}
