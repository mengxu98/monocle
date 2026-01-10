#' @import Biobase
#' @import plyr
#' @import ggplot2
#' @import Matrix
#' @import methods
#' @import slam
#' @import irlba
#' @import DDRTree
#' @import Rtsne
#' @import grid
#' @import pheatmap
#' @import VGAM
#' @import slam
#' @import HSMMSingleCell

#' @importFrom stats median
#' @importFrom Rcpp sourceCpp
#' @importFrom igraph V E V<-
#' @importFrom dplyr %>%
#' @importFrom BiocGenerics sizeFactors<- estimateSizeFactors estimateDispersions
#' @importFrom thisutils log_message

utils::globalVariables(
  c(
    "rowname", "Size_Factor", "next_node",
    "use_for_ordering", "from", "to",
    "pseudocount", "Branch", "CellType", "Pseudotime",
    "ids", "prin_graph_dim_1", "prin_graph_dim_2", "State",
    "feature_label", "expectation", "colInd", "rowInd", "value",
    "source_prin_graph_dim_1", "source_prin_graph_dim_2",
    "component"
  )
)
