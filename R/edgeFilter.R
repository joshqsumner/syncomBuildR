#' Function to filter edges from networks generated by \code{asvNet}.
#'
#' @param net Object returned from \link{asvNet}.
#' @param filter Value to filter edges for. If non-NULL then only edges with edgeWeight greater than
#' this value are kept. This can be a character vector or a numeric.
#' @param edge Optional weighting for edges. Must be present in the "edges" of net. Default of NULL
#' will show equal size edges between all connected nodes.
#' @importFrom stats quantile
#' @return A modified version of net with filtered edges (and nodes if any were now isolated).
#'
#' @examples
#'
#' # a<-qc(); b<-cal(a); c<-thresh(b); d<-asvDist(a) ; e<-net(d, thresh = c)
#' print(load("~/scripts/SINC/sincUtils/syncomBuilder/net_output.rdata"))
#' dim(net_data$edges)
#' net_data2 <- edgeFilter(net_data, 0.6)
#' dim(net_data2$edges)
#' net_data3 <- edgeFilter(net_data, "0.6")
#' dim(net_data3$edges)
#'
#' @export
#'

edgeFilter <- function(net, filter, edge = "spearman") {
  original_nodes <- net[["nodes"]]
  original_edges <- net[["edges"]]
  if (is.character(filter)) {
    cutoff <- stats::quantile(original_edges[[edge]], probs = as.numeric(filter))
    edges <- original_edges[original_edges[[edge]] >= as.numeric(cutoff), ]
    removed_edges <- original_edges[original_edges[[edge]] < as.numeric(cutoff), ]
  } else if (is.numeric(filter)) {
    edges <- original_edges[original_edges[[edge]] >= filter, ]
    removed_edges <- original_edges[original_edges[[edge]] < filter, ]
  }
  removed_edge_names <- paste(removed_edges$from, removed_edges$to, sep = "|")
  nodes <- original_nodes[original_nodes$asv %in% unique(c(edges$from, edges$to)), ]
  net[["edges"]] <- edges
  net[["nodes"]] <- nodes
  removed_nodes <- setdiff(original_nodes$asv, nodes$asv)
  net[["graph"]] <- igraph::delete_edges(net[["graph"]], removed_edge_names)
  net[["graph"]] <- igraph::delete_vertices(net[["graph"]], removed_nodes)
  return(net)
}
