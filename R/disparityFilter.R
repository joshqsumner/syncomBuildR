#' Disparity Filtering for microbial networks
#'
#' @param net Object returned from \link{asvNet}.
#' @param weights Optional weights specification. This defaults to NULL in which case "weight" is used
#' as calculated by igraph. This can also be a numeric vector of edge weights
#' (ordered for \code{net$edges}) or a column name from \code{net$edges}.
#' @param alpha the significance level to filter edges for. Defaults to 0.05 for no serious reason.
#' @param cores Number of cores to run in parallel, defaults to 1 if "mc.cores" option is unset.
#' @param replot Logical, should nodes be rearranged to represent the network better visually?
#' Defaults to TRUE.
#' @importFrom stats quantile
#' @return A modified version of net with filtered edges (and nodes if any were now isolated).
#'
#' @examples
#'
#' taxa <- c(
#'   "Bacteria", "Proteobacteria", "Betaproteobacteria", "Burkholderiales",
#'   "Burkholderiaceae", "Paraburkholderia", NA
#' )
#' taxa <- matrix(rep(taxa, 10), nrow = 10, byrow = TRUE)
#' colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#' rownames(taxa) <- paste0("ASV", 1:10)
#' # taxonomy data if used should have ASV names explicitly as a column
#' taxa_df <- as.data.frame(taxa)
#' taxa_df$asv <- rownames(taxa_df)
#'
#' sp_dist <- asvDist(asv, method = "spearman", clr_transform = TRUE, edgeFilter = 0.5)
#' net_data <- asvNet(sp_dist, taxa_df, edge = "spearman_similarity")
#' # we could use spearman_similarity, but not specifying will use the graph's default weights
#' net_data_disp_filt <- dispFilter(net_data, alpha = 0.4)
#'
#' @export
#'

dispFilter <- function(net, weights = NULL, alpha = 0.05, cores = getOption("mc.cores", 1),
                       replot = TRUE) {
  # these could be arguments, but I don't see a reason to change them ever
  node_metric <- "degree"
  node_metric_minimum <- 1
  # calculated values
  weights <- .check_weights(net, weights)
  net_met <- .check_node_metric(net, node_metric)
  net$edges$disp_filt_weight <- as.numeric(weights)
  # grab edge data for convenience
  edges <- net$edges
  edges$p_value <- NA
  # get ids of nodes to loop over
  all_nodes <- igraph::as_ids(igraph::V(net$graph))
  use_node <- net_met > node_metric_minimum
  # apply disparity filtering over nodes
  edges <- do.call(rbind, parallel::mclapply(seq_along(all_nodes), function(i) {
    nd <- all_nodes[i]
    edges_sub <- edges[which(edges$from == nd | edges$to == nd), ]
    if (!use_node[i]) {
      return(edges_sub)
    }
    w <- sum(as.numeric(edges_sub$disp_filt_weight)) # total weight coming into this node
    k <- net_met[i] # degree for this node
    edges_sub <- do.call(
      rbind,
      lapply(
        igraph::as_ids(
          igraph::ego(net$graph, order = 1, nodes = nd)[[1]][-1]
        ),
        function(v) {
          j <- which(edges_sub$from == nd & edges_sub$to == v)
          edges_sub$p_value[j] <- (1 - edges_sub$disp_filt_weight[j] / w)^(k - 1)
          return(edges_sub[j, ])
        }
      )
    )
    return(edges_sub)
  }, mc.cores = cores))
  # invert p-values so that edgeFilter can be used
  edges$inverted_p <- 1 - edges$p_value
  edges$inverted_p <- ifelse(is.na(edges$inverted_p), 0, edges$inverted_p)
  net$edges <- edges
  # call edgeFilter
  net <- edgeFilter(net, filter = 1 - alpha, edge = "inverted_p", replot = replot)
  net$edges <- net$edges[, -which(colnames(net$edges) %in% c("disp_filt_weight", "inverted_p"))]
  return(net)
}

#' @keywords internal
#' @noRd

.check_node_metric <- function(net_data, node_metric) {
  if (is.null(node_metric)) {
    node_metric <- igraph::degree(net_data$graph)
  } else if (length(node_metric) == 1 && is.character(node_metric)) {
    node_metric <- net_data$nodes[[node_metric]]
  } else if (!is.numeric(node_metric)) {
    stop("node_metric must be a numeric vector, NULL, or a single column name from the nodes data")
  }
  return(node_metric)
}

#' @keywords internal
#' @noRd
#'
.check_weights <- function(net_data, weights) {
  if (is.null(weights)) {
    weights <- igraph::E(net_data$graph)$weight
  } else if (length(weights) == 1 && is.character(weights)) {
    weights <- as.numeric(net_data$edges[[weights]])
  } else if (!is.numeric(weights)) {
    stop("weights must be a numeric vector, NULL, or a single column name from the edges data")
  }
  return(weights)
}
