#' Function to make distance/dissimilarity data from an ASV table.
#'
#'
#' @param asvTab ASV table to use. Should have rows as samples and columns representing ASVs.
#' @param asvCols Either column names, positions, or a single character string to use as regex to
#' find all asv columns.
#' @param method A distance method compatible with \code{vegan::vegdist}, a correlation
#' coefficient for use with \code{Hmisc::rcorr} (spearman or pearson), or one of "blomqvist" or
#' "gaussrank", which are implemented here in R and therefore will be slower.
#' @param parallel How many cores to use in parallel. The parallelized component only comes into play
#' with correlation coefficient methods and is not particularly heavy, so more cores will have a
#' minimal speed improvement.
#' @param clr_transform Should ASV counts be transformed with a centered log ratio?
#' This defaults to TRUE to alleviate the compositional nature of the data.
#' @param edgeFilter A filter for the edges returned. This defaults to NULL where no filtering is done.
#' Optionally this can be numeric, in which case distances/correlations below that number are removed,
#' or it can be a number as a string in which case that quantile is used
#' (ie, "0.5" filters for above the median.).
#' @param plot Logical, should the data be plotted? Defaults to FALSE.
#' @param returnASVs Should ASVs or samples be considered the experimental unit?
#' Defaults to TRUE in which case distances/correlations between ASVs across samples are returned.
#' If FALSE then samples are compared by their ASV composition. Note that different methods are
#' appropriate depending on whether samples or ASVs are being considered "nodes".
#' @keywords asv
#' @importFrom vegan vegdist
#' @importFrom parallel mclapply
#' @import data.table
#' @import ggplot2
#' @importFrom stats setNames median qnorm
#' @importFrom viridis scale_fill_viridis
#' @importFrom Hmisc rcorr
#' @importFrom compositions clr
#' @return A dataframe showing pairwise correlations between individual ASVs/samples. If plot is TRUE
#' then a list is returned with the dataframe as the "data" element and the plot as the "plot" element.
#'
#' @examples
#'
#' sp_dist <- asvDist(asv, method = "spearman", clr_transform = TRUE, edgeFilter = 0.5)
#'
#' # note how distance is calculated:
#' x <- seq(-1, 1, 0.0001) # correlation range
#' y <- sqrt(2 * (1 - x)) # cosine transformation to euclidean distance
#' y2 <- 1 / sqrt(2 * (1 - x)) # final transformation to similarity
#' plot(x, y, type = "l")
#' plot(x, y2, type = "l")
#'
#' @details
#' If "spearman" or "pearson" methods are used then euclidean distance is also calculated from those
#' correlations using the cosine theorem. That distance is bounded on [0, 2] corresponding to
#' correlations of -1 and 1, respectively. That distance is then used to calculate similarity for
#' consistent network edges. See examples.
#'
#' @export

asvDist <- function(asvTab, asvCols = NULL, method = "spearman",
                    parallel = getOption("mc.cores", 1), clr_transform = TRUE,
                    edgeFilter = NULL, plot = FALSE, returnASVs = TRUE) {
  #* `calculated values`

  if (is.null(asvCols)) {
    asvCols <- "ASV"
  }
  if (length(asvCols) == 1) {
    asvCols <- colnames(asvTab)[grepl(asvCols, colnames(asvTab))]
  }
  vegan_distances <- c(
    "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower",
    "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis",
    "chisq", "chord", "aitchison", "robust.aitchison"
  )
  correlation_coefficients <- c("spearman", "pearson")
  scb_correlations <- c("blomqvist", "gaussrank")

  #* *pull asvs as matrix* [same for all options]
  if (returnASVs) {
    mat <- as.matrix(t(asvTab[, asvCols]))
  } else {
    mat <- as.matrix(asvTab[, asvCols])
  }
  old_dim <- dim(mat)
  mat <- mat[rowSums(mat) != 0, colSums(mat) != 0]
  new_dim <- dim(mat)
  if (any(old_dim != new_dim)) {
    warning("Rows and Columns that summed to 0 have been removed.")
  }

  #* *apply transformations* [CLR makes most sense, but there could be others in the future]

  if (clr_transform) {
    mat <- compositions::clr(mat)
  }

  #* *make long data* [from similarity matrices]

  if (method %in% vegan_distances) {
    dist_mat <- as.matrix(vegan::vegdist(mat, method = method))
    sim_mat <- as.data.frame(1 / dist_mat)
    sim_mat$c1 <- rownames(sim_mat)
    ldf <- data.table::melt(data.table::as.data.table(sim_mat),
      id.vars = c("c1"),
      variable.name = "c2", value.name = method
    )
  } else if (method %in% correlation_coefficients) {
    M <- Hmisc::rcorr(t(mat), type = method)
    M[[paste0(method, "_distance")]] <- sqrt(2 * (1 - M[["r"]])) # turn correlation into a distance
    M[[paste0(method, "_similarity")]] <- 1 / sqrt(2 * (1 - M[["r"]])) # turn distance into similarity
    method <- paste0(method, "_similarity")
    ldf <- do.call(rbind, parallel::mclapply(which(unlist(lapply(M, methods::is, "matrix"))),
                                             function(m) {
      x <- as.data.frame(M[[m]])
      x$rowname <- rownames(x)
      x$trait <- names(M)[[m]]
      data.table::melt(data.table::as.data.table(x), id.vars = c("rowname", "trait"))
    }, mc.cores = parallel))
    ldf <- as.data.frame(data.table::dcast(ldf, ... ~ trait))
    colnames(ldf)[1:2] <- c("c1", "c2")
    ldf$c2 <- as.character(ldf$c2)
  } else if (method %in% scb_correlations) {
    matched_fun <- get(paste0(".scb_", method))
    ldf <- do.call(rbind, parallel::mclapply(asvCols, function(ac1) {
      remaining_asv_cols <- c(ac1, asvCols[-c(1:which(asvCols == ac1))])
      do.call(rbind, lapply(remaining_asv_cols, function(ac2) {
        data.frame(c1 = ac1, c2 = ac2, d = matched_fun(
          as.numeric(mat[ac1, ]),
          as.numeric(mat[ac2, ])
        ))
      }))
    }))
    ldf <- stats::setNames(ldf, c("c1", "c2", method))
  }

  #* *plotting*
  if (plot) {
    # consider grouping things here maybe with hclust.
    p <- ggplot2::ggplot(ldf, ggplot2::aes(
      x = .data[["c1"]],
      y = .data[["c2"]],
      fill = .data[[method]]
    )) +
      ggplot2::geom_tile(color = NA, linewidth = 0) +
      viridis::scale_fill_viridis(na.value = "grey100") +
      ggplot2::labs(fill = method) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank()
      ) +
      ggplot2::theme_minimal()
  }

  #* *filter long data* [based on strength or P value or both]
  #* First remove NA/Inf values since those are bad for us

  ldf[[method]][is.na(ldf[[method]]) | is.infinite(ldf[[method]])] <- NA
  ldf <- ldf[!is.na(ldf[[method]]), ]
  if (!is.null(edgeFilter)) {
    if (is.character(edgeFilter)) {
      cutoff <- quantile(ldf[[method]], probs = as.numeric(edgeFilter))
      ldf <- ldf[ldf[[method]] >= as.numeric(cutoff), ]
    } else if (is.numeric(edgeFilter)) {
      ldf <- ldf[ldf[[method]] >= edgeFilter, ]
    } else {
      stop("edgeFilter must be character or numeric, see ?asvDist for details.")
    }
  }
  if (plot) {
    ldf <- list("data" = ldf, "plot" = p)
  }
  #* *return data*
  return(ldf)
}

#' @keywords internal
#' @noRd

.scb_blomqvist <- function(x, y) {
  med_x <- stats::median(x, na.rm = TRUE)
  med_y <- stats::median(y, na.rm = TRUE)
  x1 <- sign(x - med_x)
  y1 <- sign(y - med_y)
  out <- cor(x1, y1, method = "pearson")
  return(out)
}

#' @keywords internal
#' @noRd

.scb_gaussrank <- function(x, y) {
  out <- cor(stats::qnorm(rank(x) / (length(x) + 1)),
    stats::qnorm(rank(y) / (length(y) + 1)),
    method = "pearson"
  )
  return(out)
}
