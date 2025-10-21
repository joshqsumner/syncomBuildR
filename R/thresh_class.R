#' Class \code{thresh} for output from \code{syncomBuildR::thresh} and related functions.
#'
#' Class for changepoint models as used in syncomBuildR.
#'
#' @name thresh-class
#' @docType class
#'
#' @details
#' See \code{methods(class = "thresh")} for an overview of available methods.
#'
#' @slot intercept The intercept term (y value pre changepoint, per hinge model)
#' @slot changepoint The changepoint term (x value at which slope takes effect)
#' @slot slope The slope term
#' @slot phenotype The phenotype (outcome) that the model was fit to.
#' @slot model the model object (optionally)
#' @slot predictor the ASV/Cluster (predictor) that the model was fit using.
#' @slot data The data used to fit the model
#' @slot type The backend used to fit the model as well as the type of model. Currently "chngptm" is
#' supported.
#' @slot unit The unit/scope. Current options are "individual" for single ASVs or "cluster" for groups
#' of ASVs.
#' @slot control Additional information such as the call.
#'
#' @seealso
#'   \code{\link{thresh}}
#'
NULL

as.thresh <- function(x) {
  class(x) <- "thresh"
  return(x)
}

#' Print an \code{thresh} object.
#'
#' @aliases print.thresh
#'
#' @param x An object of class \code{thresh}.
#' @param ... further arguments, passed to print.default.
#'
#' @seealso \code{\link{summary.thresh}}
#' @method print thresh
#' @export
print.thresh <- function(x, ...) {
  return(summary.thresh(x, ...))
}


#' Summarize an \code{thresh} object.
#'
#' @aliases summary.thresh
#'
#' @param object An object of class \code{thresh}.
#' @param ... further arguments, passed to summary.default
#'
#' @method summary thresh
#' @importFrom utils head tail
#' @export

summary.thresh <- function(object, ...) {
  #* `Main Points`
  has_adj <- is.null(object$control$p.adjust.method) || object$control$p.adjust.method == "none"
  n_significant <- sum(
    object$pval < 0.05
  )
  if (n_significant == 0) {
    n_significant <- sum(
      object$post.prob > 0.95
    )
  }
  cat(
    paste0(
      "Thresh fit with ", object$type, " to ", object$unit, "s\n",
      "\tIncludes ", length(unique(object$predictor)), " predictors over ",
      length(unique(object$phenotype)), " phenotypes using ", nrow(object$data), " observations\n",
      "\t", n_significant, " significant slopes",
      ifelse(has_adj,
        "",
        paste0(" (", object$control$p.adjust.method, " adjusted)")
      )
    )
  )
  n_phenos <- length(unique(object$phenotype))
  phenos <- utils::head(unique(object$phenotype), 2)
  for (p in phenos) {
    cat(paste0("\n\n", p, "\n"))
    sub <- object[object$phenotype == p]
    if (length(sub$intercept) > 1 && object$type == "chngptm") {
      sum_mat <- cbind(
        sapply(
          sub[c("intercept", "changepoint")],
          function(x) {
            summary(unlist(x))
          }
        ),
        summary(sub$slope),
        summary(sub$pval)
      )
      colnames(sum_mat)[3:4] <- c("slope", "pval")
      print(sum_mat)
    } else if (length(sub$post_means) > 1) {
      sum_mat <- apply(do.call(rbind, sub$post_means), 2, summary)
      sum_mat <- cbind(sum_mat, summary(sub$post_probs))
      colnames(sum_mat)[ncol(sum_mat)] <- "post_prob"
      print(sum_mat)
    }
  }
  if (n_phenos > 2) {
    cat("\n...")
  }
  return(invisible(object))
}

#' Bracket subsetting for thresh objects
#' @aliases `[.thresh`
#'
#' Thresh objects have variable lengths and it's common to want to subset them based on
#' phenotypes or predictors. Specialized bracket methods are implemented to support that.
#' @param x thresh object
#' @param i indexing, logical and numeric vectors are handled specially to subset the thresh object.
#' Other options are handled per a normal list
#' @examples
#' asv$biomass_z <- rnorm(nrow(asv))
#' tm <- thresh(asv, "biomass_z")
#' tm["predictor"] # pulls list
#' tm[1:2] # subsets lists to make a smaller thresh object
#' @method [ thresh
#' @export

"[.thresh" <- function(x, i) {
  if (!is.logical(i) && !is.numeric(i)) {
    return(NextMethod())
  }
  # apply subsetting to subsettable slots
  x[x$control$subsettable] <- lapply(x[x$control$subsettable], function(j) {
    j[i]
  })
  # other slots don't change
  return(x)
}


#' Getting intersections of thresh objects
#'
#' @param x thresh object
#' @param y thresh object
#' @param index.return Logical, should the index be returned as an element of the output list
#' (TRUE, default) or should only subset thresh objects be returned (FALSE)?
#'
#' @export

intersect_thresh <- function(x, y, index.return = TRUE) {
  x_ids <- mapply(
    function(i, thresh) {
      paste(unlist(thresh[i][c("phenotype", "predictor")]), collapse = "_")
    },
    seq_along(x$phenotype),
    MoreArgs = list(thresh = x),
    SIMPLIFY = FALSE
  )
  y_ids <- mapply(
    function(i, thresh) {
      paste(unlist(thresh[i][c("phenotype", "predictor")]), collapse = "_")
    },
    seq_along(y$phenotype),
    MoreArgs = list(thresh = y),
    SIMPLIFY = FALSE
  )
  conserved_ids <- intersect(x_ids, y_ids)
  x_index <- which(x_ids %in% conserved_ids)
  y_index <- which(y_ids %in% conserved_ids)
  # a little weird to return a list, but I don't have a clearer idea yet about what the
  out <- list( # output could be other than a list. Normally a unique vector for intersect, but since
    x[x_index],  # the thresh objects have different slopes/p-values/changepoints that we do want
    y[y_index] # to keep I don't think that would make a lot of sense.
  )
  if (index.return) {
    out$ix <- list(x = x_index, y = y_index)
  }
  return(out)
}
