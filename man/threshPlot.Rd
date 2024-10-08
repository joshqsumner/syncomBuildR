% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/threshPlot.R
\name{threshPlot}
\alias{threshPlot}
\title{Function to plot results of a changepoint model with data returned from \code{thresh}}
\usage{
threshPlot(
  thresh,
  asv,
  asvCols = NULL,
  phenotype = NULL,
  unit = "asv",
  net = NULL
)
}
\arguments{
\item{thresh}{Output from \code{\link{thresh}}}

\item{asv}{The asv table used in making \code{thresh} object.}

\item{asvCols}{A vector of ASV column names. Defaults to NULL in which case all columns containing
"ASV" are used and a list of ggplots is returned.}

\item{phenotype}{A vector of phenotype names in \code{thresh}. Defaults to NULL where all phenotypes
are used and a list of plots is returned per ASV.}

\item{unit}{The unit or scale of the changepoint models. This defaults to "asv" for use with thresh
and "cluster" should be used with netThresh output.}

\item{net}{The asvNet object if netThresh output is being plotted.}
}
\value{
A ggplot or list of ggplots showing changepoint models against some set of phenotypes.
}
\description{
Function to plot results of a changepoint model with data returned from \code{thresh}
}
\examples{

asv$biomass_z <- rnorm(nrow(asv))
tm <- thresh(asv, "biomass_z")
threshPlot(tm, asv, "ASV9")

}
\keyword{changepoint,}
\keyword{ggplot}
\keyword{phenotype,}
\keyword{regression,}
\keyword{threshold,}
