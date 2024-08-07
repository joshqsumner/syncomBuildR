% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcadf.R
\name{pcadf}
\alias{pcadf}
\title{Function to run a PCA, plot and return the data with PC coordinates}
\usage{
pcadf(
  df = NULL,
  cols = NULL,
  color = NULL,
  returnData = TRUE,
  ncp = NULL,
  umap = FALSE
)
}
\arguments{
\item{df}{Dataframe to ordinate}

\item{cols}{columns to reduce dimensions of. Can be specified with names or positions. Defaults to
all column names containing "ASV".}

\item{color}{column name used to color points in the pca plot.}

\item{returnData}{Logical, should data be returned? Defaults to TRUE where data and a ggplot are
eturned.}

\item{ncp}{Optional, number of principal components to return attached to dataframe if data is
returned. Defaults to all.}

\item{umap}{Logical, should a UMAP also be performed? Defaults to FALSE. If TRUE then UAMP1 and 2
columns will be added and a list of 2 ggplots will be returned.}
}
\value{
a ggplot, optionally a dataframe and FactoMineR::PCA output if returnData is TRUE.
}
\description{
Function to run a PCA, plot and return the data with PC coordinates
}
\examples{
print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/cal_output.rdata"))
asv <- are_c[[1]]
zinbCalibrated <- are_c[[2]][are_c[[2]]$model == "ZINB", "asv"]
asv <- are_c[[1]][, c("tissue", "plot", "row", "genotype", "biomass", "sd", zinbCalibrated)]
x <- pcadf(df = asv, cols = NULL, color = c("tissue", "genotype"), returnData = T, ncp = NULL)
x$plot

}
\keyword{pca}
\keyword{umap}
