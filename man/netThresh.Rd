% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netThresh.R
\name{netThresh}
\alias{netThresh}
\title{Function to aggregate clustered networks overall phenotypic impacts.}
\usage{
netThresh(
  net,
  asvTab,
  asvCols = NULL,
  clusterCol = NULL,
  cluster = NULL,
  phenoCols = NULL,
  model = "hinge",
  cores = getOption("mc.cores", 1),
  calibratePheno = NULL
)
}
\arguments{
\item{net}{Object returned from \code{\link{asvNet}} with clusters added using \code{\link{netClust}}}

\item{asvTab}{An asv table with phenotypes joined.}

\item{asvCols}{A character vector of columns representing microbes (nodes).
Defaults NULL where all column names containing the string "ASV" will be used.}

\item{clusterCol}{The column name to use for clusters. If NULL then the first column name containing
the string "cluster" is used.}

\item{cluster}{A vector of clusters to be regressed against. By default this is NULL and all
clusters will be used.}

\item{phenoCols}{A vector of phenotype columns to use in changepoint regression.
These should be column names in the ASV table}

\item{model}{Type of changepoint model in chngpt::chngptm labeling convention. Currently hinge,
upperhinge, and segmented are supported. See Figure 2.1 of the chngpt
\link{[vignette](https://cran.r-project.org/web/packages/chngpt/vignettes/chngpt-vignette.pdf)}}

\item{cores}{Number of cores to run in parallel, defaults to 1 if "mc.cores" option is not set.}

\item{calibratePheno}{An optional vector of column names to calibrate the phenotypes by.
This should generally correspond to those used in `cal` if the ASV table has been calibrated or
just represent confounders that you wish to remove effects from in the changepoint regression.
\code{lm} is used to get residuals of the phenotype after these effects are regressed out.}
}
\value{
A named list (same as net) with three elements:
\itemize{
   \item{"Nodes" is a dataframe of nodes and their metadata}
   \item{"Edges" is a dataframe of edges connecting nodes.}
   \item{"graph" is the igraph object used to generate the dataframes.}
}
}
\description{
Function to aggregate clustered networks overall phenotypic impacts.
}
\examples{

# a<-qc(); b<-cal(a); c<-thresh(b); d<-asvDist(a) ; e<-net(d, thresh = c)
print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/net_output_clustered.rdata"))
print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/cal_output.rdata"))
asv <- are_c[[1]]
zinbCalibrated <- are_c[[2]][are_c[[2]]$model == "ZINB", "asv"]
asv <- are_c[[1]][, c("tissue", "plot", "row", "genotype", "biomass", "sd", zinbCalibrated)]

netThresh_output <- netThresh(net_data,
  asvTab = asv, asvCols = NULL,
  clusterCol = "kmeans_cluster", cluster = NULL, phenoCols = "biomass",
  model = "hinge", calibratePheno = "genotype"
)

}
\keyword{changepoint}
\keyword{network,}
