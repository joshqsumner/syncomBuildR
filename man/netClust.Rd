% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netClust.R
\name{netClust}
\alias{netClust}
\title{Function to cluster networks generated from \code{asvNet}.}
\usage{
netClust(net, method = "components", ...)
}
\arguments{
\item{net}{Object returned from \code{\link{asvNet}}}

\item{method}{Method to use for clustering. This can be a method from "components", "dbscan",
and "kmeans" or the output from pullNode, in which case that node and it's connections are labelled
as a cluster.}

\item{...}{Additional arguments passed to function picked by method.}
}
\value{
A named list (same as net) with three elements, same as \code{\link{net}}:
\itemize{
   \item{"Nodes" is a dataframe of nodes and their metadata}
   \item{"Edges" is a dataframe of edges connecting nodes.}
   \item{"graph" is the igraph object used to generate the dataframes.}
}
}
\description{
Function to cluster networks generated from \code{asvNet}.
}
\details{
Each method will use a different function to cluster data according to the layout
         in the graph. Note that layouts in the graph are determined by \code{igraph::layout_nicely}.
\itemize{
   \item{"component" uses \code{igraph::components} to cluster data and
         requires no additional arguments.}
   \item{"dbscan" uses \code{dbscan::dbscan} to cluster data.
         This requires at least that the eps argument is set. See \code{?dbscan::dbscan}.}
   \item{"kmeans" uses \code{stats::kmeans} to cluster data.
         This requires at least that the centers argument is set.}
}
}
\examples{

# a<-qc(); b<-cal(a); c<-thresh(b); d<-asvDist(a) ; e<-net(d, thresh = c)
print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/net_output.rdata"))
table(netClust(net = net_data, "component")[["nodes"]]$component_cluster)
table(netClust(net = net_data, "dbscan", eps = 3)[["nodes"]]$dbscan_cluster)
table(netClust(net = net_data, "kmeans", centers = 3)[["nodes"]]$kmeans_cluster)

net_data <- netClust(net = net_data, "component")
net_data <- netClust(net = net_data, "dbscan", eps = 3)
net_data <- netClust(net = net_data, "kmeans", centers = 3)

}
\keyword{changepoint}
\keyword{network,}
