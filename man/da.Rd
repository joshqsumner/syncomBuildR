% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/da.R
\name{da}
\alias{da}
\title{Function to assess differential abundance of microbe counts or pseudo counts}
\usage{
da(df, col, predictors, zi_cutoff = 0.1)
}
\arguments{
\item{df}{ASV table with metadata and counts}

\item{col}{ASV column to model}

\item{predictors}{Column names of predictors, can be discrete or continuous.}

\item{zi_cutoff}{Proportion of zeros to consider the distribution to be zero inflated,
defaults to 0.1.}
}
\value{
A data frame of model results.
}
\description{
Function to assess differential abundance of microbe counts or pseudo counts
}
\examples{

# Using uncalibrated data

print(load("~/Desktop/stargate/SINC/sincUtils/syncomBuilder/qc_output.rdata"))
asv[1:10, 1:10]

df <- asv[, 1:20]
sapply(df[, 3:20], function(c) mean(c == 0))

col <- "ASV9"
predictors <- "tissue"
zi_cutoff <- 0.1

df[[predictors]] <- relevel(factor(df[[predictors]]), "AS")

da(df, col, predictors, zi_cutoff)
}
