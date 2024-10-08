% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demultiplex.R
\name{demultiplex}
\alias{demultiplex}
\title{Function to demultiplex sequences either to single files or to multiple files for use with dada2.}
\usage{
demultiplex(
  reads,
  barcodes,
  fwd = "FWD",
  rev = "REV",
  name = c("Name", "Sample_name"),
  mode = c("dada", "general"),
  cores = getOption("mc.cores", 1),
  writeOut = TRUE,
  stat = FALSE
)
}
\arguments{
\item{reads}{A ShortReadQ object or a file path to a fastq file that will be read in with
\code{ShortRead::readFastq}.}

\item{barcodes}{A dataframe of barcodes and metadata for naming or a file path to such a flat file
readable by \code{utils::read.delim}.}

\item{fwd}{The column name of the barcodes dataframe representing forward barcodes.
Defaults to "FWD".}

\item{rev}{The column name of the barcodes dataframe representing reverse barcodes.
Defaults to "REV".}

\item{name}{One or more column names from the barcodes dataframe to be used in naming samples.}

\item{mode}{One of "dada" or "general". This controls how files are written out.
If mode is "dada" then a fastq file will be written out for
each row of the barcodes dataframe and the files will be named using the name argument.
If mode is "general" then each input fastq file will yield one
output fastq file with sample names in that file corresponding to the name argument.
Defaults to "dada".}

\item{cores}{Optionally the number of cores to run in parallel. This defaults to 1 if the "mc.cores"
option is not set.}

\item{writeOut}{Logical, should fastq/stats files be written out. Defaults to TRUE.
This can also be a path where data should be written to. Note that files are written as mode "w"
and this will throw an error if the file already exists.}

\item{stat}{Logical, should stats files be generated if writeOut is TRUE. Defaults to FALSE.}
}
\value{
Writes out fq files and summary statistics depending on writeOut and stat arguments.
If both are TRUE then summary stats are returned as a dataframe.
}
\description{
Function to demultiplex sequences either to single files or to multiple files for use with dada2.
}
\examples{
# these examples are not run because they require specific fastq files.
\dontrun{
library(ShortRead)
library(Biostrings)
library(parallel)
setwd("~/scripts/SINC/sincUtils/syncomBuilder/nastya96WellEx")

barcodes <- read.delim("barcode_tab.tsv")
name <- c("Name", "Well")
reads <- ShortRead::readFastq("Plate11_R1_001.fastq.gz")

demultiplex(reads, barcodes[1:10, ],
  fwd = "FWD", rev = "REV",
  name = name, mode = "dada", cores = 10, writeOut = "tests", stat = FALSE
)

path <- paste0(
  "/run/user/1000/gvfs/smb-share:server=nas01,share=research/bart_lab/",
  "Previous_people/Mingsheng_Qi/research_dr/SINC/limiting_dilution/20210314/p2p11_rawseq"
)
barcodes <- read.csv(paste0(path, "/barcodetable1.csv"))
barcodes$exp <- "exp"
reads <- ShortRead::readFastq(paste0(path, "/210314S05P10_merge.fq"))

x <- demultiplex(reads, barcodes,
  fwd = "FWD_primer", rev = "REV_primer", name = c("exp"), mode = "general",
  writeOut = TRUE, stat = TRUE, cores = 10
)
x
# benchmarking against:
# "plate","totalReads","groupedReads","ratio"
# "P10",234509,154507,0.6589
}

}
\keyword{dada2}
\keyword{demultiplex,}
