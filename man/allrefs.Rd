% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ops.R
\name{allrefs}
\alias{allrefs}
\title{Get all the datasets from each of the six celltypes in the ABRIG data.}
\usage{
allrefs(con, pfiles, genelocs, gene)
}
\arguments{
\item{con}{A duckdb connection}

\item{pfiles}{list of absolute paths to the data for each tissue.}

\item{genelocs}{GenomicRanges instance}

\item{gene}{A gene  from the human genome.}
}
\value{
A list of data table by type for the given gene.
}
\description{
Get all the datasets from each of the six celltypes in the ABRIG data.
}
\examples{
utils::data("gloc_hg19", package = "GGIpack2")
con = DBI::dbConnect(duckdb::duckdb())
allrefs( con =con, gene = 'DSP', pfiles =ABRIGparquet_paths(), genelocs = gloc_hg19)
DBI::dbDisconnect(con)
}
