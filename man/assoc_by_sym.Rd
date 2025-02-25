% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assoc_by_sym.R
\name{assoc_by_sym}
\alias{assoc_by_sym}
\title{produce data.frame for GTExresource content filtered for a specific gene}
\usage{
assoc_by_sym(gtexres, sym = "ORMDL3")
}
\arguments{
\item{gtexres}{GTExresource instance}

\item{sym}{character(1) gene symbol}
}
\description{
produce data.frame for GTExresource content filtered for a specific gene
}
\examples{
lu = ggi_gtex_cache("lungpl05.parquet")
con = DBI::dbConnect(duckdb::duckdb())
lures = GTExresource(con, tisstag="lung", pfile=lu)
aa = assoc_by_sym( lures )
DBI::dbDisconnect(con)
}
