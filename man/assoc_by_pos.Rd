% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assoc_by_sym.R
\name{assoc_by_pos}
\alias{assoc_by_pos}
\title{produce data.frame for GTExresource content filtered for a specific position/radius}
\usage{
assoc_by_pos(gtexres, chr = "17", pos = 3.9e+07, radius = 1e+06)
}
\arguments{
\item{gtexres}{GTExresource instance}

\item{chr}{character(1) seqnames value}

\item{pos}{numeric(1) central coordinate of view}

\item{radius}{numeric(1) flanking view region}
}
\description{
produce data.frame for GTExresource content filtered for a specific position/radius
}
\examples{
lu = ggi_gtex_cache("lungpl05.parquet")
con = DBI::dbConnect(duckdb::duckdb())
lures = GTExresource(con, tisstag="lung", pfile=lu)
bb = assoc_by_pos( lures, chr="17", pos=39e6, radius=2e5 )
ggplot2::ggplot(bb, ggplot2::aes(x=start, y=-log10(score), 
   colour=molecular_trait_id)) + ggplot2::geom_point()
if (interactive()) {
 dd = ggplot2::ggplot(dplyr::mutate(bb, 
  rsidg=paste(rsid,e2s(molecular_trait_id),sep=":")), 
    ggplot2::aes(x=start, y=-log10(score), 
    colour=molecular_trait_id, text=rsidg)) + ggplot2::geom_point()
   plotly::ggplotly(dd)
 }
DBI::dbDisconnect(con)
}
