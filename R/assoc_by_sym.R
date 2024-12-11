#' produce data.frame for GTExresource content filtered for a specific gene
#' @param gtexres GTExresource instance
#' @param sym character(1) gene symbol
#' @examples
#' lu = ggi_gtex_cache("lungpl05.parquet")
#' con = DBI::dbConnect(duckdb::duckdb())
#' lures = GTExresource(con, tisstag="lung", pfile=lu)
#' aa = assoc_by_sym( lures )
#' DBI::dbDisconnect(con)
#' @export
assoc_by_sym = function( gtexres, sym = "ORMDL3" ) {
 # assumes number of associations for a gene is manageable as data.frame
 stopifnot(is(gtexres, "GTExresource"))
 data("gloc_hg19", package="GGIpack2")
 mid = which(gloc_hg19$symbol == sym)
 stopifnot(length(mid)>=1)
 mid = names(gloc_hg19[mid[1]]) # should message if length(mid)>1
 tmp = gtexres@tbl |> dplyr::filter(molecular_trait_id == mid) |> as.data.frame()
 dplyr::mutate(tmp, genesym=sym)
}


#' produce data.frame for GTExresource content filtered for a specific position/radius
#' @param gtexres GTExresource instance
#' @param chr character(1) seqnames value
#' @param pos numeric(1) central coordinate of view
#' @param radius numeric(1) flanking view region
#' @examples
#' lu = ggi_gtex_cache("lungpl05.parquet")
#' con = DBI::dbConnect(duckdb::duckdb())
#' lures = GTExresource(con, tisstag="lung", pfile=lu)
#' bb = assoc_by_pos( lures, chr="17", pos=39e6, radius=2e5 )
#' ggplot2::ggplot(bb, ggplot2::aes(x=start, y=-log10(score), 
#'    colour=molecular_trait_id)) + ggplot2::geom_point()
#' if (interactive()) {
#'  dd = ggplot2::ggplot(dplyr::mutate(bb, 
#'   rsidg=paste(rsid,e2s(molecular_trait_id),sep=":")), 
#'     ggplot2::aes(x=start, y=-log10(score), 
#'     colour=molecular_trait_id, text=rsidg)) + ggplot2::geom_point()
#'    plotly::ggplotly(dd)
#'  }
#' DBI::dbDisconnect(con)
#' @export
assoc_by_pos = function( gtexres, chr="17", pos=39e6, radius=1e6) {
 # assumes number of associations for a gene is manageable as data.frame
 stopifnot(is(gtexres, "GTExresource"))
 data("gloc_hg19", package="GGIpack2")
 tmp = gtexres@tbl |> dplyr::filter(seqnames == chr & start >= (pos-radius) &
          end <= (pos+radius)) |> as.data.frame()
 dplyr::mutate(tmp, pos=pos, radius=radius)
}


#' transform ENSG ids to symbols where feasible
#' @param e character() vector of ENSG
#' @export
e2s = function(e) {
  ok = which(e %in% names(gloc_hg19)) 
  ans = e
  ans[ok] = gloc_hg19[e[ok]]$symbol
  ans
}

