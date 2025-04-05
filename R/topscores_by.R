
#' get scores of selected types from a GGIresource instance
#' @param res GGIresource instance
#' @param scope character(1) one of c("chr", "gene", "rsid")
#' @param val a value in the selected `scope`
#' @param n numeric(1) used to filter results through `head()`
#' @examples
#' # we rely on duckdb for query resolution
#' blex = system.file("parquet", "testbl.parquet", package="GGIpack2")
#' pp = sprintf("read_parquet(%s)", sQuote(blex, q=FALSE))
#' con = DBI::dbConnect(duckdb::duckdb())
#' tb = dplyr::tbl(con, pp)
#' nres = new("GTExresource", space = "hg19", tbl = tb)
#' topscores_by(nres, scope="chr", val="1", n=5)
#' topscores_by(nres, scope="gene", val="FO538757.2", n=5)
#' topscores_by(nres, scope="rsid", val="rs368811019", n=5)
#' DBI::dbDisconnect(con)
#' @export
topscores_by = function(res, scope="chr", val, n=100) {
   if (scope=="chr")
          dat = slot(res, "tbl") |> dplyr::filter(seqnames == as.character(local(val))) |>
                  dplyr::arrange(score) |>
                  head(n) |> as.data.frame()
   else if (scope == "gene") {
          g2e = setup_ens_gnames()
          print(val)
          stopifnot(val %in% names(g2e))
          dat = slot(res, "tbl") |> dplyr::filter(molecular_trait_id == as.character(local(g2e[val]))) |>
                  dplyr::arrange(score) |>
                  head(n) |> as.data.frame()
          }
   else if (scope == "rsid") {
          dat = slot(res, "tbl") |> dplyr::filter(rsid == as.character(local(val))) |>
                  dplyr::arrange(score) |>
                  head(n) |> as.data.frame()
          }
   else stop("unrecognized value for scope")
   dat
}


