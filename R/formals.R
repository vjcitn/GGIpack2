# any ggiResource should explicitly have seqnames, score,
# start, end, space (Genome tag, like hg19)

#' manage a GGI resource which must have start, end, seqnames, space, score
#' it also can have a tbl element
#' @export
setClass("ggiResource", slots=c(start="numeric",
 end="numeric", seqnames="ANY", space="ANY",
 score="numeric", tbl="ANY"))

setMethod("show", "ggiResource", function(object) {
 callNextMethod()
})

#' ABRIGresource is tailored to ABRIG data resources
#' @export
setClass("ABRIGresource", contains="ggiResource")

setMethod("show", "ABRIGresource", function(object) {
 cat("CDNM ggiResource\n")
 print(slot(object, "tbl"))
})


#' GTEx resource is tallored to the  wholeblpl05.parquet and lungpl05.parquet  resources.
#' @export
setClass("GTExresource", contains="ggiResource")

#' GTEx resource is tallored to the  wholeblpl05.parquet and lungpl05.parquet  resources.
#' @export
setClass("GTExresource", contains="ggiResource")

#' constructor for GTEx examples
#' @param con DBI connection
#' @param space character(1) must indicate build
#' @param tisstag character(1) string for tissue name
#' @param pfile character(1) path to parquet file
#' @note 'score' is p-value
#' @examples
#' fo = ggi_gtex_cache("lungpl05.parquet")
#' if (nchar(fo)>0) {
#'   lungpa = fo
#'   con = DBI::dbConnect(duckdb::duckdb())
#'   lu = GTExresource(con, tisstag="lung", pfile=lungpa)
#'   lu
#' # please do dbDisconnect(con) when done with 'lu'
#' }
#' @export
GTExresource = function (con, space = "hg19", tisstag, pfile) 
{
    pp = sprintf("read_parquet(%s)", sQuote(pfile, q=FALSE))
    tb = tbl(con, pp)
    ans = dplyr::mutate(tb, score = pvalue, seqnames = chromosome, tissue=tisstag,
       start=position, end=position) |> select(tissue, variant, rsid, molecular_trait_id, score,
            maf, beta, se, seqnames, ref, alt, start, end)
    new("GTExresource", space = space, tbl = ans)
}

