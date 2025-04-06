#' convert the tabular information from a ggiResource instance into GRanges
#' @import GenomicRanges
#' @param tbl resource with seqnames, start, end
#' @param genome character(1) should come from resource 'space' element
#' @param maxn.msg numeric(1) if number of records for conversion to GRanges exceeds this, a message is given
#' @export 
toGRanges = function(tbl, genome=NA, maxn.msg=5e4) {
  nrec = tbl |> dplyr::count() |> as.data.frame() |> unlist() |> as.numeric()
  if (nrec > maxn.msg) message(sprintf("you are asking for %d records but maxn is %d\n", nrec, maxn))
  df = as.data.frame(tbl) # expensive?
  g = GenomicRanges::GRanges(seqnames=df$seqnames, IRanges::IRanges(start=df$start,
    end=df$end))
  mcols(g) = df |> dplyr::select(!c(seqnames, start, end))
  genome(g) = genome
  g
}
 
