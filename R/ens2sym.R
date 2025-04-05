#' map naked ensembl gene ids to symbols when available in gmap,
#' otherwise just pass back the ENSG tag
#' @param e character() ENSG symbols
#' @param gmap a character vector with gene symbols as element values and ENSG tags as element names
#' @examples
#' ensg = setup_ens_gnames()
#' gsyms = names(ensg)
#' names(gsyms) = as.character(ensg)
#' head(ens2sym(head(ensg), gsyms))
#' @export
 ens2sym = function(e, gmap) {
   ind = match(e, names(gmap))
   dr = which(is.na(ind))
   ans = e
   if (length(dr)>0)
      ans[-dr] = gmap[ind[-dr]]
   else ans = gmap[ind]
   ans
 }
