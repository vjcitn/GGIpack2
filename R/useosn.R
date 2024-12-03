#wget https://mghp.osn.xsede.org/bir190004-bucket01/BiocGGIData/

#' add or retrieve selected parquet file from BiocFileCache
#' @param pfname character(1) name of a parquet file assumed present in Bioc Open Storage Network
#' @param cache BiocFileCache-type cache
#' @examples
#' ggi_gtex_cache()
#' @export
ggi_gtex_cache = function(pfname = "lungpl05.parquet", cache=BiocFileCache::BiocFileCache()) {
   url = paste0("https://mghp.osn.xsede.org/bir190004-bucket01/BiocGGIData/", pfname)
   inf = BiocFileCache::bfcquery(cache, url)
   if (nrow(inf)>0) {
     res = inf[nrow(inf),]
     message(sprintf("resource %s already in cache from %s\n", res$rid, url))
     return(res$rpath)
     } 
   if (!url_ok(url)) stop(sprintf("HEAD for %s does not return status code 200\n", url))
   BiocFileCache::bfcadd(cache, rname = basename(url), fpath=url, rtype="web")
}


#' check that a URL can get a 200 for a HEAD request
#' @importFrom httr status_code HEAD
#' @param url character(1)
#' @return logical(1)
url_ok = function(url) {
  httr::status_code(httr::HEAD(url)) == 200
}

