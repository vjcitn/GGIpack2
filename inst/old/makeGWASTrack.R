
#' make GWASTrack for shiny app. 
#' @param dat the data table exported from the find_data 
#' @param name the name of the track in the shiny app
#' @return a gwasTrack read to graph by igvshiny
#' @note
#' It is required that the column  header need to be the following:
#' c("SNP", "CHR", "BP", "A1", "A2", "gene", "geneId", "statistic", 
#'  "P", "FDR", "BETA", "SE", "MAF")
#'  It is  advised if the desired data does not have these column headers that the column headers need to be changed to the above. 
#' @examples
#'con = DBI::dbConnect(duckdb::duckdb())
#' nn = make_data_frame_from_tissue_and_gene(con, "BAL", "DSP")
#' makeGWASTrack(dat = nn)
#' @export
makeGWASTrack = function( name="NA", dat) {
  if (!requireNamespace("igvShiny")) stop("install igvShiny to use this function")
  ndat = names(dat)
  pindex = which(ndat == "P")
  bpindex = which(ndat == "BP")
  chrindex = which(ndat == "CHR")
  igvShiny::GWASTrack(trackName = name, data = dat, chrom.col=chrindex, pos.col = bpindex, 
                      pval.col = pindex)
}

