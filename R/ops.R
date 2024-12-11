#' define constructor for ABRIGresource 
#' @import dplyr
#' @param con is a DBI connection (typically duckdb)
#' @param tissue character(1)
#' @param space character(1) e.g., "hg19"
#' @param pfiles list of absolute paths to the data for each tissue. 
#' @examples
#' con = DBI::dbConnect(duckdb::duckdb())
#' ll = ABRIGresource( con, "BAL" , pfiles= ABRIGparquet_paths())
#' print(ll)
#' DBI::dbDisconnect(con)
#' @return ABRIGresource instance for the identified tissue
#' @export
ABRIGresource = function(con, tissue, space="hg19", pfiles) {
   ttypes = c("BAL", "BroncEpiBrush", "CD4Stim", "CD4Unstim",
                  "AlvMacphage", "PaxRNA")
   stopifnot(tissue %in% ttypes)
   ans = dplyr::tbl(con, sprintf("read_parquet(%s)", sQuote(pfiles[tissue])), q=FALSE) |> dplyr::mutate(score=FDR, seqnames=CHR, start=BP, molecular_trait_id=gene, trait2=geneId)
   new("ABRIGresource", space=space, tbl=ans)
}


#' filter a GGI (ABRIGresource) instance by range
#' @rawNamespace import(GenomicRanges, except=c(intersect, union, setdiff))
#' @param res ABRIGresource instance
#' @param ggr GenomicRanges instance
#' @param tag character(1) value in `ggr_field` used for filtering, e.g., a gene
#' @param radius numeric(1) flanking region size
#' @param ggr_field character(1) metadata element in mcols(ggr) for filtering
#' @examples
#' con = DBI::dbConnect(duckdb::duckdb())
#' ll = ABRIGresource( con, "BAL" , pfiles= ABRIGparquet_paths())
#' utils::data("gloc_hg19", package = "GGIpack2")
#' kk <- filterByRange(ll, gloc_hg19, "DSP", ggr_field="gene_name")
#' print(kk)
#' DBI::dbDisconnect(con)
#' @return filterByRange instance for the identifies gene_name or gene_id 
#' @export
filterByRange = function(res, ggr, tag, radius=1e5, ggr_field="gene_name") {
  stopifnot(inherits(res, "ABRIGresource"))
  ok = ggr[which(GenomicRanges::mcols(ggr)[[ggr_field]] == tag)]
  stopifnot(width(ok)>0)
  anac = function(x) as.numeric(as.character(x)) # for Rle
  ans = methods::slot(res, "tbl") |> dplyr::filter(CHR == local(anac(GenomicRanges::seqnames(ok)[1])), BP >= local(IRanges::start(ok)-radius), BP<= local(IRanges::end(ok)+radius))
  new("ABRIGresource", space = methods::slot(res, "space"), tbl=ans)
}

#' Takes in a path to a data table and coheres the data into the proper format 
#' @param path The complete path to the DN8 file.
#' @examples
#' path = system.file("extdata/Alveolar_Macrophages_IS.MICA_ILMN_3241692.CAU.meta", package="GGIpack2")
#' head(checkData(path = path))
#' @return a data set that can be used to graph in JBrowseR.
#' @export
checkData = function(path){
    df = utils::read.table(path, header = TRUE)
    needCols = c("BP", "P", "CHR")
    stopifnot(needCols %in% colnames(df))
    data =dplyr::mutate(df, start = BP, end =BP)
    return(data)
}


#' make GWASTrack for shiny app. 
#' @importFrom igvShiny GWASTrack
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
  ndat = names(dat)
  pindex = which(ndat == "P")
  bpindex = which(ndat == "BP")
  chrindex = which(ndat == "CHR")
  igvShiny::GWASTrack(trackName = name, data = dat, chrom.col=chrindex, pos.col = bpindex, 
                      pval.col = pindex)
}

#' this makes a data.frame
#' @param con is a DBI connection (typically duckdb)
#' @param tissue character(1)
#' @param gene character(1)
#' @return a data frame filters by tissue and gene.
#' @examples
#' con = DBI::dbConnect(duckdb::duckdb())
#' nn = make_data_frame_from_tissue_and_gene(con, "BAL", "DSP")
#' DBI::dbDisconnect(con)
#' @export
make_data_frame_from_tissue_and_gene = function(con, tissue, gene) {
  # add code to validate tissue and gene
  ll = ABRIGresource( con, tissue , pfiles= ABRIGparquet_paths())
  utils::data("gloc_hg19", package = "GGIpack2")
  kk <- filterByRange(ll, gloc_hg19, gene, ggr_field="gene_name")
  tmp = as.data.frame(kk@tbl)
  stopifnot(inherits(tmp, "data.frame"))
  tmp
}


#' Takes in a metaABRIG data table and rounds the SE MAF BETA FDR and statistic to 3 digits. For P (pvalue) that variable is turned into  the scientific notation then rounded to 3 digits. 
#' @import dplyr
#' @param mydf  metaABRIG data table
#' @return a metaABRIG data table with rounded values . 
#'@examples
#'con = DBI::dbConnect(duckdb::duckdb())
#' nn = make_data_frame_from_tissue_and_gene(con, "BAL", "DSP")
#' dorounds(mydf= nn )
#' DBI::dbDisconnect(con)
#' @export
dorounds = function(mydf) {
  mydf$P = formatC(mydf$P, format = "e", digits= 3)
  mydf$SE = round(mydf$SE, 3)
  mydf$MAF = round(mydf$MAF, 3)
  mydf$BETA = round(mydf$BETA, 3)
  mydf$FDR= round(mydf$FDR, 3)
  mydf$statistic= round(mydf$statistic, 3)
  mydf |> dplyr::select(-score, -seqnames, -SNP)
}



#' Get all the datasets from each of the six celltypes in the ABRIG data.
#' @param gene A gene  from the human genome.
#' @param con  A duckdb connection
#' @param genelocs  GenomicRanges instance
#' @param pfiles list of absolute paths to the data for each tissue.
#' @return A list of data table by type for the given gene.
#' @examples
#' utils::data("gloc_hg19", package = "GGIpack2")
#' con = DBI::dbConnect(duckdb::duckdb())
#' allrefs( con =con, gene = 'DSP', pfiles =ABRIGparquet_paths(), genelocs = gloc_hg19)
#' DBI::dbDisconnect(con)
#' @export
allrefs = function(con, pfiles, genelocs, gene){
  ttypes = c("BAL", "BroncEpiBrush", "CD4Stim", "CD4Unstim",
            "AlvMacphage", "PaxRNA")
  allres = lapply(ttypes, function(x) ABRIGresource(con, x, pfiles=pfiles))
  allfilt = lapply(allres, function(x) filterByRange(x,
                                                     genelocs, gene, ggr_field="gene_name"))
  names(allfilt) = ttypes
  return(allfilt)
}



