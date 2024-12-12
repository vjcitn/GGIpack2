#' A vector of gene names.
#' @docType data
#' @format vector
"geneNames"



#' A instance of GenomicRanges human genome 19
#' @docType data
#' @format data object
"gloc_hg19"

#' A map from symbols to Ensembl IDs
#' @docType data
#' @note this was made primarily to allow
#' mapping of GTEx v7 molecular_trait_ids to 
#' symbols.  
#' @format data object
"ensg"

#' map from symbols to v79 Ensembl IDs
#' @docType data
#' @note previous 'ensg' vector had lots of missing
#' symbol mappings.  This combines a more recent
#' ensembl catalog with 'obsolete' mappings that
#' were present in ensg.
#' @format named character vector
"newensg"
