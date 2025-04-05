#' use newensg data element to provide a mapping between symbol and ENSG id
#' @examples
#' head(setup_ens_gnames())
#' @export
setup_ens_gnames = function() {
  data("newensg", package="GGIpack2")
  newensg[order(names(newensg))] # preserve old name
}

