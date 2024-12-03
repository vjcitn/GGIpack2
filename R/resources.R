#' provide paths to all ABRIG parquet
#' @export
ABRIGparquet_paths = function() {
pfiles =
c("/udd/remcr/abrig/ca_ba.parquet", "/udd/remcr/abrig/ca_brep.parquet",
"/udd/remcr/abrig/cc_s.parquet", "/udd/remcr/abrig/cc_unstim.parquet",
"/udd/remcr/abrig/ch_am.parquet", "/udd/remcr/abrig/comb.parquet")

names(pfiles) = c("BAL", "BroncEpiBrush", "CD4Stim", "CD4Unstim",
                  "AlvMacphage", "PaxRNA")
pfiles
}

