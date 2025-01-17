library(testthat)

con = DBI::dbConnect(duckdb::duckdb())

test_that("ABRIGparquet_paths", {
  test = ABRIGparquet_paths() 
  answer = c("/udd/remcr/abrig/ca_ba.parquet", "/udd/remcr/abrig/ca_brep.parquet",
             "/udd/remcr/abrig/cc_s.parquet", "/udd/remcr/abrig/cc_unstim.parquet",
             "/udd/remcr/abrig/ch_am.parquet", "/udd/remcr/abrig/comb.parquet")
  
  names(answer) = c("BAL", "BroncEpiBrush", "CD4Stim", "CD4Unstim",
                    "AlvMacphage", "PaxRNA")
  
  expect_identical(test, answer)
  
})
 

test_that("filterByRange", {
 ll = ABRIGresource( con, "BAL" , pfiles= ABRIGparquet_paths())
  utils::data("gloc_hg19", package = "GGIpack2")
  BAL_DSP <- filterByRange(ll, gloc_hg19, "DSP", ggr_field="gene_name")
  BAL_DSP <- BAL_DSP@tbl |> as.data.frame()
  answerPath <- system.file("extdata", "BAL_DSP.rds", package = "GGIpack2")
  answer <- readRDS(file= answerPath)
  expect_equal(BAL_DSP, answer)
 #expect_match(BAL_DSP, answer)
})


test_that("ABRIGresource",
          {BAL_ABRIGresource = ABRIGresource( con, "BAL" , pfiles= ABRIGparquet_paths())
          answerPath <- system.file("extdata", "BAL_ABRIGresource.rds", package = "GGIpack2")
          answer <- readRDS(file= answerPath)
          expect_equal(BAL_ABRIGresource, answer)
          })


DBI::dbDisconnect(con, shutdown=TRUE)


test_that("Symbol mapping is OK", {
 data("gloc_hg19", package="GGIpack2")
 data("newensg", package="GGIpack2")
 ensg = newensg[order(names(newensg))] # preserve old name

 gsyms = names(ensg)
 names(gsyms) = as.character(ensg)
 e2g = function(e) {
   ind = match(e, names(gsyms))
   dr = which(is.na(ind))
   ans = e
   if (length(dr)>0)
   ans[-dr] = gsyms[ind[-dr]]
   else ans = gsyms[ind]
   ans
 }
 expect_equal(e2g("ENSG00000073605"), "GSDMB")
})
