
test_that("new file mgt works",{
 library(GGIpack)
 lu = ggi_gtex_cache("lungpl05.parquet")
 expect_true(nchar(lu)>0)

 con = DBI::dbConnect(duckdb::duckdb())
 lures = GTExresource(con, tisstag="lung", pfile=lu)
 expect_true( (lures@tbl |> dplyr::count() |> dplyr::collect()) == 16959106L )
 orm = lures@tbl |> dplyr::filter(molecular_trait_id == "ENSG00000172057")
 expect_true( (orm |> dplyr::count() |> dplyr::collect()) == 1555L )
 DBI::dbDisconnect(con)
} )
