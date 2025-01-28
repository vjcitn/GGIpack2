#!/usr/bin/env -S Rscript --vanilla

library(GGIpack2)

ggi_gtex_cache("lungpl05.parquet")  # may ask for creation of cache folder
ggi_gtex_cache("wholeblpl05.parquet")
options(useFancyQuotes=FALSE)
options(bitmapType="cairo")
# needed for docker
options(shiny.port = 8090)
# needed for docker
options(shiny.host = "0.0.0.0")

# Sys.getenv returns a DList...?
if (is.na(Sys.getenv()["GGIPACK2_SHINY_MODE"]))
{
	message("GGIPACK2_SHINY_MODE env var not found")
} else if (identical(unname(unlist(as.list(Sys.getenv()["GGIPACK2_SHINY_MODE"]))), "gtexapp"))
{
	message("starting gtexapp")
	gtexapp()
} else if (identical(unname(unlist(as.list(Sys.getenv()["GGIPACK2_SHINY_MODE"]))), "abrigapp2"))
{
	message("starting abrigapp2")
	abrigapp2()
} else
{
	message(paste("unrecognized GGIPACK2_SHINY_MODE value: '", unname(Sys.getenv()["GGIPACK2_SHINY_MODE"]), "', not starting the app"))
}