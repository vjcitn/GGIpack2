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
gtexapp()