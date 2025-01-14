#!/usr/bin/env -S Rscript --vanilla

.libPaths()
# installed.packages()

options(HTTPUserAgent = sprintf("R/%s R (%s)", rver, paste(rver, R.version$platform, R.version$arch, R.version$os)))

rver <- getRversion()
distro <- system2('lsb_release', '-sc', stdout = TRUE)
str(rver)
str(distro)

#pak::repo_add(binaries = "https://bioc.r-universe.dev/bin/linux/focal/4.4")
#pak::repo_add(r_universe = "https://bioc.r-universe.dev/")

# pak::repo_status()

install.packages('Rhtslib', repos = c('https://bioc.r-universe.dev'), type = "binary")


#str(BiocManager::version())
#str(BiocManager::containerRepository(version = BiocManager::version(), type = "binary"))
#str(BiocManager:::BINARY_BASE_URL)

#binaries = "https://bioc.r-universe.dev/bin/linux/noble/4.4",
#source = "https://r-lib.r-universe.dev"

#message("checking")
#checkResults <- devtools::check(error_on = c("never"), env_vars = c(`_R_CHECK_TESTS_NLINES_` = '0', `CI` = 'true'))
message("pkg_deps_tree")
pak::pkg_deps_tree("local::.", upgrade = FALSE, dependencies = TRUE)

message("pkg_install")
pak::pkg_install("local::.", upgrade = FALSE, dependencies = TRUE)
