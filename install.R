#!/usr/bin/env -S Rscript --vanilla

# requirements:
#	ubuntu 24.04 "noble"
#		with lsb_release utility installed
#	r >= 4.4.1
#	r "pak" library -> https://github.com/r-lib/pak/

library(pak)

# suppresses missing pillar package message
options(pak.no_extra_messages = TRUE)

.libPaths()

# config pak repos based on host
source("https://changit.bwh.harvard.edu/raw/jenkins/standard-pipelines/rekgm-ci-testing/resources/R/CI/pak/config_pak_repos.R?token=GHSAT0AAAAAAAAAAL35ZK5HEAB6TLZBBZNSZ4ZDIAQ")

pak::repo_status()


message("pkg_deps_tree")
pak::pkg_deps_tree("local::.", upgrade = FALSE, dependencies = TRUE)

# could be a flag to control this
#checkResults <- devtools::check(error_on = c("never"), env_vars = c(`_R_CHECK_TESTS_NLINES_` = '0', `CI` = 'true'))

# use devtools build/install instead...?
# can i install the package from the source package made by devtools::check?
# does this run tests when building?
message("pkg_install")
install_results <- pak::pkg_install("local::.", upgrade = FALSE, dependencies = TRUE)
#str(install_results)

#library(jsonlite)

# flag for this
library(jsonlite)

cat(toJSON(install_results, pretty = TRUE, force = TRUE))




# only need to run this to make the docker image smaller
pak::pak_cleanup(force = TRUE)

installed.packages()
