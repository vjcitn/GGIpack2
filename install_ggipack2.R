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

rver <- getRversion()
distro <- system2('lsb_release', '-sc', stdout = TRUE)

# so if the ubuntu version is "noble", then rig won't set up the PPPM repo because "unsupported", which is bogus
# and if the ubuntu distro is *not* "noble", then pak can't find the bioconductor binaries at r-universe, 
# because r-universe only has the binaries for "noble"
# so we add this repo back in, even though rig really should have put it in already
linux_binaries <- sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", distro)

bioc_binaries <- sprintf('%s/bin/linux/%s/%s', "https://bioc.r-universe.dev", distro, substr(rver, 1, 3))

# put this back in...?
#pak::repo_add(r_universe_bioc = "https://bioc.r-universe.dev/") # cran...?
pak::repo_add(r_universe_bioc_binaries = bioc_binaries)
pak::repo_add(p3m_retrofit = linux_binaries)
pak::repo_status()


#message("checking")
#checkResults <- devtools::check(error_on = c("never"), env_vars = c(`_R_CHECK_TESTS_NLINES_` = '0', `CI` = 'true'))
message("pkg_deps_tree")
pak::pkg_deps_tree("local::.", upgrade = FALSE, dependencies = TRUE)

# use devtools build/install instead...?
message("pkg_install")
pak::pkg_install("local::.", upgrade = FALSE, dependencies = TRUE)

installed.packages()