#!/usr/bin/env -S Rscript --vanilla

# requirements:
#	ubuntu 24.04 "noble"
#		with lsb_release utility installed
#	r >= 4.4.1
#	r "pak" library -> https://github.com/r-lib/pak/

library(pak)

str(tempdir())

# suppresses missing pillar package message
options(pak.no_extra_messages = TRUE)

.libPaths()



rver <- getRversion()
distro <- system2('lsb_release', '-sc', stdout = TRUE)

# so if the ubuntu version is "noble", then rig won't set up the PPPM repo because "unsupported", which is bogus
# and if the ubuntu distro is *not* "noble", then pak can't find the bioconductor binaries at r-universe, 
# because r-universe only has the binaries for "noble" and R 4.4 + 4.5
# so we add this repo back in, even though rig really should have put it in already
# ONLY WORKS FOR AMD64
p3m_binaries <- sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", distro)

# seems to understand arm64 and amd64, but only R >= 4.4
r_universe_bioc_binaries <- sprintf('%s/bin/linux/%s/%s', "https://bioc.r-universe.dev", distro, substr(rver, 1, 3))

# save this for later
#rspm_binaries <- sprintf("https://packagemanager.rstudio.com/cran/__linux__/%s/latest", distro)

# NEED THIS FOR ARM (super duper hardcoded, unofficial, may go away without warning)
# https://forum.posit.co/t/arm64-binary-packages-for-linux-in-posit-public-package-manager/178514/4
# https://github.com/r-hub/repos?tab=readme-ov-file#ubuntu-2404--r-release-on-aarch64-ubuntu-2404-aarch64-r44
github_binaries <- sprintf("https://raw.githubusercontent.com/r-hub/repos/main/ubuntu-24.04-aarch64/%s", substr(rver, 1, 3))

pak::repo_add(r_universe_bioc_binaries = r_universe_bioc_binaries)

machine <- Sys.info()["machine"]
# "aarch64" or "x86_64"

if (identical(unname(machine), "aarch64"))
{
	message("adding repo for aarch64 (arm64) binaries")
	pak::repo_add(github_binaries = github_binaries)
} else if (identical(unname(machine), "x86_64"))
{
	message("adding repo for x86_64 binaries")
	pak::repo_add(p3m_binaries = p3m_binaries)
} else 
{ 
	message( paste("unrecognized binary architecture: ", machine, ", enjoy compiling from source"))
}

#pak::repo_add(rspm_binaries = rspm_binaries)
pak::repo_status()




#message("checking")
#checkResults <- devtools::check(error_on = c("never"), env_vars = c(`_R_CHECK_TESTS_NLINES_` = '0', `CI` = 'true'))
message("pkg_deps_tree")

#pak::pkg_deps_tree("local::.", upgrade = FALSE, dependencies = TRUE)
pak::pkg_deps_tree("cli", upgrade = FALSE, dependencies = TRUE)

# use devtools build/install instead...?
message("pkg_install")
#pak::pkg_install("local::.", upgrade = FALSE, dependencies = TRUE)
pak::pkg_install("cli", upgrade = FALSE, dependencies = TRUE)

# only need to run this to make the docker image smaller
pak::pak_cleanup(force = TRUE)

installed.packages()

# clean up tmpdir, for some reason this needs to be run after installed.packages() is run
#unlink(tempdir(), recursive = T)
