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
rver_short <- substr(rver, 1, 3)
machine <- Sys.info()["machine"]
# "aarch64" or "x86_64"
codename <- system2('lsb_release', '-sc', stdout = TRUE)
release <- system2('lsb_release', '-sr', stdout = TRUE)
lsb_id <- system2('lsb_release', '-si', stdout = TRUE)
lsb_desc <- system2('lsb_release', '-sd', stdout = TRUE)

message("rver: ", rver)
message("machine: ", machine)
message("codename: ", codename)
message("release: ", release)
message("lsb_id: ", lsb_id)
message("lsb_desc: ", lsb_desc)

# save this for later
#rspm_binaries <- sprintf("https://packagemanager.rstudio.com/cran/__linux__/%s/latest", codename)

# the goldilocks config here is ubuntu noble with r 4.4
# compensating for pak not seeming to understand ubuntu noble yet
if ( (identical(unname(machine), "aarch64")) 
	& (identical(unname(tolower(lsb_id)), "ubuntu"))
	& ((identical(unname(release), "22.04")) | (identical(unname(release), "24.04"))) 
	& (identical(unname(rver_short), "4.4")) )
{	
	# super duper hardcoded, unofficial, may go away without warning, supports jammy(22.04) and noble (24.04)
	# https://github.com/r-hub/repos
	# https://forum.posit.co/t/arm64-binary-packages-for-linux-in-posit-public-package-manager/178514/4
	# https://github.com/r-hub/repos?tab=readme-ov-file#ubuntu-2404--r-release-on-aarch64-ubuntu-2404-aarch64-r44
	message( paste("adding cran(like) repo for aarch64 (arm64) binaries for:", lsb_id, release, ", r:", rver_short))
	github_cran_binaries <- sprintf("https://raw.githubusercontent.com/r-hub/repos/main/%s-%s-%s/%s", tolower(lsb_id), release, machine, rver_short)
	pak::repo_add(github_binaries = github_cran_binaries)
		
} else if (identical(unname(machine), "x86_64"))
{
	# unclear what version(s) of r are required
	message( paste ("adding cran(like) repo for x86_64 (amd64) binaries for:", lsb_id, release, ", r:", rver_short))
	p3m_cran_binaries <- sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", codename)
	pak::repo_add(p3m_binaries = p3m_cran_binaries)
	
} else 
{ 
	message( paste("architecture:", machine, ", distro:", lsb_id, release, ", r:", rver_short, "unrecognized for cran binary package availability, compiling cran dependencies from source"))
}


# seems to understand arm64 and amd64, but only r >= 4.4
if 		( ((identical(unname(machine), "aarch64")) | (identical(unname(machine), "x86_64")))
	& 	 ((identical(unname(substr(rver, 1, 3)), "4.4")) | (identical(unname(rver_short), "4.5"))))
{
	message( paste("adding bioconductor repo for", machine, "binaries for:", lsb_id, release, ", r:", rver_short))
	r_universe_bioc_binaries <- sprintf('%s/bin/linux/%s/%s', "https://bioc.r-universe.dev", codename, rver_short)
	pak::repo_add(r_universe_bioc_binaries = r_universe_bioc_binaries)
	
} else
{
	message( paste("architecture:", machine, ", distro:", lsb_id, release, ", r:", rver_short, "unrecognized for bioconductor binary package availability, compiling bioconductor dependencies from source"))
}


pak::repo_status()

message("pkg_deps_tree")
pak::pkg_deps_tree("local::.", upgrade = FALSE, dependencies = TRUE)

# could be a flag to control this
#checkResults <- devtools::check(error_on = c("never"), env_vars = c(`_R_CHECK_TESTS_NLINES_` = '0', `CI` = 'true'))

# use devtools build/install instead...?
# can i install the package from the source package made by devtools::build?
message("pkg_install")
pak::pkg_install("local::.", upgrade = FALSE, dependencies = TRUE)


# only need to run this to make the docker image smaller
pak::pak_cleanup(force = TRUE)

installed.packages()
