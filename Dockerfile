# syntax=docker.io/docker/dockerfile:1.7-labs
ARG UBUNTU_RELEASE=focal
FROM ubuntu:$UBUNTU_RELEASE

ARG R_VERSION=4.4.0

ENV DEBIAN_FRONTEND=noninteractive
ENV R_CLI_NUM_COLORS=256
ENV BIOCONDUCTOR_USE_CONTAINER_REPOSITORY=TRUE

RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
	rm -f /etc/apt/apt.conf.d/docker-clean && \
    apt update && apt upgrade && apt install -y curl apt-utils 
	
RUN curl -L https://rig.r-pkg.org/deb/rig.gpg -o /etc/apt/trusted.gpg.d/rig.gpg

RUN sh -c 'echo "deb http://rig.r-pkg.org/deb rig main" > /etc/apt/sources.list.d/rig.list'

# keep this around - non-x86_64 platforms won't be getting the binaries of these required packages
# so we may need to revive this for alternative platforms sicne everything is going to have to get
# compiled from the primordial ooze
# also, we should probably make this a separate project so we can try to make multi-arch images?
#RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
#	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
#	rm -f /etc/apt/apt.conf.d/docker-clean && \
#    apt update && apt install -y r-rig build-essential qpdf tree libssl-dev pkg-config libcurl4-openssl-dev

RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
	rm -f /etc/apt/apt.conf.d/docker-clean && \ 
	apt update && apt install -y r-rig build-essential pkg-config qpdf && rig install ${R_VERSION}

# if dependencies=true, devtools installs biocmanager
# this seems to take up a quarter gigabyte of space for some reason?
# is is pretty cool, though
#RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
#	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
#	rm -f /etc/apt/apt.conf.d/docker-clean && \ 
#	Rscript --vanilla -e 'pak::pkg_deps_tree(c("pillar", "pkgcache", "devtools", "desc", "jsonlite"), upgrade = FALSE)'

#RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
#	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
#	rm -f /etc/apt/apt.conf.d/docker-clean && \ 
#	Rscript --vanilla -e 'pak::pkg_install(c("pillar", "pkgcache", "devtools", "desc", "jsonlite"), upgrade = FALSE)'	

#RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
#	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
#	rm -f /etc/apt/apt.conf.d/docker-clean && \ 
#	Rscript --vanilla -e 'pak::pkg_install(c("BiocManager"), upgrade = FALSE)'	
	
#RUN Rscript --vanilla -e 'installed.packages()'

#RUN Rscript --vanilla -e 'pak::cache_summary()'
# then clean pak cache...?

COPY --exclude=test_path/ --exclude=test_path.R . /tmp/GGIpack2
WORKDIR /tmp/GGIpack2
RUN chmod +x install_ggipack2.R
RUN ./install_ggipack2.R
COPY test_path/ /tmp/GGIpack2
COPY test_path.R /tmp/GGIpack2

