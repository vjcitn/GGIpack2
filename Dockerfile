#ARG UBUNTU_RELEASE=noble
#FROM ubuntu:$UBUNTU_RELEASE

# might need the latest tag
# FROM ghcr.io/r-lib/rig/ubuntu-24.04
FROM ubuntu:noble

ARG R_VERSION=4.4.1

ENV DEBIAN_FRONTEND=noninteractive
ENV R_CLI_NUM_COLORS=256

# for some reason, libcurl4-openssl-dev is needed when building the arm version of this image
# maybe libcurl4-openssl-dev isn't included in the arm version of this image
RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
	rm -f /etc/apt/apt.conf.d/docker-clean && \
    apt update && \
	apt upgrade -y && \
	apt install -y lsb-release libcurl4-openssl-dev apt-utils curl build-essential pkg-config qpdf && \
	curl -L https://rig.r-pkg.org/deb/rig.gpg -o /etc/apt/trusted.gpg.d/rig.gpg && \
	sh -c 'echo "deb http://rig.r-pkg.org/deb rig main" > /etc/apt/sources.list.d/rig.list' && \
	apt update && \
	apt install -y r-rig && \
	rig install ${R_VERSION} && \
	rm -rf /tmp/rig	
	
# the pak installation of the arm version of the base image seems to be partially broken, so reinstalling it
# RUN Rscript --vanilla -e 'pak::pak_update(force = TRUE)'
	
# RUN Rscript --vanilla -e 'pak::pak_install_extra()'
# then clean pak cache...?

# https://devops.stackexchange.com/questions/13446/how-can-i-get-the-docker-target-platform-inside-the-build-environment-dockerfi
# need to clean /tmp/?
COPY . /opt/GGIpack2
WORKDIR /opt/GGIpack2
RUN chmod +x install_ggipack2.R
RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
	./install_ggipack2.R && \
	rm -rf /tmp/*
RUN chmod +x start_gtexapp.R
# cleanup tmp.  i think pak uses its own tmp directory in addition to the session tmp directory?
# RUN rm -rf /tmp/*

ENTRYPOINT [ "./start_gtexapp.R" ]