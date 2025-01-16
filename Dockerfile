#ARG UBUNTU_RELEASE=noble
#FROM ubuntu:$UBUNTU_RELEASE

# might need the latest tag
FROM ghcr.io/r-lib/rig/ubuntu-24.04

ARG R_VERSION=4.4.1

ENV DEBIAN_FRONTEND=noninteractive
ENV R_CLI_NUM_COLORS=256

# for some reason, libcurl4-openssl-dev is needed when building the arm version of this image
# maybe libcurl4-openssl-dev isn't included in the arm version of this image
RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
	rm -f /etc/apt/apt.conf.d/docker-clean && \
    apt update && apt upgrade -y && apt install -y lsb-release libcurl4-openssl-dev && rig install ${R_VERSION}

# the pak installation of the arm version of tis base image seems to be broken, so reinstalling it
RUN Rscript --vanilla -e 'pak::pak_update(force = TRUE)'
	
# RUN Rscript --vanilla -e 'pak::pak_install_extra()'
# then clean pak cache...?
# pak seems to be bad on the arm64 image installation for some reason and must be reinstalled
# https://devops.stackexchange.com/questions/13446/how-can-i-get-the-docker-target-platform-inside-the-build-environment-dockerfi
COPY . /tmp/GGIpack2
WORKDIR /tmp/GGIpack2
RUN chmod +x install_ggipack2.R
RUN ./install_ggipack2.R
RUN chmod +x start_gtexapp.R

ENTRYPOINT [ "./start_gtexapp.R" ]