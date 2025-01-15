#ARG UBUNTU_RELEASE=noble
#FROM ubuntu:$UBUNTU_RELEASE

# might need the latest tag
FROM ghcr.io/r-lib/rig/ubuntu-24.04

ARG R_VERSION=4.4.1

ENV DEBIAN_FRONTEND=noninteractive
ENV R_CLI_NUM_COLORS=256

RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
	rm -f /etc/apt/apt.conf.d/docker-clean && \
    apt update && apt upgrade -y && apt install -y lsb-release && rig install ${R_VERSION}

# RUN Rscript --vanilla -e 'pak::pak_install_extra()'
# then clean pak cache...?

COPY . /tmp/GGIpack2
WORKDIR /tmp/GGIpack2
RUN chmod +x install_ggipack2.R
RUN ./install_ggipack2.R

