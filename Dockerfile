ARG UBUNTU_RELEASE=noble
FROM ubuntu:$UBUNTU_RELEASE

ARG R_VERSION=4.4.2

ENV DEBIAN_FRONTEND=noninteractive
ENV R_CLI_NUM_COLORS=256
ENV _R_SHLIB_STRIP_=1

# stripping out the debug symbols from all the builtin and installed libraries saves ~20% space?
# also, i don't think _R_SHLIB_STRIP_ is doing much if anything, but that might be because we're installing alot of binaries'

RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
	rm -f /etc/apt/apt.conf.d/docker-clean && \
    apt update && \
	apt upgrade -y && \
	apt install -y apt-utils build-essential curl libcurl4-openssl-dev lsb-release pkg-config qpdf && \
	curl -L https://rig.r-pkg.org/deb/rig.gpg -o /etc/apt/trusted.gpg.d/rig.gpg && \
	sh -c 'echo "deb http://rig.r-pkg.org/deb rig main" > /etc/apt/sources.list.d/rig.list' && \
	apt update && \
	apt install -y r-rig && \
	rig install ${R_VERSION} && \
#	find / -name "*.so*" | xargs strip --strip-debug || true && \
#	find / -name "*.o" | xargs strip --strip-debug || true && \
	rm -rf /tmp/*	

COPY . /opt/GGIpack2
WORKDIR /opt/GGIpack2
	
RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
	Rscript --vanilla -e 'pak::repo_add(chanart_prod = "https://chanart.bwh.harvard.edu/artifactory/cran-prod-local") ; pak::pkg_install("cdnmsdlcutils@0.0.4")' && \
	Rscript --vanilla -e 'cdnmsdlcutils::sdlc_install_support()' && \
	mkdir -p /opt/target && \
	Rscript --vanilla -e 'cdnmsdlcutils::sdlc_build(output_dir="/opt/target", do_install=TRUE, fail_on_check_error=TRUE)' && \
#	find / -name "*.so*" | xargs strip --strip-debug --verbose || true && \
#	find / -name "*.o" | xargs strip --strip-debug --verbose || true && \
	rm -rf /tmp/* && \
	chmod +x start_shinyapp.R

ENTRYPOINT [ "./start_shinyapp.R" ]