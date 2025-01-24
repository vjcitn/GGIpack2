ARG UBUNTU_RELEASE=noble
FROM ubuntu:$UBUNTU_RELEASE

ARG R_VERSION=4.4.2

ENV DEBIAN_FRONTEND=noninteractive
ENV R_CLI_NUM_COLORS=256

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
	rm -rf /tmp/*	

COPY . /opt/GGIpack2
WORKDIR /opt/GGIpack2
	
# install my utils package
# RUN Rscript --vanilla -e 'pak::pkg_install("url::https://chanart.bwh.harvard.edu:443/artifactory/cran-dev-local/src/contrib/cdnmsdlcutils_0.0.4.tar.gz")'
# RUN Rscript --vanilla -e 'install.packages("cdnmsdlcutils", repos=c("https://chanart.bwh.harvard.edu/artifactory/cran-dev-local"))'

# pak is self-updating
RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
	Rscript --vanilla -e 'pak::repo_add(chanart_prod = "https://chanart.bwh.harvard.edu/artifactory/cran-prod-local") ; pak::pkg_install("cdnmsdlcutils")' && \
	Rscript --vanilla -e 'cdnmsdlcutils::sdlc_install_support()' && \
	rm -rf /tmp/*	

RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
	mkdir -p /opt/target && \
	Rscript --vanilla -e 'cdnmsdlcutils::sdlc_build(output_dir="/opt/target", do_install=TRUE, fail_on_check_error=TRUE)' && \
	rm -rf /tmp/* && \
	chmod +x start_shinyapp.R

ENTRYPOINT [ "./start_shinyapp.R" ]