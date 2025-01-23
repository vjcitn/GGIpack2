FROM r-support-docker

COPY . /opt/GGIpack2
WORKDIR /opt/GGIpack2
RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
	/opt/sdlc/build.R --do_install TRUE --fail_on_check_error TRUE && \
	rm -rf /tmp/*
RUN chmod +x start_shinyapp.R

ENTRYPOINT [ "./start_shinyapp.R" ]