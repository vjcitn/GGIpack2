FROM chanart.bwh.harvard.edu:5002/jenkins/r-support-docker:latest

# currently your rig options for this image are:
# rig default 4.0.5
# rig default 4.1.3
# rig default 4.2.3
# rig default 4.3.3
# rig default 4.4.2
# rig default release  (same as 4.4.2)

COPY . /opt/GGIpack2
WORKDIR /opt/GGIpack2
RUN --mount=type=cache,id=apt_cache,target=/var/cache/apt,sharing=locked \
	--mount=type=cache,id=apt_lib,target=/var/lib/apt,sharing=locked \
	rig default release && \
	rig ls && \
	/opt/sdlc/build.R --do_install TRUE --fail_on_check_error TRUE && \
	rm -rf /tmp/*
RUN chmod +x start_shinyapp.R
