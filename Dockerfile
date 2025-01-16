FROM mambaorg/micromamba:1.5.8-bookworm
# This only installs dinosaur
LABEL description="Additional stuff and dinosaur which does not work in biocontainer due to lack of fontconfig"

ARG NEW_MAMBA_USER=mambauser
ARG NEW_MAMBA_USER_ID=1
ARG NEW_MAMBA_USER_GID=1
USER root

# to have envsubst, ps
RUN apt update && apt install -y gettext-base procps

# for dinosaur
RUN micromamba install -y -n base -c conda-forge -c bioconda dinosaur=1.2.0 pyarrow=18.1.0
