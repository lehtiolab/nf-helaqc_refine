FROM mambaorg/micromamba:1.5.8-bookworm
# This only installs dinosaur
LABEL description="DIA-NN, additional stuff and dinosaur which does not work in biocontainer due to lack of fontconfig"

ARG NEW_MAMBA_USER=mambauser
ARG NEW_MAMBA_USER_ID=1
ARG NEW_MAMBA_USER_GID=1
USER root

RUN apt update && apt upgrade -y
# to have envsubst, ps
RUN apt install -y gettext-base procps
# to get DIA-NN
RUN apt install -y wget unzip libgomp1 locales

# for dinosaur
RUN micromamba install -y -n base -c conda-forge -c bioconda dinosaur=1.2.0 pyarrow=18.1.0

# For DIA-NN
# Configure locale to avoid runtime errors
RUN sed -i 's/# en_US.UTF-8/en_US.UTF-8/' /etc/locale.gen && \
    locale-gen && \
    update-locale LANG=en_US.UTF-8

# Set environment variables for locale
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# Download DIA-NN version <version>
#RUN wget https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.3.1-Academia-Linux.zip -O /tmp/diann.Linux.zip
RUN wget https://github.com/vdemichev/DiaNN/releases/download/1.9.2/diann-1.9.2.Linux_update_2024-10-31.zip -O /tmp/diann.Linux.zip
# Unzip the DIA-NN package
RUN cd / && unzip /tmp/diann.Linux.zip
RUN rm /tmp/diann.Linux.zip

# Set appropriate permissions for the DIA-NN folder
RUN chmod -R 775 /diann-1.9.2
