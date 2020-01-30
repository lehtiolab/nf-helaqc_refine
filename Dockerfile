FROM nfcore/base
LABEL description="Docker image containing all requirements for lehtiolab/helaqc pipeline"

COPY environment.yml /

RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/helaqc-2.0/bin:$PATH

RUN apt update && apt install -y fontconfig && apt clean -y
