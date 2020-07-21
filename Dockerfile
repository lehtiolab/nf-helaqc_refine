FROM nfcore/base
LABEL description="Docker image containing all requirements for lehtiolab/helaqc pipeline"

COPY environment.yml /

RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/helaqc-2.1/bin:$PATH
