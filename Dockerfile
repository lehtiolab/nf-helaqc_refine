FROM nfcore/base
LABEL description="Docker image containing all requirements for lehtiolab/helaqc pipeline"

COPY environment.yml /

RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/helaqc-2.0/bin:$PATH
RUN git clone https://github.com/glormph/msstitch /msstitch 
RUN cd /msstitch && git pull && git checkout a36b08b689fb33276bbc269315fb5bb0561f9f6c && pip install -e /msstitch
