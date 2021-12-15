FROM continuumio/miniconda3:4.10.3p0

RUN conda install mamba -n base -c conda-forge
COPY environment.yml /
RUN mamba env create -f /environment.yml

RUN echo "conda activate vsd-1.0" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

COPY docker-entrypoint.sh /
ENTRYPOINT ["/docker-entrypoint.sh"]
