FROM continuumio/miniconda3:4.10.3p0

RUN conda install mamba -n base -c conda-forge
COPY environment.yml /
RUN mamba env create -f /environment.yml

ENV PATH /opt/conda/envs/vsd-1.0/bin:$PATH

