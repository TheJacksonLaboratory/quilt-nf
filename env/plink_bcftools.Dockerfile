FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>


COPY plink_bcftools.yml .
RUN \
   conda env update -n root -f plink_bcftools.yml \
&& conda clean -a

RUN  apt-get --allow-releaseinfo-change update \
    && apt-get install -y procps \
    ssh \
    bash \
    pkg-config \
    libglpk-dev \
    libz-dev \
    tk \
    libxml2 \
    libxml2-dev \
    libbz2-dev \
    liblzma-dev \
    libblas-dev \
    libssh2-1-dev \
    libgit2-dev \
    zlib

RUN conda install -c bioconda plink
RUN conda install -c bioconda bcftools
