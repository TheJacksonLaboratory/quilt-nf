FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>


COPY umi_tools.yml .
RUN \
   conda env update -n root -f umi_tools.yml \
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
    libgit2-dev

RUN conda install -c bioconda -c conda-forge umi_tools
