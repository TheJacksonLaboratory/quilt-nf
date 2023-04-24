FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>


COPY QUILT.yml .
RUN \
   conda env update -n root -f QUILT.yml \
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
    xterm \
    x11-utils \
    libcairo2-dev \
    libblas-dev \
    libssh2-1-dev \
    libgit2-dev

#RUN R -e "install.packages('parallel', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('Rcpp', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('RcppArmadillo', repos='http://cran.us.r-project.org')"
RUN R -e "library('QUILT')"
