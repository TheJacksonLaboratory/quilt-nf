FROM rocker/rstudio:4.1.0
LABEL Sam Widmayer <sjwidmay@gmail.com>

#FROM continuumio/miniconda
#COPY stat_markdown.yml .
#RUN \
#   conda env update -n root -f stat_markdown.yml \
#&& conda clean -a

RUN  apt-get --allow-releaseinfo-change update \
    && apt-get install -y procps \
    ssh \
    bash \
    pkg-config \
    libglpk-dev \
    libjpeg62-dev \
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
    
RUN R -e "install.packages('vroom', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('purrr', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('dplyr', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('tidyr', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('knitr', repos='http://cran.us.r-project.org', dependencies = TRUE)"
RUN R -e "install.packages('png', repos='http://cran.us.r-project.org', dependencies = TRUE)"
RUN R -e "install.packages('rmarkdown', repos='http://cran.us.r-project.org', dependencies = TRUE)"
RUN R -e "install.packages('pandoc', repos='http://cran.us.r-project.org', dependencies = TRUE)"
