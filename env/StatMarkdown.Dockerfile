FROM rocker/r-ver:4.0.0
LABEL Sam Widmayer <sjwidmay@gmail.com>

#FROM continuumio/miniconda
#COPY stat_markdown.yml .
#RUN \
#   conda env update -n root -f stat_markdown.yml \
#&& conda clean -a

RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps
RUN apt-get install libcurl4-openssl-dev libxml2-dev libssh-dev x11-apps -y
RUN R -e "install.packages('vroom', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('purrr', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('dplyr', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('tidyr', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('knitr', repos='http://cran.us.r-project.org', dependencies = TRUE)"
RUN R -e "install.packages('png', repos='http://cran.us.r-project.org', dependencies = TRUE)"
RUN R -e "install.packages('rmarkdown', repos='http://cran.us.r-project.org', dependencies = TRUE)"
RUN R -e "install.packages('pandoc', repos='http://cran.us.r-project.org', dependencies = TRUE)"