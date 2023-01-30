FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>

COPY stat_markdown.yml .
RUN \
   conda env update -n root -f stat_markdown.yml \
&& conda clean -a

RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps 

RUN R -e "install.packages('vroom', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('purrr', repos='http://cran.us.r-project.org')"