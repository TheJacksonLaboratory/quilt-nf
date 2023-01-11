FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>


COPY STITCH.yml .
RUN \
   conda env update -n root -f STITCH.yml \
&& conda clean -a

RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps  

RUN R -e "install.packages('parallel', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('Rcpp', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('RcppArmadillo', repos='http://cran.us.r-project.org')"
RUN R -e "library('STITCH')"