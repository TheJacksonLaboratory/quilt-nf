FROM rocker/r-ver:4.0.0
LABEL Sam Widmayer <sjwidmay@gmail.com>

RUN R -e "install.packages('rmarkdown', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('vroom', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('purrr', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('dplyr', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('tidyr', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('plotly', repos='http://cran.us.r-project.org')"