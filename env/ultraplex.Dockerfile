FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>

COPY ultraplex.yml .
RUN \
   conda env update -n root -f ultraplex.yml \
&& conda clean -a

RUN  apt-get --allow-releaseinfo-change update \
    && apt-get install -y g++

RUN pip install ultraplex