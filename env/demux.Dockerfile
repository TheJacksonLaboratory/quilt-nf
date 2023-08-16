FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>

COPY demux.yml .
RUN \
   conda env update -n root -f demux.yml \
&& conda clean -a

RUN  apt-get --allow-releaseinfo-change update \
    && apt-get install -y g++

RUN pip install demultiplex
