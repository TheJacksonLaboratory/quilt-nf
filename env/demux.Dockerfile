FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>


RUN  apt-get --allow-releaseinfo-change update \
    && apt-get install -y g++ \
    python3-pip

RUN git clone https://github.com/jfjlaros/demultiplex
RUN cd demultiplex
RUN pip install .
