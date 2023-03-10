FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>

#RUN apt-get install -y autotools-dev autoconf

RUN wget https://catchenlab.life.illinois.edu/stacks/source/stacks-2.64.tar.gz --no-check-certificate
RUN tar xfvz stacks-2.64.tar.gz
RUN cd stacks-2.64
RUN sh configure
RUN make

#RUN conda install -c bioconda stacks
#RUN conda install -c "bioconda/label/main" stacks-2.61
RUN conda install -c bioconda seqtk
RUN conda install -c "bioconda/label/cf201901" seqtk
