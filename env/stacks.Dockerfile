FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>

ADD https://anaconda.org/bioconda/stacks/2.61/download/linux-64/stacks-2.61-hd03093a_1.tar.bz2 ./
RUN tar xfvz stacks-2.64.tar.gz \
        cd stacks-2.64 \
        ./configure \
        make

#RUN conda install -c bioconda stacks
#RUN conda install -c "bioconda/label/main" stacks-2.61
RUN conda install -c bioconda seqtk
RUN conda install -c "bioconda/label/cf201901" seqtk
