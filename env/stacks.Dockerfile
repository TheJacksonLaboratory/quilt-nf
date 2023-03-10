FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>

#RUN conda install -c bioconda stacks
RUN conda install -c "bioconda/label/main" stacks-2.61
RUN conda install -c bioconda seqtk
RUN conda install -c "bioconda/label/cf201901" seqtk
