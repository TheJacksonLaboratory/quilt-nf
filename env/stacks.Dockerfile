FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>

RUN conda install -c "bioconda/label/cf201901" stacks
RUN conda install -c bioconda seqtk
RUN conda install -c "bioconda/label/cf201901" seqtk
