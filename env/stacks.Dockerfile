FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>

#COPY bbtools.yml .
#RUN \
#   conda env update -n root -f bbtools.yml \
#&& conda clean -a

RUN conda install -c bioconda stacks
RUN conda install -c "bioconda/label/cf201901" stacks

RUN conda install -c bioconda seqtk
RUN conda install -c "bioconda/label/cf201901" seqtk