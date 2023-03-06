FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>

#COPY bbtools.yml .
#RUN \
#   conda env update -n root -f bbtools.yml \
#&& conda clean -a

RUN conda install -c bioconda bbmap
RUN conda install -c "bioconda/label/cf201901" bbmap
