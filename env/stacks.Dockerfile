FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>

COPY stacks.yml .
RUN \
   conda env update -n root -f stacks.yml \
&& conda clean -a

RUN conda install -c bioconda seqtk
RUN conda install -c "bioconda/label/cf201901" seqtk
