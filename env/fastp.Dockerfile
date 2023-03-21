FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>
RUN  apt-get --allow-releaseinfo-change update \
    && apt-get install -y g++ \
    procps \
    automake \
    autoconf \
    libpcre3-dev \
    libssl-dev \
    make \
    zlib1g-dev

RUN wget http://opengene.org/fastp/fastp
RUN chmod a+x ./fastp
