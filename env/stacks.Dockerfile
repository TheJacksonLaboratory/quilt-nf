FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>


RUN  apt-get --allow-releaseinfo-change update \
    && apt-get install -y g++ \
    apt-get install -y procps \
    automake \
    autoconf \
    libpcre3-dev \
    libssl-dev \
    make \
    zlib1g-dev

RUN wget https://catchenlab.life.illinois.edu/stacks/source/stacks-2.64.tar.gz --no-check-certificate
RUN tar xfvz stacks-2.64.tar.gz

WORKDIR /stacks-2.64
RUN chmod 775 ./configure && ./configure && make
