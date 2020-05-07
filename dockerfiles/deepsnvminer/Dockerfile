FROM ubuntu:20.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        zlib1g-dev \
        libncurses5-dev \
        gcc \
        git \
        curl \
        ca-certificates \
        make && \
    apt-get clean

ENV BWA_VERSION 0.7.17
RUN curl -SL http:\//downloads.sourceforge.net/project/bio-bwa/bwa-$BWA_VERSION.tar.bz2 \
    | tar -xjC /tmp/ \
    && cd /tmp/bwa-$BWA_VERSION \
    && make \
    && mv /tmp/bwa-$BWA_VERSION/bwa /usr/bin

ENV SAMTOOLS_VERSION 1.3.1
RUN curl -SL https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2 \
    | tar -xjC /tmp/ \
    && cd /tmp/samtools-$SAMTOOLS_VERSION \
    && make \
    && make install

RUN rm -rf /tmp/bwa-$BWA_VERSION /tmp/samtools-$SAMTOOLS_VERSION \
    && apt-get clean \
    && apt-get remove --yes --purge build-essential gcc-multilib apt-utils zlib1g-dev

 RUN git clone https://github.com/mattmattmattmatt/DeepSNVMiner.git && \
    cd DeepSNVMiner