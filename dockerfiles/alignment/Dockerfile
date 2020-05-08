FROM ubuntu:20.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        zlib1g-dev \
        libncurses5-dev \
        librust-bzip2-dev \
        libcurl4-openssl-dev \
        pkg-config \
        liblzma-dev \
        gcc \
        git \
        curl \
        ca-certificates \
        make && \
    apt-get clean

ENV SAMTOOLS_VERSION 1.10
RUN curl -SL https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2 \
    | tar -xjC /tmp/ \
    && cd /tmp/samtools-$SAMTOOLS_VERSION \
    && ./configure \
    && make \
    && make install

ENV BWA_VERSION 0.7.17
RUN curl -SL http:\//downloads.sourceforge.net/project/bio-bwa/bwa-$BWA_VERSION.tar.bz2 \
    | tar -xjC /tmp/ \
    && cd /tmp/bwa-$BWA_VERSION \
    && make \
    && mv /tmp/bwa-$BWA_VERSION/bwa /usr/bin


ENV SEQTK_VERSION 1.3
RUN curl -SL https://github.com/lh3/seqtk/archive/v${SEQTK_VERSION}.tar.gz \
    | tar -xzC /tmp/ \
    && cd /tmp/seqtk-${SEQTK_VERSION} \
    && make