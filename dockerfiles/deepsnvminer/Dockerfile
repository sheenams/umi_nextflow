FROM ubuntu:20.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        zlib1g-dev \
        libncurses5-dev \
        librust-bzip2-dev \
        liblzma-dev \
        gcc \
        git \
        curl \
        ca-certificates \
        make && \
    apt-get clean

ENV SAMTOOLS_VERSION 0.1.19
RUN curl -SL https://sourceforge.net/projects/samtools/files/samtools/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    | tar -xjC /tmp/ \
    && cd /tmp/samtools-$SAMTOOLS_VERSION \
    && make \
    && mv samtools /usr/bin 


ENV BWA_VERSION 0.7.17
RUN curl -SL http:\//downloads.sourceforge.net/project/bio-bwa/bwa-$BWA_VERSION.tar.bz2 \
    | tar -xjC /tmp/ \
    && cd /tmp/bwa-$BWA_VERSION \
    && make \
    && mv /tmp/bwa-$BWA_VERSION/bwa /usr/bin

RUN git clone https://github.com/mattmattmattmatt/DeepSNVMiner.git 

# CMD [ "DeepSNVMiner/run_deepseq.pl" ]