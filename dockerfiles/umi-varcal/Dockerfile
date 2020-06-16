FROM ubuntu:20.04

RUN apt-get update && \
    apt-get install -y git \
    python3 \
    libbz2-dev \
    python3-pip && \
    apt-get clean && apt-get purge
   
RUN python3 -m pip install msgpack==0.5.6 \
    cython \
    pyfasta \
    pysam \
    scipy \
    statsmodels

RUN mkdir /home/umi-varcal &&\
    git clone --recursive https://gitlab.com/vincent-sater/umi-varcal_v3.git /home/umi-varcal && \
    cd /home/umi-varcal/functions && \
    rm -r functions.so functions.c* &&\
    python3 setup.py build_ext --inplace

WORKDIR /home/umi-varcal

CMD [ "python3", "./umi-varcal.py" ]