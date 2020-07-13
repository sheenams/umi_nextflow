FROM ubuntu:20.04

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev git \
  && apt-get clean

WORKDIR /usr/local/bin
RUN ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip \
  && pip install pandas pysam

CMD [ "/bin/bash" ]
