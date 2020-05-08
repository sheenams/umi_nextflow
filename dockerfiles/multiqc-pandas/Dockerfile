FROM ubuntu:20.04

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev git \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip \
  apt-get clean

RUN pip install pandas

#RUN pip install multiqc
RUN pip install git+https://github.com/nkrumm/multiqc.git@fgbio-relative-percent-view

CMD [ "/bin/bash" ]
