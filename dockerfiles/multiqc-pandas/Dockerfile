FROM ubuntu:18.04

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

RUN pip install pandas

#RUN pip install multiqc
RUN apt-get install -y git
RUN pip install git+https://github.com/nkrumm/multiqc.git@fgbio-relative-percent-view

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

ENTRYPOINT [ ]
CMD [ "/bin/bash" ]
