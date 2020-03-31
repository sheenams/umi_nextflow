FROM nkrumm/alignment

WORKDIR /

RUN apt-get install make

ADD https://github.com/alexpreynolds/sample/archive/v1.0.3.tar.gz sample.v1.0.3
RUN cd /sample.v1.0.3/sample-1.0.3 \
	&& make \
	&& cp /sample.v1.0.3/sample-1.0.3/sample /home/biodocker/bin/