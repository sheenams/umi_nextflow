# Teeny-tiny matplotlib image based on alpine
FROM ewels/multiqc:latest

RUN apk add --no-cache --update \
    python3 python3-dev gcc \
    gfortran musl-dev g++ \
    libffi-dev openssl-dev \
    libxml2 libxml2-dev \
    libxslt libxslt-dev \
    libjpeg-turbo-dev zlib-dev

RUN pip install --upgrade pip

RUN pip install pandas

ENTRYPOINT [ "multiqc" ]
CMD [ "." ]