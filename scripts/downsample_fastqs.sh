#!/bin/bash
# Downsamples fastqs in a directory to 1,400,000 reads each

RUN_DIR=$1
OUTPUT_DIR=$2

if [ ! -d "$RUN_DIR" ]; then
    echo "Did not provide valid run directory, exiting."
    exit 1
fi

if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Did not provide valid run directory, exiting."
    exit 1
fi

for fastq_file in $RUN_DIR/*.fastq.gz; do
    echo $fastq_file
    fname=$(basename $fastq_file)
    echo $fname
    seqtk sample -s100 $fastq_file 1400000 > $OUTPUT_DIR/downsample_$fname
done
