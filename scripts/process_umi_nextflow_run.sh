#!/bin/bash

CONFIG_DIR=$1
NEXTFLOW_PATH=/mnt/disk10/users/bricegc/umi_nextflow/nextflow
NF_FILE=/mnt/disk10/users/bricegc/umi_nextflow/main.nf

if [ ! -d "$CONFIG_DIR" ]; then
    echo "Did not provide valid directory, exiting."
    exit 1
fi

for config_file in $CONFIG_DIR/*.config; do
  (  NF_CMD="$NEXTFLOW_PATH run $NF_FILE -c $config_file -profile singularity"
    echo $NF_CMD 
    $NF_CMD ) &
done
