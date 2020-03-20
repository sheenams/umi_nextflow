#!/nin/bash

bcl2fastq="singularity exec --bind $(pwd) --pwd $(pwd) /mnt/disk2/com/container-images/bcl2fastq-v2.20.simg bcl2fastq "

#READ_STRUCTURE='101T8B9M8B101T'
INPUT_DIR=$1;
SAMPLESHEET=$2;
OUTPUT_DIR=$3;

#Run info has 101,17,8,101
#SampleSheet does NOT have N in barcode
#Current testing:
#191212_NB0330_MiniOncoKAPA263R-MONCv1-bcl2fastq-1-30 
#200103_HA0717_MiniOncoKAPA266R-MONCv1-bcl2fastq-1-30

echo $INPUT_DIR $SAMPLESHEET $OUTPUT_DIR
 $bcl2fastq --runfolder-dir $INPUT_DIR \
            --output-dir $OUTPUT_DIR \
            --sample-sheet $SAMPLESHEET \
            --barcode-mismatches 0 \
 	   --processing-threads 40 \
            --with-failed-reads \
 	   --use-bases-mask Y101,I8Y9,I8,Y101 \
            --ignore-missing-bcls \
 	   --mask-short-adapter-reads 0 \
            --no-lane-splitting ;

for sample in `cut -f1 -d',' $SAMPLESHEET |awk '{if(NR>2)print}'`;do
    R1=$sample*_R1_*
    UMI=$sample*_R2_*
    R2=$sample*_R3_*
    paste \
        <(zcat ${OUTPUT_DIR}/*/$R1) \
        <(zcat ${OUTPUT_DIR}/*/$UMI) \
        | awk 'NR%4==1{readname=$1}
               NR%4==2{seq=$1; umi=$2}
               NR%4==0 {print readname " RX:Z:"umi"\n"seq"\n+\n"$1;}' \
		   | gzip > ${OUTPUT_DIR}/${sample}.1.fastq.gz ;
    paste \
        <(zcat ${OUTPUT_DIR}/*/$R2) \
        <(zcat ${OUTPUT_DIR}/*/$UMI) \
        | awk 'NR%4==1{readname=$1}
               NR%4==2{seq=$1; umi=$2}
               NR%4==0 {print readname " RX:Z:"umi"\n"seq"\n+\n"$1;}' \
		   | gzip > ${OUTPUT_DIR}/${sample}.2.fastq.gz

done
