#!/usr/bin/env nextflow

// Setup the various inputs, defined in nexflow.config
fastq_pair_ch = Channel.fromFilePairs(params.input_folder + '*{1,2}.fastq.gz', flat: true)
                       .last(2)
                       .view()
                       .into{align_input; fastqc_ch}

// Assay specific files
picard_targets = file(params.picard_targets)
picard_baits = file(params.picard_baits)
bed_file = file(params.bed)

// Reference genome is used multiple times
reference_fasta = file(params.ref_fasta)
reference_index = Channel.fromPath(params.ref_index).into { 
  bwa_ref_index;
  bwa_realign_ref_index;
  picard_ref_index;
  qc_ref_index;
  filter_con_ref_index
}


process bwa {
   // Align fastqs
   label 'bwa'
   tag "${sample_id}"

   input:
     file(reference_fasta) from reference_fasta
     file("*") from bwa_ref_index.collect()
     tuple sample_id, file(fastq1), file(fastq2) from align_input
 
   output:
     tuple val(sample_id), file('*.bam') into align_ch
     tuple val(sample_id), val("standard"), file('*.standard.bam'), file('*.bai') into qc_standard_bam

   publishDir params.output, overwrite: true

   cpus 8
   
   script:
   // bwa mem options:
   // -K seed, -C pass tags from FASTQ -> alignment, -Y recommended by GATK?
   """
   bwa mem \
     -R'@RG\\tID:${sample_id}\\tSM:${sample_id}' \
     -K 10000000 \
     -C \
     -Y \
     -t${task.cpus}  \
     ${reference_fasta} ${fastq1} ${fastq2} 2> log.txt \
   | samtools sort -t${task.cpus} -m4G - -o ${sample_id}.standard.bam
   
   samtools index ${sample_id}.standard.bam
   """
} 

process sort_bam {
   //  Sort alignment by query name
   label 'picard'
   tag "${sample_id}"

   input: 
     tuple val(sample_id), file(bam) from align_ch

   output:
     tuple val(sample_id), file('*.sorted.bam') into set_mate_ch
     tuple val(sample_id), val("sorted"), file('*.sorted.bam') into temp_qc_sorted_bam

   memory "32GB"

   script:
   """
   picard -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./ \
   SortSam \
   I=${bam} \
   O=${sample_id}.sorted.bam \
   SORT_ORDER=queryname
   """
 }

process fgbio_setmateinformation{
  //  Adds and/or fixes mate information on paired-end reads
  // tuples the MQ (mate mapping quality), MC (mate cigar string), 
  // ensures all mate-related flag fields are tuple correctly, 
  // and that the mate reference and mate start position are correct.
   label "fgbio"
   tag "${sample_id}"

   input: 
    tuple sample_id, path(bam) from set_mate_ch

   output:
    tuple sample_id, "${sample_id}.mateinfo.bam" into mate_info_bam_ch
    tuple val(sample_id), val("set_mate"), file('*.bam') into temp_qc_mate_bam

   memory "32G"

   publishDir params.output, overwrite: true

   script:
   """
   fgbio -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./\
   SetMateInformation \
   --input ${bam} \
   --output ${sample_id}.mateinfo.bam
   """
 }

process fgbio_group_umi {
   // Groups reads together that appear to have come from the same original molecule.
   label "fgbio"
   tag "${sample_id}"

   input: 
    tuple sample_id, path(bam) from mate_info_bam_ch

   output:
     tuple val(sample_id), file('*.grpumi.bam') into grp_umi_bam_ch
     file('*.grpumi.histogram') into histogram_ch

   memory "32G"

   publishDir params.output, overwrite: true

   script:
   """
   fgbio -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./\
   GroupReadsByUmi \
   --input=${bam} \
   --output=${sample_id}.grpumi.bam \
   --family-size-histogram=${sample_id}.grpumi.histogram \
   --strategy=adjacency
   """
 }

process fgbio_callconsensus{
   //   Combined each set of reads to generate consensus reads
   //   1. base qualities are adjusted
   //   2. consensus sequence called for all reads with the same UMI, base-by-base.
   //   3. consensus raw base quality is modified by incorporating the probability of an error prior to
   //   calls each end of a pair independently, and does not jointly call bases that overlap within a pair. Insertion or deletion
   //   errors in the reads are not considered in the consensus model.integrating the unique molecular tags
   label 'fgbio'
   tag "${sample_id}"

   input: 
    tuple sample_id, path(bam) from grp_umi_bam_ch

   output:
     tuple val(sample_id), file('*.consensus.bam') into (consensus_bam_ch, qc_consensus_bam)
  
   memory "32G"

   publishDir params.output, overwrite: true

   script:
   """
   fgbio -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./\
   CallMolecularConsensusReads \
   --input=${bam} \
   --output=${sample_id}.consensus.bam \
   --min-reads=1 \
   --read-name-prefix=${sample_id} \
   --read-group-id=${sample_id} \
   --error-rate-pre-umi=45 \
   --error-rate-post-umi=40 \
   --min-input-base-quality=10 \
   --output-per-base-tags=true \
   --sort-order=queryname
   """
 }

process fgbio_filterconsensus{
   label "fgbio"
   tag "${sample_id}"

   input: 
     file(reference_fasta) from reference_fasta
     file("*") from filter_con_ref_index.collect()
     tuple val(sample_id), file(bam) from consensus_bam_ch

   output:
     tuple val(sample_id), file('*.filtered_consensus.bam') into (filter_consensus_bam_ch, consensus_fastq_ch, qc_filtered_consensus_bam)

   memory "32G"

   publishDir params.output, overwrite: true

   script:
   """
   fgbio -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./\
   FilterConsensusReads \
   --input=${bam} \
   --output=${sample_id}.filtered_consensus.bam \
   --ref=${reference_fasta} \
   --min-reads=2 \
   --max-read-error-rate 0.05 \
   --min-base-quality 10 \
   --max-base-error-rate 0.1 \
   --max-no-call-fraction 0.1
   """
 }

process sort_filter_bam {
   // Sort alignment by query name
   label "picard"
   tag "${sample_id}"

   input: 
    tuple sample_id, path(bam) from filter_consensus_bam_ch

   output:
    //tuple sample_id, "${sample_id}.sorted_filtered.bam" into sorted_filter_consensus_ch 
    tuple sample_id, "${sample_id}.sorted_consensus.bam" into (sorted_consensus_ch, sorted_consensus_fastq_ch,temp_qc_sorted_filter_bam) 
  
   publishDir params.output, overwrite: true
   
   memory "32G"

   script:
   """
   picard -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./ \
   SortSam \
   I=${bam} \
   O=${sample_id}.sorted_consensus.bam \
   SORT_ORDER=queryname
   """
 }

process bam_to_fastqs {
   label 'picard'
   tag "${sample_id}"

   input:
    tuple sample_id, path(bam) from sorted_consensus_fastq_ch

   output:
    tuple sample_id, "${sample_id}.fastq" into consensus_fastq
  
   memory "32G"

   script:
   """
   picard -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./ \
   SamToFastq \
   I=${bam} \
   FASTQ=${sample_id}.fastq \
   INTERLEAVE=true
   """
 }

 process realign_consensus {
   //-p Assume the first input query file is interleaved paired-end FASTA/Q.
   //-Y use soft clipping for supplementary alignment
   //-K process INT input bases in each batch regardless of nThreads (for reproducibility)
   label "bwa"
   tag "${sample_id}"

   input:
     file(reference_fasta) from reference_fasta
     file("*") from bwa_realign_ref_index.collect()
     tuple val(sample_id), file(fastq) from consensus_fastq

   output:
    tuple sample_id, "${sample_id}.realigned.bam" into realign_ch
    tuple val(sample_id), val("realigned"), file('*realigned.bam'), file('*.bai') into qc_sorted_final_bam

   cpus 8 

   script:
   """
   bwa mem \
   -R "@RG\\tID:${sample_id}\\tSM:${sample_id}" \
   -K 10000000 \
   -p \
   -Y \
   -t ${task.cpus} \
   ${reference_fasta} ${fastq} 2> log.txt \
   | samtools sort -t${task.cpus} -m4G - -o ${sample_id}.realigned.bam
   
   samtools index ${sample_id}.realigned.bam
   """
 }

 process sort_realign_bam {
   //  Sort alignment by query name
   label 'picard'
   tag "${sample_id}"

   input: 
    tuple sample_id, path(bam) from realign_ch

   output:
    tuple sample_id, "${sample_id}.sorted.bam" into (sorted_realign_consensus_ch, temp_qc_realign_bam)
  
   memory "32G"
  
   script:
   """
   picard -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./ \
   SortSam \
   I=${bam} \
   O=${sample_id}.sorted.bam \
   SORT_ORDER=queryname
   """
 }

 merge_ch = sorted_realign_consensus_ch.join(sorted_consensus_ch, remainder: true)
 process final_bam {
  //Merge consensus bam (unaligned) with aligned bam
   label 'picard'
   tag "${sample_id}"

   input:
    file(reference_fasta) from reference_fasta
    file("*") from picard_ref_index.collect()
    tuple sample_id, path(sorted_bam), path(sorted_filtered_bam) from merge_ch

   output:
     tuple val(sample_id), val("final"), file('*.final.bam') into qc_final_bam
  
   publishDir params.output, overwrite: true
   memory "32G"

   script:
   """
   picard -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./ \
   MergeBamAlignment \
   UNMAPPED=${sorted_filtered_bam} \
   ALIGNED=${sorted_bam} \
   O=${sample_id}.final.bam \
   R=${reference_fasta} \
   VALIDATION_STRINGENCY=SILENT \
   SORT_ORDER=coordinate
   """
 }

// ######### QC ######### //

// Combine channels for quality here
// into a single quality_ch for processing below.

//qc_consensus_bam,
//qc_filtered_consensus_bam

qc_sorted_final_bam.mix(qc_standard_bam)
            .into{ hs_metrics_ch; mosdepth_qc_ch } 


process quality_metrics {
   label 'picard'
   tag "${sample_id}"

   input:
     file(picard_targets) from picard_targets
     file(picard_baits) from picard_baits
     file(reference_fasta) from reference_fasta
     file("*") from qc_ref_index.collect()
     tuple val(sample_id), val(bam_type), file(bam), file(bai) from hs_metrics_ch

   output:
     path("${sample_id}.${bam_type}.hs_metrics") into hs_metrics_out_ch
     path("${sample_id}.${bam_type}.insert_size_metrics") into insert_size_metrics_ch

   publishDir params.output, overwrite: true
   
   memory "32G"

   script:
   """
 
   picard -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./ \
   CollectHsMetrics \
   TARGET_INTERVALS=${picard_targets} \
   BAIT_INTERVALS=${picard_baits} \
   REFERENCE_SEQUENCE=${reference_fasta} \
   INPUT=${bam} \
   OUTPUT=${sample_id}.${bam_type}.hs_metrics 

   picard -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./ \
   CollectInsertSizeMetrics \
   INCLUDE_DUPLICATES=true \
   INPUT=${bam} \
   OUTPUT=${sample_id}.${bam_type}.insert_size_metrics \
   HISTOGRAM_FILE=${sample_id}.${bam_type}.insert_size_histogram.pdf
   """
}

process fastqc {
  label 'fastqc'
  tag "${sample_id}"

  input:
    tuple sample_id, file(fastq1), file(fastq2) from fastqc_ch

  output:
    path "fastqc/*", type:"dir" into fastqc_report_ch

  cpus 2

  memory '8 GB'

  publishDir params.output, pattern: "*.html", mode: "copy", overwrite: true

  script:
  fastqc_path = "fastqc/${sample_id}/"
  """
  mkdir -p ${fastqc_path}
  zcat ${fastq1} ${fastq2} | fastqc --quiet -o ${fastqc_path} stdin:${sample_id}
  """
}

process mosdepth {
   label 'mosdepth'
   tag "${sample_id}"

   input:
      file(bed) from bed_file
      tuple val(sample_id), val(bam_type), file(bam), file(bai) from mosdepth_qc_ch
   output:
      file "${sample_id}.${bam_type}.regions.bed.gz"
      file "${sample_id}.${bam_type}.mosdepth.region.dist.txt" into mosdepth_out_ch
 
   memory '4 GB'
 
   cpus 4 // per docs, no benefit after 4 threads
 
   publishDir params.output

   script:
   """
   mosdepth -t ${task.cpus} --by ${bed} --no-per-base --fast-mode ${sample_id}.${bam_type} ${bam}
   """
}


process multiqc {
  label 'multiqc'

  input:
     path('*') from fastqc_report_ch.flatMap().collect()
     path('*') from hs_metrics_out_ch.flatMap().collect()
     path('*') from insert_size_metrics_ch.flatMap().collect()
     path('*') from histogram_ch.flatMap().collect()
     path("*") from mosdepth_out_ch.flatMap().collect()

  output:
     file "multiqc_report.${params.run_id}/multiqc_data.json"
     file "qc_summary.${params.run_id}.mqc.tsv"

  memory '4 GB'
  cpus 4

  publishDir params.output, saveAs: {f -> "multiqc/${f}"}, mode: "copy", overwrite: true

  script:
  """
  multiqc -v -d --filename "multiqc_report_pre.${params.run_id}.html" .
  python preprocess_qc.py picard multiqc_report_pre.${params.run_id}/multiqc_data.json qc_summary.${params.run_id}.mqc.tsv
  multiqc -v -d --filename "multiqc_report.${params.run_id}.html" .
  """
}


/*
   VarDict \
   -G ${reference_fasta} \
   -C  \
   -F 0 \
   -f 0.000000000001 \
   -N ${sample_id} \
   -b ${bam} \
   -c 1 \
   -S 2 \
   -E 3 \
   -g 4 \
   -r 1 \
   -q 1 \
   -VS SILENT \
   -th ${task.cpus} \
   ${bed_file} |
   teststrandbias.R | 
   var2vcf_valid.pl \
   -N ${sample_id} \
   -f 0.000000000001 \
   > ${sample_id}.vardict.vcf
   """
 }
 
*/