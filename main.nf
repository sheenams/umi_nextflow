#!/usr/bin/env nextflow

// Setup the various inputs, defined in nexflow.config
fastq_pair_ch = Channel.fromFilePairs(params.input_folder + '*{1,2}.fastq.gz', flat: true)
                       .take(1)
                       .view()
                       .into{align_input; fastqc_ch}

// Assay specific files
picard_bed_file = file(params.picard_bed)
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
     set sample_id, file(fastq1), file(fastq2) from align_input
 
   output:
     set val(sample_id), file('*.bam'), file('*.bai') into (align_ch, qc_standard_bam)
     file("*.bai")
 
   publishDir params.output, overwrite: true
   
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
   | samtools sort -t@${task.cpus} -m4G - -o ${sample_id}.bam
   
   samtools index ${sample_id}.bam
   """
} 

process sort_sam {
   //  Sort alignment by query name
   label 'picard'
   tag "${sample_id}"

   input: 
     set val(sample_id), file(bam), file(bai) from align_ch

   output:
     set val(sample_id), file('*.sorted.bam') into (temp_qc_initial_bam, set_mate_ch)
   
   memory "32GB"

   script:
   """
   picard -Xmx${task.memory.toGiga()}g \
   SortSam \
   I=${bam} \
   O=${sample_id}.sorted.bam \
   SORT_ORDER=queryname
   """
 }

process fgbio_setmateinformation{
  //  Adds and/or fixes mate information on paired-end reads
   label 'fgbio'
   tag "${sample_id}"

   input: 
     set val(sample_id), file(bam) from set_mate_ch

   output:
     set val(sample_id), file('*.mateinfo.bam') into mate_info_bam_ch

   memory "32GB"

   script:
   """
   fgbio -Xmx${task.memory.toGiga()}g \
   SetMateInformation \
   --input ${bam} \
   --output ${sample_id}.mateinfo.bam
   """
 }

process fgbio_group_umi {
   // Groups reads together that appear to have come from the same original molecule.
   label 'fgbio'
   tag "${sample_id}"

   input: 
     set val(sample_id), file(bam) from mate_info_bam_ch

   output:
     set val(sample_id), file('*.grpumi.bam') into grp_umi_bam_ch
     file('*.grpumi.histogram') into histogram_ch

   memory "32G"
   publishDir params.output, overwrite: true

   script:
   """
   fgbio -Xmx${task.memory.toGiga()}g \
   GroupReadsByUmi \
   --input=${bam} \
   --output=${sample_id}.grpumi.bam \
   --family-size-histogram=${sample_id}.grpumi.histogram \
   --strategy=adjacency
   """
 }

process fgbio_callconsensus{
   /* 
      Combined each set of reads to generate consensus reads
      1. base qualities are adjusted
      2. consensus sequence called for all reads with the same UMI, base-by-base.
      3. consensus raw base quality is modified by incorporating the probability of an error prior to
      calls each end of a pair independently, and does not jointly call bases that overlap within a pair. Insertion or deletion
      errors in the reads are not considered in the consensus model.integrating the unique molecular tags
   */
   label 'fgbio'
   tag "${sample_id}"

   input: 
     set val(sample_id), file(bam) from grp_umi_bam_ch

   output:
     set val(sample_id), file('*.consensus.bam') into (consensus_bam_ch, qc_consensus_bam)
  
   memory "32G"
   publishDir params.output, overwrite: true

   script:
   """
   fgbio -Xmx${task.memory.toGiga()}g \
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
  //  --reverse-per-base-tags

   label 'fgbio'
   tag "${sample_id}"

   input: 
     file(reference_fasta) from reference_fasta
     file("*") from filter_con_ref_index.collect()
     set val(sample_id), file(bam) from consensus_bam_ch

   output:
     set val(sample_id), file('*.filtered_consensus.bam') into (filter_consensus_bam_ch, consensus_fastq_ch, qc_filtered_consensus_bam)

   memory "32G"
   publishDir params.output, overwrite: true

   script:
   """
   fgbio -Xmx${task.memory.toGiga()}g \
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
   label 'picard'
   tag "${sample_id}"

   input: 
     set val(sample_id), file(bam) from filter_consensus_bam_ch

   output:
     set val(sample_id), file('*.sorted_filtered.bam') into sorted_filter_consensus_ch 
  
   publishDir params.output, overwrite: true
   memory "32G"

   script:
   """
   picard -Xmx${task.memory.toGiga()}g \
   SortSam \
   I=${bam} \
   O=${sample_id}.sorted_filtered.bam \
   SORT_ORDER=queryname
   """
 }

process bam_to_fastqs {
   label 'picard'
   tag "${sample_id}"

   input:
     set val(sample_id), file(bam) from consensus_fastq_ch

   output:
     set val(sample_id), file('*.fastq') into consensus_fastq
  
   memory "32G"

   script:
   """
   picard -Xmx${task.memory.toGiga()}g \
   SamToFastq \
   I=${bam} \
   FASTQ=${sample_id}.fastq \
   INTERLEAVE=true
   """
 }

 process realign_consensus {
   label 'bwa'
   tag "${sample_id}"

   input:
     set val(sample_id), file(fastq) from consensus_fastq
     file(reference_fasta) from reference_fasta
     file("*") from bwa_realign_ref_index.collect()

   output:
     set val(sample_id), file('*.realign.sam') into realign_ch
  
   cpus 8 

   script:
   """
   bwa mem \
   -R "@RG\\tID:${sample_id}\\tSM:${sample_id}" \
   -K 10000000 \
   -p \
   -Y \
   -t ${task.cpus} \
   ${reference_fasta} \
   ${fastq} > ${sample_id}.realign.sam
   """
 }

 process sort_realign_sam {
   //  Sort alignment by query name
   label 'picard'
   tag "${sample_id}"

   input: 
     set val(sample_id), file(sam) from realign_ch

   output:
     set val(sample_id), file('*.sorted.bam') into sorted_realign_consensus_ch 
  
   memory "32 GB"
  
   script:
   """
   picard -Xmx${task.memory.toGiga()}g \
   SortSam \
   I=${sam} \
   O=${sample_id}.sorted.bam \
   SORT_ORDER=queryname
   """
 }

 process final_bam {
   //Merge consensus bam (unaligned) with aligned bam
   label 'picard'
   tag "${sample_id}"

   input:
     set val(sample_id), file(consensus_bam) from sorted_filter_consensus_ch
     set val(sample_id), file(sorted_bam) from sorted_realign_consensus_ch
     file(reference_fasta) from reference_fasta
     file("*") from picard_ref_index.collect()
    
   output:
     set val(sample_id), file('*.final.bam'), file('*.bai') into qc_final_bam
     file("*.bai")
  
   publishDir params.output, overwrite: true
   memory "32G"

   script:
   """
   picard -Xmx${task.memory.toGiga()}g \
   MergeBamAlignment \
   UNMAPPED=${consensus_bam} \
   ALIGNED=${sorted_bam} \
   O=${sample_id}.final.bam \
   R=${reference_fasta} \
   CREATE_INDEX=true \
   VALIDATION_STRINGENCY=SILENT \
   SORT_ORDER=coordinate
   """
 }


// ######### QC ######### //

// Combine channels for quality here
// into a single quality_ch for processing below.

//qc_consensus_bam,
//qc_filtered_consensus_bam

qc_final_bam.mix(qc_standard_bam)
            .into{ hs_metrics_ch; mosdepth_qc_ch } 


process quality_metrics {
   label 'picard'
   tag "${sample_id}"

   input:
     file(bed_file) from picard_bed_file
     set val(sample_id), file(bam), file(bai) from hs_metrics_ch
     file(reference_fasta) from reference_fasta
     file("*") from qc_ref_index.collect()

   output:
     path('*.hs_metrics') into qc_metrics
  
   publishDir params.output, overwrite: true
   memory "32G"

   script:
   """
   picard -Xmx${task.memory.toGiga()}g \
   CollectHsMetrics \
   BAIT_SET_NAME=${bed_file} \
   BAIT_INTERVALS=${bed_file} \
   TARGET_INTERVALS=${bed_file} \
   REFERENCE_SEQUENCE=${reference_fasta} \
   INPUT=${bam} \
   OUTPUT=${sample_id}.hs_metrics
   """
}

process fastqc {
  label 'fastqc'
  tag "${sample_id}"
  container 'quay.io/biocontainers/fastqc:0.11.8--1'
  tag "${sample_id}"
  cpus 2
  memory '4 GB'

  publishDir params.output, pattern: "*.html", mode: "copy", overwrite: true

  input:
    set sample_id, file(fastq1), file(fastq2) from fastqc_ch
  output:
    path "fastqc/*", type:"dir" into fastqc_report_ch

  script:
  fastqc_path = "fastqc/${sample_id}/"
  """
  mkdir -p ${fastqc_path}
  zcat ${fastq1} ${fastq2} | fastqc --quiet -o ${fastqc_path} stdin:${sample_id}
  """
}

process mosdepth {
   label 'mosdepth_qc'
   tag "${sample_id}"
   container "quay.io/biocontainers/mosdepth:0.2.9--hbeb723e_0"
   memory '4 GB'
   cpus 4 // per docs, no benefit after 4 threads
 
   input:
      file(bed) from bed_file
      set val(sample_id), file(bam), file(bai) from mosdepth_qc_ch
   output:
      file "${sample_id}.region.bed.gz"
      file "${sample_id}.mosdepth.regions.dist.txt" into mosdepth_out_ch

   publishDir params.output

   script:
   """
   mosdepth -t ${task.cpus} --by ${bed} --no-per-base --fast-mode ${sample_id} ${bam}
   """
}


process multiqc {
  label 'multiqc'
  tag "${sample_id}"
  container 'quay.io/biocontainers/multiqc:1.7--py_4'
  memory '4 GB'
  cpus 4
  input:
     path('fastqc.*') from fastqc_report_ch.flatMap().collect()
     path('*.hs_metrics') from qc_metrics.flatMap().collect()
     path('*.grpumi.histogram') from histogram_ch.flatMap().collect()
     path("*.mosdepth.global.dist.txt") from mosdepth_out_ch.flatMap().collect()

  output:
     file "multiqc_report.${params.run_id}.html"
  publishDir params.output, saveAs: {f -> "multiqc/${f}"}, mode: "copy", overwrite: true
  script:
  """
  multiqc -v -d --filename "multiqc_report.${params.run_id}.html" .
  """
}