#!/usr/bin/env nextflow

// Setup the various inputs, defined in nexflow.config
fastq_pair_ch = Channel.fromFilePairs(params.input_folder + '*{1,2}.fastq.gz', flat: true).view()

// Assay specific files
picard_bed_file = Channel.fromPath(params.picard_bed)
bed_file = Channel.fromPath(params.bed)
reference_fasta = Channel.fromPath(params.ref_fasta)
reference_index = Channel.fromPath(params.ref_index)

 // Reference genome is used multiple times
reference_fasta.into { bwa_ref; bwa_realign_ref; picard_ref; qc_ref; filter_con_ref }
reference_index.into { bwa_ref_index; bwa_realign_ref_index; picard_ref_index; qc_ref_index; filter_con_ref_index }


process bwa {
  // Align fastqs
  label 'bwa'
  tag "${sample_id}"
  input:
    file(reference_fasta) from bwa_ref
    file("*") from bwa_ref_index.collect()
    set sample_id, file(fastq1), file(fastq2) from fastq_pair_ch

  output:
    set val(sample_id), file('*.bam') into align_ch
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
     set val(sample_id), file(bam) from align_ch

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

 // ToDo: create alignment image with bwa and picard 

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
     set val(sample_id), file('*.grpumi.histogram') into histogram_ch

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
  //ToDo: implement munge html parser 
 /* 
 umi_html,umi_table = e.Command(
     target=['$pfxout/${specimen}.umi_metrics.html',
             '$pfxout/${specimen}.umi_metrics.csv'],
     source=grp_umi_histo,
     action=('$sing_exec_genome '
             '$IMAGES/python-2.7-new.simg '
             '$GENOMES/src/munge-umi/munge plot_umi '
             '$SOURCE '
             '${TARGETS[0]} '
             '${TARGETS[1]} '))
 e.Depends(umi_html, grp_umi_histo)"""
      output files: grpumi.bam, .grpumi.histogram, umi_metrics.html, logs/fgbio_groupreadsbyumi.log
      input: mapped_bam
 } */

 process fgbio_callconsensus{
   /*  Combined each set of reads to generate consensus reads
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
     set val(sample_id), file('*.consensus.bam') into consensus_bam_ch
  
   memory "32G"
   publishDir params.output, overwrite: true
   script:
   """
   fgbio -Xmx${task.memory.toGiga()}g \
   CallMolecularConsensusReads \
   --input=${bam} \
   --output=${sample_id}.consensus.bam \
   --min-reads=2 \
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
     file(reference_fasta) from filter_con_ref
     file("*") from filter_con_ref_index.collect()
     set val(sample_id), file(bam) from consensus_bam_ch

   output:
     set val(sample_id), file('*.filtered_consensus.bam') into (filter_consensus_bam_ch, consensus_fastq_ch)

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
     file(reference_fasta) from bwa_realign_ref
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
     file(reference_fasta) from picard_ref
     file("*") from picard_ref_index.collect()
    
   output:
     set val(sample_id), file('*.final.bam') into (qc_final_bam, vardict_final_bam_ch)
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

  process quality_metrics {
   /* tools: picard, munge 
    files: #files:.hs_metrics,logs/HsMetrics.log,.Quality_Analysis.txt
    inputs mapped_bam, final_bam
    */
   label 'picard'
   tag "${sample_id}"

   input:
     file(bed_file) from picard_bed_file
     set val(sample_id), file(bam) from qc_final_bam
     file(reference_fasta) from qc_ref
     file("*") from qc_ref_index.collect()

   output:
     set val(sample_id), file('*.hs_metrics') into qc_metrics
  
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

 