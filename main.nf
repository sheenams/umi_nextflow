#!/usr/bin/env nextflow
def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --input sample.tsv -profile singularity --assay MONCv1
""".stripIndent()
}
//ToDo: specify assay information 
params.genome = 'GRCh37'
params.assay = 'MONCv1'

// Read in the genome fasta and index
ref_fasta = file(params.genomes[params.genome].ref_fasta, checkIfExists: true)
ref_fasta_fai = file(params.genomes[params.genome].ref_fasta_fai, checkIfExists: true)

bwa_index = Channel.fromPath(params.genomes[params.genome].bwa_index, checkIfExists: true).collect()
ref_fasta_dict = Channel.fromPath(params.genomes[params.genome].ref_fasta_dict, checkIfExists: true).collect()
picard_intervals = file(params.assays[params.assay].picard_intervals, checkIfExists: true)
bed_file = file(params.assays[params.assay].bed_file, checkIfExists: true)

// tupleup the various inputs, defined in nexflow.config
fastq_pair_ch = Channel.fromFilePairs(params.input + '*{1,2}.fastq.gz', flat: true, checkIfExists: true)
//sample_info = file(params.sample_metadata, checkIfExists: true)

//print(sample_info.text)

process bwa {
   //-Y use soft clipping for supplementary alignment
   //-K process INT input bases in each batch regardless of nThreads (for reproducibility)
   //-C append FASTA/FASTQ comment to SAM output
  label "bwa"

  tag "${sample_id}"

  input:
    path ref_fasta
    path bwa_index
    tuple sample_id, fastq1, fastq2 from fastq_pair_ch    
    // The path qualifier should be preferred over file to handle process input and output files when using Nextflow 19.10.0 or later. 
    //tuple val(sample_id), path(fastq1), path(fastq2) from fastq_pair_ch


  output:
    tuple sample_id, "${sample_id}.sam" into align_ch
    //tuple sample_id, path("${sample_id}.sam") into align_ch
    //tuple val(sample_id), path("${sample_id}.sam") into align_ch
    //tuple val(sample_id), path("*.sam") into align_ch

  publishDir params.output, overwrite: true

  cpus 16
  
  script:
  """ 
  bwa mem \
  -R "@RG\\tID:${sample_id}\\tSM:${sample_id}" \
  -K 10000000 \
  -C \
  -Y \
  -t ${task.cpus} \
  ${ref_fasta} \
  ${fastq1} ${fastq2} > ${sample_id}.sam
  """
} 


 process sort_sam {
  //  Sort alignment by query name
   label "picard"

   tag "${sample_id}"

   input: 
    tuple sample_id, path(sam) from align_ch

   output:
    tuple sample_id, "${sample_id}.sorted.bam" into (temp_qc_initial_bam, mate_ch)
  
   memory "32G"
   
   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
   SortSam \
   I=${sam} \
   O=${sample_id}.sorted.bam \
   SORT_ORDER=queryname
   """
 }

 // ToDo: create alignment image with bwa and picard 

 process fgbio_tuplemateinformation{
  //  Adds and/or fixes mate information on paired-end reads
  // tuples the MQ (mate mapping quality), MC (mate cigar string), 
  // ensures all mate-related flag fields are tuple correctly, 
  // and that the mate reference and mate start position are correct.
   label "fgbio"

   tag "${sample_id}"

   input: 
    tuple sample_id, path(bam) from mate_ch

   output:
    tuple sample_id, "${sample_id}.mateinfo.bam" into (mate_info_bam_ch, temp_qc_mate_bam)

   memory "32G"
   publishDir params.output, overwrite: true

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/fgbio-1.1.0.jar \
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
    tuple sample_id, "${sample_id}.grpumi.bam" into (grp_umi_bam_ch, temp_qc_grp_umi_bam)
    tuple sample_id, "${sample_id}.grpumi.histogram" into histogram_ch

   memory "32G"

   publishDir params.output, overwrite: true

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/fgbio-1.1.0.jar \
   GroupReadsByUmi \
   --input=${bam} \
   --output=${sample_id}.grpumi.bam \
   --family-size-histogram=${sample_id}.grpumi.histogram \
   --strategy=adjacency
   """
 }
  
 process fgbio_callconsensus{
   //  Combined each tuple of reads to generate consensus reads
   //    1. base qualities are adjusted
   //   2. consensus sequence called for all reads with the same UMI, base-by-base.
   //    3. consensus raw base quality is modified by incorporating the probability of an error prior to
   //   calls each end of a pair independently, and does not jointly call bases that overlap within a pair. Insertion or deletion
   //     errors in the reads are not considered in the consensus model.
   label "fgbio"

   tag "${sample_id}"

   input: 
    tuple sample_id, path(bam) from grp_umi_bam_ch

   output:
    tuple sample_id, "${sample_id}.consensus.bam" into (consensus_bam_ch, temp_qc_consensus_bam)
  
   memory "32G"
  
   publishDir params.output, overwrite: true

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/fgbio-1.1.0.jar \
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

   label "fgbio"

   tag "${sample_id}"

   input: 
    path ref_fasta
    path ref_fasta_dict
    tuple sample_id, path(bam) from consensus_bam_ch

   output:
    tuple sample_id, "${sample_id}.filtered_consensus.bam" into (filter_consensus_bam_ch, consensus_fastq_ch,temp_qc_filter_bam)

   memory "32G"
  
   publishDir params.output, overwrite: true

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/fgbio-1.1.0.jar \
   FilterConsensusReads \
   --input=${bam} \
   --output=${sample_id}.filtered_consensus.bam \
   --ref=${ref_fasta} \
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
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
   SortSam \
   I=${bam} \
   O=${sample_id}.sorted_consensus.bam \
   SORT_ORDER=queryname
   """
 }

 process bam_to_fastqs {
  // Create interleaved fastq 
   label "picard"

   tag "${sample_id}"

   input:
    //tuple sample_id, path(bam) from consensus_fastq_ch
    tuple sample_id, path(bam) from sorted_consensus_fastq_ch
   output:
    tuple sample_id, "${sample_id}.fastq" into consensus_fastq
  
   memory "32G"

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
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
    path ref_fasta
    path bwa_index
    tuple sample_id, path(fastq) from consensus_fastq
    
   output:
    tuple sample_id, "${sample_id}.realign.sam" into realign_ch
  
   cpus 8 

   script:
   """
   bwa mem \
   -R "@RG\\tID:${sample_id}\\tSM:${sample_id}" \
   -K 10000000 \
   -p \
   -Y \
   -t ${task.cpus} \
   ${ref_fasta} \
   ${fastq} > ${sample_id}.realign.sam
   """
 }

 process sort_realign_sam {
  //  Sort alignment by query name
   label "picard"

   tag "${sample_id}"

   input: 
    tuple sample_id, path(sam) from realign_ch

   output:
    tuple sample_id, "${sample_id}.sorted.bam" into (sorted_realign_consensus_ch, temp_qc_realign_bam)
  
   memory "32G"
  
   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
   SortSam \
   I=${sam} \
   O=${sample_id}.sorted.bam \
   SORT_ORDER=queryname
   """
 }

 merge_ch = sorted_realign_consensus_ch.join(sorted_consensus_ch, remainder: true)
 process final_bam {
  //Merge consensus bam (unaligned) with aligned bam
   label "picard"
   
   tag "${sample_id}"

   input:
    path ref_fasta
    path ref_fasta_dict
    tuple sample_id, path(sorted_bam), path(sorted_filtered_bam) from merge_ch
        
   output:
    tuple sample_id, "${sample_id}*.final.bam" into (qc_final_bam, final_bam_ch)
  
   publishDir params.output, overwrite: true
   memory "32G"

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
   MergeBamAlignment \
   UNMAPPED=${sorted_filtered_bam} \
   ALIGNED=${sorted_bam} \
   O=${sample_id}.final.bam \
   R=${ref_fasta} \
   CREATE_INDEX=true \
   VALIDATION_STRINGENCY=SILENT \
   SORT_ORDER=coordinate
   """
 }
 
  process hs_metrics {
  // hybrid-selection (HS) metrics
   label "picard"

   tag "${sample_id}"

   input:
    path ref_fasta
    path ref_fasta_fai
    path bwa_index
    
    path picard_intervals
    tuple sample_id, path(bam) from qc_final_bam

   output:
    tuple sample_id, "${sample_id}.final_hs_metrics" into qc_metrics
  
   publishDir params.output, overwrite: true
   
   memory "32G"

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
   CollectHsMetrics \
   BAIT_SET_NAME=${picard_intervals} \
   BAIT_INTERVALS=${picard_intervals} \
   TARGET_INTERVALS=${picard_intervals} \
   REFERENCE_SEQUENCE=${ref_fasta} \
   INPUT=${bam} \
   OUTPUT=${sample_id}.final_hs_metrics
  """
  }

process initial_hs_metrics {
  // hybrid-selection (HS) metrics
   label "picard"

   tag "${sample_id}"

   input:
    path ref_fasta
    path ref_fasta_fai
    path bwa_index
    
    path picard_intervals
    tuple sample_id, path(bam) from temp_qc_initial_bam

   output:
    tuple sample_id, "${sample_id}.initial_hs_metrics" 
  
   publishDir params.output, overwrite: true
   
   memory "32G"

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
   CollectHsMetrics \
   BAIT_SET_NAME=${picard_intervals} \
   BAIT_INTERVALS=${picard_intervals} \
   TARGET_INTERVALS=${picard_intervals} \
   REFERENCE_SEQUENCE=${ref_fasta} \
   INPUT=${bam} \
   OUTPUT=${sample_id}.initial_hs_metrics
  """
  }
  
  process mate_hs_metrics {
  // hybrid-selection (HS) metrics
   label "picard"

   tag "${sample_id}"

   input:
    path ref_fasta
    path ref_fasta_fai
    path bwa_index
    
    path picard_intervals
    tuple sample_id, path(bam) from temp_qc_mate_bam

   output:
    tuple sample_id, "${sample_id}.mate_hs_metrics" 
  
   publishDir params.output, overwrite: true
   
   memory "32G"

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
   CollectHsMetrics \
   BAIT_SET_NAME=${picard_intervals} \
   BAIT_INTERVALS=${picard_intervals} \
   TARGET_INTERVALS=${picard_intervals} \
   REFERENCE_SEQUENCE=${ref_fasta} \
   INPUT=${bam} \
   OUTPUT=${sample_id}.mate_hs_metrics
  """
  }
  
  process grp_umi_hs_metrics {
  // hybrid-selection (HS) metrics
   label "picard"

   tag "${sample_id}"

   input:
    path ref_fasta
    path ref_fasta_fai
    path bwa_index
    
    path picard_intervals
    tuple sample_id, path(bam) from temp_qc_grp_umi_bam

   output:
    tuple sample_id, "${sample_id}.grp_umi_hs_metrics" 
  
   publishDir params.output, overwrite: true
   
   memory "32G"

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
   CollectHsMetrics \
   BAIT_SET_NAME=${picard_intervals} \
   BAIT_INTERVALS=${picard_intervals} \
   TARGET_INTERVALS=${picard_intervals} \
   REFERENCE_SEQUENCE=${ref_fasta} \
   INPUT=${bam} \
   OUTPUT=${sample_id}.grp_umi_hs_metrics
  """
  }
  
  process realign_hs_metrics {
  // hybrid-selection (HS) metrics
   label "picard"

   tag "${sample_id}"

   input:
    path ref_fasta
    path ref_fasta_fai
    path bwa_index
    
    path picard_intervals
    tuple sample_id, path(bam) from temp_qc_realign_bam

   output:
    tuple sample_id, "${sample_id}.realign_hs_metrics" 
  
   publishDir params.output, overwrite: true
   
   memory "32G"

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
   CollectHsMetrics \
   BAIT_SET_NAME=${picard_intervals} \
   BAIT_INTERVALS=${picard_intervals} \
   TARGET_INTERVALS=${picard_intervals} \
   REFERENCE_SEQUENCE=${ref_fasta} \
   INPUT=${bam} \
   OUTPUT=${sample_id}.realign_hs_metrics
  """
  }

  /*

  process filtered_hs_metrics {
  // hybrid-selection (HS) metrics
   label "picard"

   tag "${sample_id}"

   input:
    path ref_fasta
    path ref_fasta_fai
    path bwa_index
    path ref_fasta_dict
    path picard_intervals
    tuple sample_id, path(bam) from temp_qc_filter_bam

   output:
    tuple sample_id, "${sample_id}.filtered_hs_metrics" 
  
   publishDir params.output, overwrite: true
   
   memory "32G"

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
   CollectHsMetrics \
   BAIT_SET_NAME=${picard_intervals} \
   BAIT_INTERVALS=${picard_intervals} \
   TARGET_INTERVALS=${picard_intervals} \
   REFERENCE_SEQUENCE=${ref_fasta} \
   INPUT=${bam} \
   OUTPUT=${sample_id}.filtered_hs_metrics
  """
  }
  /*
 These fail with "Sequence dictionaries are not the same size (0, 84)"
  process consensus_hs_metrics {
  // hybrid-selection (HS) metrics
   label "picard"

   tag "${sample_id}"

   input:
    path ref_fasta
    path ref_fasta_fai
    path bwa_index
    
    path picard_intervals
    tuple sample_id, path(bam) from temp_qc_consensus_bam

   output:
    tuple sample_id, "${sample_id}.consensus_hs_metrics"
   publishDir params.output, overwrite: true
   
   memory "32G"

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
   CollectHsMetrics \
   BAIT_SET_NAME=${picard_intervals} \
   BAIT_INTERVALS=${picard_intervals} \
   TARGET_INTERVALS=${picard_intervals} \
   REFERENCE_SEQUENCE=${ref_fasta} \
   INPUT=${bam} \
   OUTPUT=${sample_id}.consensus_hs_metrics
  """
  }
  process sorted_filtered_hs_metrics {
  // hybrid-selection (HS) metrics
   label "picard"

   tag "${sample_id}"

   input:
    path ref_fasta
    path ref_fasta_fai
    path bwa_index
    
    path picard_intervals
    tuple sample_id, path(bam) from temp_qc_sorted_filter_bam

   output:
    tuple sample_id, "${sample_id}.sorted_filtered_hs_metrics" 
  
   publishDir params.output, overwrite: true
   
   memory "32G"

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
   CollectHsMetrics \
   BAIT_SET_NAME=${picard_intervals} \
   BAIT_INTERVALS=${picard_intervals} \
   TARGET_INTERVALS=${picard_intervals} \
   REFERENCE_SEQUENCE=${ref_fasta} \
   INPUT=${bam} \
   OUTPUT=${sample_id}.sorted_filtered_hs_metrics
  """
  }
  */
  
/*
  process mark_duplicates {

   label "picard"

   tag "${sample_id}"

   input:
    tuple sample_id, path(bam) from dup_final_bam_ch

   output:
    tuple sample_id, "${sample_id}.final_dup_metrics" into dup_metrics
    tuple sample_id, path("${sample_id}.final_rmdup.bam") into rmdup_bam
  
   publishDir params.output, overwrite: true
   
   memory "32G"

   script:
   """
   java -Xmx${task.memory.toGiga()}g -jar /opt/picard.jar \
   MarkDuplicates \
   INPUT=${bam} \
   OUTPUT=${sample_id}.final_rmdup.bam \
   METRICS_FILE=${sample_id}.final_dup_metrics \
   VALIDATION_STRINGENCY=SILENT \
   CREATE_INDEX=true 
  """
  }
 

 /* process vardict {
   
  // -C  #Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
  //  -F 0  #The hexical to filter reads using samtools. Default: 0x500 (filter 2nd alignments and duplicates).  Use -F 0 to turn it off
  //  -f 0.000000000001   #The threshold for allele frequency, default: 0.05 or 5%
  //  -N ${specimen}  #The sample name to be used directly
  //  -b ${SOURCES[0]} 
  //  -c 1  #The column for chromosome
  //  -S 2  #The column for the region start, e.g. gene start
  //  -E 3  #The column for the region end, e.g. gene end
  // -g 4  #The column for gene name, or segment annotation
  //  -r 1  #The minimum # of variant reads, default 2
  // -q 10  #The phred score for a base to be considered a good call.  Default: 25 (for Illumina), tuple to 10 to match
  
   label "vardict"

   tag "${sample_id}"
 
   input:
    tuple sample_id, path(bam) from vardict_final_bam_ch
    path bed_file
    path ref_fasta
    path ref_fasta_dict

   output:
    tuple sample_id, path("*.vardict.vcf") into vardict_vcf_ch

   publishDir params.output, overwrite: true

   memory "32G"

   cpus 8

   script:
   """
   VarDict \
   -G ${ref_fasta} \
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
 
 

//  process picard_mark_duplicates {
//       For validation, run mark duplicates on initial mapped_bam and final_bam
//  }

//   process coverage {
//      #tools: bedtools 
//      #files: MeanCoverageBED.txt
//      base_cov,coverage_metrics=SConscript("bedcov.scons', exports='e assay_items final_bam')
//      e.Depends([base_cov,coverage_metrics], final_bam)
//  } 


//   */
 