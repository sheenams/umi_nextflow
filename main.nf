#!/usr/bin/env nextflow

// Setup the various inputs, defined in nexflow.config
params.input_folder = '/mnt/disk10/users/sheenams/umi_nextflow/test/'
fastq_pair_ch = Channel.fromFilePairs(params.input_folder + '*{1,2}.fastq.gz', flat: true)

reference_fasta = Channel.fromPath("/mnt/disk2/com/Genomes/gatk-bundle/human_g1k_v37.fasta")
reference_index = Channel.fromPath("/mnt/disk2/com/Genomes/gatk-bundle/human_g1k_v37.fasta.{amb,ann,bwt,pac,sa,dict}")

 // Reference genome is used multiple times
reference_fasta.into { bwa_ref; bwa_realign_ref; picard_ref }
reference_index.into { bwa_ref_index; bwa_realign_ref_index; picard_ref_index}

/* process simple_bwa {
  echo true 

  input:
    file(reference_fasta) from reference_fasta
    file("*") from reference_index.collect()

    set sample_id, file(fastq1), file(fastq2) from fastq_pair_ch
  output:
    set val(sample_id), file('*.sam') 
    // OR set val(sample_id), file('*.sorted.bam.bai') into mapped_bai_ch
  cpus 8
    
  """ 
  
  bwa mem -R "@RG\\tID:${sample_id}\\tSM:${sample_id}" -Y -t ${task.cpus} ${reference_fasta} ${fastq1} ${fastq2}  > ${sample_id}.sam
  
  """
} */

process bwa {
  // Align fastqs
  label 'bwa'
  input:
    file(reference_fasta) from bwa_ref
    file("*") from bwa_ref_index.collect()
    set sample_id, file(fastq1), file(fastq2) from fastq_pair_ch

  output:
    set val(sample_id), file('*.sam') into align_ch

  cpus 8
  """ 
  bwa mem \
  -R "@RG\\tID:${sample_id}\\tSM:${sample_id}" \
  -K 10000000 \
  -C \
  -Y \
  -t ${task.cpus} \
  ${reference_fasta} \
  ${fastq1} ${fastq2} > ${sample_id}.sam
  """
} 

process sort_sam {
  // Sort alignment by query name
  label 'picard'
  input: 
    set val(sample_id), file(sam) from align_ch

  output:
    set val(sample_id), file('*.sorted.bam') into set_mate_ch 

  """
  java -Xmx4g -jar /opt/picard.jar \
  SortSam \
  I=${sam} \
  O=${sample_id}.sorted.bam \
  SORT_ORDER=queryname
  """
}

// ToDo: create alignment image with bwa and picard 

process fgbio_setmateinformation{
  // Adds and/or fixes mate information on paired-end reads
  label 'fgbio'
  input: 
    set val(sample_id), file(bam) from set_mate_ch

  output:
    set val(sample_id), file('*.mateinfo.bam') into mate_info_bam_ch
  """
  java -Xmx4g -jar /opt/fgbio-1.1.0.jar \
  SetMateInformation \
  --input ${bam} \
  --output ${sample_id}.mateinfo.bam
  """
}

process fgbio_group_umi {
  // Groups reads together that appear to have come from the same original molecule.
  label 'fgbio'
  input: 
    set val(sample_id), file(bam) from mate_info_bam_ch

  output:
    set val(sample_id), file('*.grpumi.bam') into grp_umi_bam_ch
    set val(sample_id), file('*.grpumi.histogram') into histogram_ch

  """
  java -Xmx4g -jar /opt/fgbio-1.1.0.jar \
  GroupReadsByUmi \
  --input=${bam} \
  --output=${sample_id}.grpumi.bam \
  --family-size-histogram=${sample_id}.grpumi.histogram \
  --strategy=adjacency
  """
}
// ToDo: implement munge html parser 
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
    // output files: grpumi.bam, .grpumi.histogram, umi_metrics.html, logs/fgbio_groupreadsbyumi.log
    // input: mapped_bam
} */

process fgbio_callconsensus{
    //Combined each set of reads to generate consensus reads
    // 1. base qualities are adjusted
    // 2. consensus sequence called for all reads with the same UMI, base-by-base.
    // 3. consensus raw base quality is modified by incorporating the probability of an error prior to
    // calls each end of a pair independently, and does not jointly call bases that overlap within a pair. Insertion or deletion
    // errors in the reads are not considered in the consensus model.integrating the unique molecular tags

  label 'fgbio'
  input: 
    set val(sample_id), file(bam) from grp_umi_bam_ch

  output:
    set val(sample_id), file('*.consensus.bam') into (consensus_bam_ch, consensus_fastq_ch)

  """
  java -Xmx4g -jar /opt/fgbio-1.1.0.jar \
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
  """
}

// Should we filter? http://fulcrumgenomics.github.io/fgbio/tools/latest/FilterConsensusReads.html

process bam_to_fastqs {
  label 'picard'
  input:
    set val(sample_id), file(bam) from consensus_fastq_ch

  output:
    set val(sample_id), file('*.fastq') into consensus_fastq

  """
  java -Xmx4g -jar /opt/picard.jar \
  SamToFastq \
  I=${bam} \
  FASTQ=${sample_id}.fastq \
  INTERLEAVE=true
  """
}

process realign_consensus {
  label 'bwa'
  input:
    set val(sample_id), file(fastq) from consensus_fastq
    file(reference_fasta) from bwa_realign_ref
    file("*") from bwa_realign_ref_index.collect()

  output:
    set val(sample_id), file('*.realign.sam') into realign_ch

  cpus 8 
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
  // Sort alignment by query name
  label 'picard'
  input: 
    set val(sample_id), file(sam) from realign_ch

  output:
    set val(sample_id), file('*.sorted.bam') into sorted_realign_consensus_ch 

  """
  java -Xmx4g -jar /opt/picard.jar \
  SortSam \
  I=${sam} \
  O=${sample_id}.sorted.bam \
  SORT_ORDER=queryname
  """
}
process final_bam {
// Merge consensus bam (unaligned) with aligned bam
  label 'picard'
  input:
    set val(sample_id), file(consensus_bam) from consensus_bam_ch
    set val(sample_id), file(sorted_bam) from sorted_realign_consensus_ch
    file(reference_fasta) from picard_ref
    file("*") from picard_ref_index.collect()
    
  output:
    set val(sample_id), file('*.final.bam') into final_bam

  """
  java -Xmx4g -jar /opt/picard.jar \
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
 /*

  SortSam \
  I=${sam} \
  O=${sample_id}.sorted.bam \
  SORT_ORDER=queryname
process remap_consensus{
    // input: callconsensus.bam
    // tools: samtools, bwa mem, picard
    // files: .cns_remap.bam, .cns_final.bam, logs/sampe.log, remapped.log
      label fgbio
  input: 
    set val(sample_id), file(bam) from mate_info_bam_ch

  output:
    set val(sample_id), file('*.grpumi.bam') into grp_umi_bam_ch
    set val(sample_id), file('*.grpumi.histogram') into histogram_ch

  """
  java -Xmx4g -jar /opt/fgbio-1.1.0.jar GroupReadsByUmi --input=${bam} --output=${sample_id}.grpumi.bam --family-size-histogram=${sample_id}.grpumi.histogram --strategy=adjacency
  """
}

cns_remap_bam, log = e.Command(
    target=['$pfxout/${specimen}.cns_remap.bam',
            '$logs/sampe.log'],
    source=[cns_bam,
            '$REF_FASTA'],
    action=('$sing_exec_genome $IMAGES/samtools-1.9.simg samtools fastq ${SOURCES[0]} | ' #convert bam to fastq for input to bwa mem, which does not accept bam files
            '$sing_exec_genome $IMAGES/bwa-0.7.17.simg bwa mem -K 100000000 -t $ncores -p ${SOURCES[1]} - -R "@RG\\tID:${specimen}\\tPL:ILLUMINA\\tPU:NA\\tLB:${project}\\tSM:${specimen}\\t" -Y 2>> ${TARGETS[1]} | '
            '$sing_exec $IMAGES/samtools-1.9.simg samtools view -Shb - | ' #input Sam, print header for SAM, output BAM sorted by query name
            '$sing_exec $IMAGES/samtools-1.9.simg samtools sort -n - -o ${TARGETS[0]} '))
e.Depends(cns_remap_bam, cns_bam)

cns_final_bam, log = e.Command(
    target=['$pfxout/${specimen}.cns_final.bam',
            '$logs/remapped.log'],
    source=['$REF_FASTA',
            cns_bam,
            cns_remap_bam],
    action=('$sing_exec_genome '
            '$IMAGES/picard-2.18.16.simg '
            'java -Xmx4g -jar /opt/picard.jar '
            'MergeBamAlignment '
            'UNMAPPED=${SOURCES[1]} '
            'ALIGNED=${SOURCES[2]} '
            'O=${TARGETS[0]} '
            'R=${SOURCES[0]} '
            'CREATE_INDEX=true '
            'VALIDATION_STRINGENCY=SILENT '
            'SORT_ORDER=coordinate '
            '2> ${TARGETS[1]} '
            '&& mv $pfxout/${specimen}.cns_final.bai $pfxout/${specimen}.cns_final.bam.bai '))
e.Depends(cns_final_bam, [cns_bam, cns_remap_bam])



process quality_metrics {
  // tools: picard, munge 
  // files: #files:.hs_metrics,logs/HsMetrics.log,.Quality_Analysis.txt
  // inputs mapped_bam, final_bam
}

process picard_mark_duplicates {
    // For validation, run mark duplicates on initial mapped_bam and final_bam
}

 process coverage {
    //#tools: bedtools 
    //#files: MeanCoverageBED.txt
    //base_cov,coverage_metrics=SConscript('bedcov.scons', exports='e assay_items final_bam')
    //e.Depends([base_cov,coverage_metrics], final_bam)
} 


 */
