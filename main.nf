#!/usr/bin/env nextflow

// Setup the various inputs, defined in nexflow.config
fastq_pair_ch = Channel.fromFilePairs(params.input_folder + '*{1,2}.fastq.gz', flat: true)
reference_fasta = file("/mnt/disk2/com/Genomes/gatk-bundle/human_g1k_v37.fasta")
reference_index = Channel.fromPath("/mnt/disk2/com/Genomes/gatk-bundle/human_g1k_v37.fasta.{amb,ann,bwt,pac,sa}")

process align {

    input:
    file(reference_fasta) from reference_fasta
    file("*") from reference_index.collect()
    set pair_name, 
          file(fastq1), 
          file(fastq2) from fastq_pair_ch
    output:
    set val(pair_name), file('*.sorted.bam') into mapped_bam_ch
    set val(pair_name), file('*.sorted.bam.bai') into mapped_bai_ch

    shell:
    """
    singularity exec --bind /mnt/disk2/com/:/mnt/disk2/com/:ro --bind $PWD --pwd $PWD /mnt/disk2/com/container-images/bwa-0.7.17.simg bwa mem -R "@RG\\tID:${params.specimen}\\tPL:ILLUMINA\\tPU:NA\\tLB:${params.project}\\tSM:${params.specimen}\\t" -Y -t 8 ${genome_file} ${fq_read1} ${fq_read2} > alignment.bam
    """
}



process fgbio_group_umi {
    // output files: grpumi.bam, .grpumi.histogram, umi_metrics.html, logs/fgbio_groupreadsbyumi.log
    // input: mapped_bam
}

process fgbio_callconsensus{}
    // input: grp_umi_bam 
    // output callconsensus.bam

process remap_consensus{
    // input: callconsensus.bam
    // tools: samtools, bwa mem, picard
    // files: .cns_remap.bam, .cns_final.bam, logs/sampe.log, remapped.log
}

process baserecal {
   // #tools: gatk, samtools
    //#files:.recal_table, bqsr.bam,final.bam, logs/gatk_bqsr.log,logs/gatk_apply_bqsr.log 
    // final_bam=SConscript('bqsr.scons', exports='e cns_final_bam')
}

process quality_metrics {
  // tools: picard, munge 
  // files: #files:.hs_metrics,logs/HsMetrics.log,.Quality_Analysis.txt
  // inputs mapped_bam, final_bam
}

process picard_mark_duplicates {
    // For validation, run mark duplicates on initial mapped_bam and final_bam
}

 process coverage {}
    //#tools: bedtools 
    //#files: MeanCoverageBED.txt
    //base_cov,coverage_metrics=SConscript('bedcov.scons', exports='e assay_items final_bam')
    //e.Depends([base_cov,coverage_metrics], final_bam)
}


