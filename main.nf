#!/usr/bin/env nextflow

// Setup the various inputs, defined in nexflow.config
if (params.input_source == "flat_folder") {
  fastq_pair_ch = Channel.fromFilePairs(params.input_folder + '*{1,2}.fastq.gz', flat: true)
                         .view()
                         .into{align_input; fastqc_ch}
} else {
  // eg "s3://uwlm-personal/umi_development/fastq/294R/demux-nf-umi/libraries/**/{1,2}.fastq.gz"
  fastq_pair_ch = Channel.fromPath(params.input_folder + "**/{1,2}.fastq.gz")
                       .map { path ->
                          def (fastq, readgroup, library_type, sample_id, rest) = path.toString().tokenize("/").reverse() 
                          return [sample_id, path]
                        }
                       .groupTuple()
                       .map{ 
                          // ensure two files present for each sample
                          key, files -> if (files.size() != 2) error "Samples must each have exactly two FASTQ files." 
                          return [key, files[0], files[1]]
                        }
                       .into{align_input; fastqc_ch}
}
// initialize optional parameters
params.downsample_reads = null
params.save_intermediate_output = false

// Assay specific files
picard_targets = file(params.picard_targets)
picard_baits = file(params.picard_baits)
bed_targets = file(params.bed_targets)
bed_targets_chr = file(params.bed_targets_chr)

// Reference genome is used multiple times
reference_fasta = file(params.ref_fasta)
reference_index = Channel.fromPath(params.ref_index).into { 
  bwa_ref_index;
  bwa_realign_ref_index;
  picard_ref_index;
  qc_ref_index;
  filter_con_ref_index;
  mpileup_ref_index;
  vardict_ref_index;
  smc2_ref_index
}
reference_fasta_chr = file(params.ref_fasta_chr)
reference_index_chr = Channel.fromPath(params.ref_index_chr).set { 
  umivarcal_ref_index
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
     //"standar" bam is coordinate sorted for use in QC metrics and IGV
   publishDir params.output, mode: 'copy', overwrite: true

   script:
   // bwa mem options:
   // -K seed, -C pass tags from FASTQ -> alignment, -Y recommended by GATK?, -p using paired end input
   // seqtk sample options:
   // -s seed
   // -2 -- use a two-pass approach to reduce memory consumption
   if ("downsample_reads" in params)
     """
     seqtk mergepe \
       <(seqtk sample -2 -s 10000000 ${fastq1} ${params.downsample_reads}) \
       <(seqtk sample -2 -s 10000000 ${fastq2} ${params.downsample_reads}) \
     | bwa mem \
       -R'@RG\\tID:${sample_id}\\tSM:${sample_id}' \
       -K 10000000 \
       -C \
       -Y \
       -t${task.cpus}  \
       ${reference_fasta} \
       -p - \
       2> log.txt \
     | samtools sort -t${task.cpus} -m4G - -o ${sample_id}.standard.bam
     
     samtools index ${sample_id}.standard.bam
     """
   else
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
   label 'sambamba'
   tag "${sample_id}"

   input: 
     tuple val(sample_id), file(bam) from align_ch

   output:
     tuple val(sample_id), file('*.sorted.bam') into set_mate_ch

   script:
   """
   sambamba sort --tmpdir=./ \
   --sort-picard \
   --nthreads ${task.cpus} \
   --memory-limit ${task.memory.toGiga()-1}GB \
   --out=${sample_id}.sorted.bam \
   ${bam}
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

   publishDir path: params.output, mode: 'copy', overwrite: true, enabled: params.save_intermediate_output

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
     tuple val(sample_id), val("grpumi"), file('*.grpumi.bam') into qc_grpumi_bam
     file('*.grpumi.histogram') into histogram_ch

   publishDir path: params.output, mode: 'copy', overwrite: true, enabled: params.save_intermediate_output

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
     tuple val(sample_id), file('*.consensus.bam') into consensus_bam_ch
     tuple val(sample_id), val("consensus"), file('*.consensus.bam') into qc_consensus_bam

   publishDir path: params.output, mode: 'copy', overwrite: true, enabled: params.save_intermediate_output

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
     tuple val(sample_id), file('*.filtered_consensus.bam') into filter_consensus_bam_ch
     tuple val(sample_id), val("filtered_consensus"), file('*.filtered_consensus.bam') into qc_filtered_consensus_bam

   publishDir path: params.output, mode: 'copy', overwrite: true, enabled: params.save_intermediate_output

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
   label "sambamba"
   tag "${sample_id}"

   input: 
    tuple sample_id, path(bam) from filter_consensus_bam_ch

   output:
    tuple sample_id, "${sample_id}.sorted_consensus.bam" into (sorted_consensus_ch, sorted_consensus_realignment_ch)
    
   publishDir path: params.output, mode: 'copy', overwrite: true, enabled: params.save_intermediate_output

   script:
   """
   sambamba sort --tmpdir=./ \
   --sort-picard \
   --nthreads ${task.cpus} \
   --memory-limit ${task.memory.toGiga()-1}GB \
   --out=${sample_id}.sorted_consensus.bam \
   ${bam}
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
     tuple sample_id, path(bam) from sorted_consensus_realignment_ch

   output:
    tuple sample_id, "${sample_id}.realigned.bam" into realign_ch
    tuple val(sample_id), val("realigned"), file('*realigned.bam'), file('*.bai') into qc_sorted_final_bam
    //"realigned" bam is coordinate sorted, for QC and IGV viewing

   script:
   """
   samtools bam2fq -n ${bam} | \
   bwa mem \
   -R "@RG\\tID:${sample_id}\\tSM:${sample_id}" \
   -K 10000000 \
   -p \
   -Y \
   -t ${task.cpus} \
   ${reference_fasta} \
   -p - 2> log.txt \
   | samtools sort -t${task.cpus} -m4G - -o ${sample_id}.realigned.bam
   
   samtools index ${sample_id}.realigned.bam
   """
 }

 process sort_realign_bam {
   //  Sort alignment by query name
   label 'sambamba'
   tag "${sample_id}"

   input: 
    tuple sample_id, path(bam) from realign_ch

   output:
    tuple sample_id, "${sample_id}.sorted.bam" into sorted_realign_consensus_ch

   script:
   """
   sambamba sort --tmpdir=./ \
   --sort-picard \
   --nthreads ${task.cpus} \
   --memory-limit ${task.memory.toGiga()-1}GB \
   --out=${sample_id}.sorted.bam \
   ${bam}
   """
 }

 merge_ch = sorted_realign_consensus_ch.join(sorted_consensus_ch, remainder: true)

 process final_bam {
  //Merge consensus bam (unaligned) with aligned bam, which is queryname sorted
   label 'picard'
   tag "${sample_id}"

   input:
    file(reference_fasta) from reference_fasta
    file("*") from picard_ref_index.collect()
    tuple sample_id, path(sorted_bam), path(sorted_filtered_bam) from merge_ch

   output:
     tuple val(sample_id), val("final"), file('*.final.bam'), file('*.bai') into (qc_final_bam, mpileup_bam, vardict_bam, umivarcal_bam, smc2_bam)
     
   publishDir params.output, mode: 'copy', overwrite: true

   script:
   """
   picard -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./ -Dpicard.useLegacyParser=false\
   MergeBamAlignment \
   -UNMAPPED ${sorted_filtered_bam} \
   -ALIGNED ${sorted_bam} \
   -O ${sample_id}.final.bam \
   -R ${reference_fasta} \
   -VALIDATION_STRINGENCY SILENT \
   -SORT_ORDER coordinate \
   -CREATE_INDEX true

   mv ${sample_id}.final.bai ${sample_id}.final.bam.bai
   """
 }

// ######### QC ######### //

// Combine channels for quality here
// into a single quality_ch for processing below.

qc_final_bam.mix(
  qc_standard_bam,
).into{ hs_metrics_ch; mosdepth_qc_ch; temp_x } 

temp_x.mix(
  qc_grpumi_bam,
  qc_consensus_bam,
  qc_filtered_consensus_bam
).map { [it[0], it[1], it[2]] } // excludes bai file which may not be present in stream
.set { simple_count_qc_ch } // use 'set' here because "into operator should be used to connect two or more target channels "

process simple_quality_metrics {
  label 'sambamba'
  tag "${sample_id}"
  
  input:
    tuple val(sample_id), val(bam_type), file(bam) from simple_count_qc_ch

  output:
    file("${sample_id}.${bam_type}.flagstats.txt") into simple_count_out_ch

  script:
  """
  sambamba flagstat ${bam} \
    2> /dev/null \
    > ${sample_id}.${bam_type}.flagstats.txt
  """

}

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

   publishDir params.output, mode: 'copy', overwrite: true
   
   script:
   """
 
   picard -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./ -Dpicard.useLegacyParser=false \
   CollectHsMetrics \
   -TARGET_INTERVALS=${picard_targets} \
   -BAIT_INTERVALS=${picard_baits} \
   -COVERAGE_CAP=100000 \
   -REFERENCE_SEQUENCE=${reference_fasta} \
   -INPUT=${bam} \
   -OUTPUT=${sample_id}.${bam_type}.hs_metrics 

   picard -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./ -Dpicard.useLegacyParser=false \
   CollectInsertSizeMetrics \
   -INCLUDE_DUPLICATES true \
   -INPUT ${bam} \
   -OUTPUT ${sample_id}.${bam_type}.insert_size_metrics \
   -H ${sample_id}.${bam_type}.insert_size_histogram.pdf
   """
}

process fastqc {
  label 'fastqc'
  tag "${sample_id}"

  input:
    tuple sample_id, file(fastq1), file(fastq2) from fastqc_ch

  output:
    path "fastqc/*", type:"dir" into fastqc_report_ch

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
      file(bed) from bed_targets
      tuple val(sample_id), val(bam_type), file(bam), file(bai) from mosdepth_qc_ch
   output:
      file "${sample_id}.${bam_type}.regions.bed.gz"
      file "${sample_id}.${bam_type}.mosdepth.region.dist.txt" into mosdepth_out_ch
 
   publishDir params.output, mode: 'copy', overwrite: true

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
     path('*') from simple_count_out_ch.flatMap().collect()
     path('*') from insert_size_metrics_ch.flatMap().collect()
     path('*') from histogram_ch.flatMap().collect()
     path("*") from mosdepth_out_ch.flatMap().collect()

  output:
     file "multiqc_report.${params.run_id}.html"
     file "multiqc_report.${params.run_id}_data/multiqc_data.json"
     file "qc_summary.${params.run_id}_mqc.csv"

  publishDir params.output, saveAs: {f -> "multiqc/${f}"}, mode: "copy", overwrite: true

  script:
  """
  preprocess_qc.py counts *.flagstats.txt --output qc_counts.${params.run_id}_mqc.csv
  rm *.flagstats.txt
  multiqc -d --filename "multiqc_report_pre.${params.run_id}.html" .
  preprocess_qc.py summary multiqc_report_pre.${params.run_id}_data/multiqc_data.json qc_summary.${params.run_id}_mqc.csv
  rm -rf multiqc_report_pre.${params.run_id}_data
  multiqc -d -e general_stats --filename "multiqc_report.${params.run_id}.html" .
  """
}

/////////////////////
// Variant Calling //
/////////////////////

process mpileup {
  label 'bwa'

  input:
    file(bed) from bed_targets
    file(reference_fasta) from reference_fasta
    file("*") from mpileup_ref_index.filter{ it.toString() =~ /fai$/ }.collect()
    tuple val(sample_id), val(bam_type), file(bam), file(bai) from mpileup_bam

  output:
    file("${sample_id}.${bam_type}.mpileup")

  publishDir params.output, mode: 'copy', overwrite: true

  script:
  //   -A, --count-orphans     do not discard anomalous read pairs                                                                                                
  //  -B, --no-BAQ            disable BAQ (per-Base Alignment Quality)                                                                                           
  //  -d, --max-depth INT     max per-file depth; avoids excessive memory usage [8000]                                                                           
  //  -E, --redo-BAQ          recalculate BAQ on the fly, ignore existing BQs                                                                                    
  //  -q, --min-MQ INT        skip alignments with mapQ smaller than INT [0]                                                                                     
  // -Q, --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]   
  // -l, --positions FILE    skip unlisted positions (chr pos) or regions (BED)                                                                                 
  """
  samtools mpileup \
  --fasta-ref ${reference_fasta} \
  --max-depth 1000000 \
  --count-orphans \
  --redo-BAQ \
  ${bam}
  > ${sample_id}.${bam_type}.mpileup
  """
}
process vardict {
  label 'vardict'

  input:
    file(bed) from bed_targets
    file(reference_fasta) from reference_fasta
    file("*") from vardict_ref_index.filter{ it.toString() =~ /fai$/ }.collect()
    tuple val(sample_id), val(bam_type), file(bam), file(bai) from vardict_bam

  output:
    file("${sample_id}.${bam_type}.vardict.vcf")

  publishDir params.output, mode: 'copy', overwrite: true

  // -f double The threshold for allele frequency, default: 0.01 or 1%
  // -c INT  The column for chromosome
  // -S INT The column for the region start, e.g. gene start
  // -E INT The column for the region end, e.g. gene end
  // -g INT The column for a gene name, or segment annotation 
  // -F bit The hexical to filter reads. Default: 0x504 (filter unmapped reads, 2nd alignments and duplicates). Use -F 0 to turn it off.
  // -q double The phred score for a base to be considered a good call. Default: 22.5 (for Illumina). For PGM, set it to ~15, as PGM tends to underestimate base quality.
  //-r minimum reads The minimum # of variance reads, default: 2
  // -k 10/1 Indicate whether to perform local realignment. Default: 1 or yes. Set to 0 to disable it.

  script:
  vardict_path = "/usr/local/share/vardict-java-1.7.0-0/bin"
  """
   ${vardict_path}/vardict-java \
   -G ${reference_fasta} \
   -f 0.00001 \
   -N ${sample_id} \
   -b ${bam} \
   -c 1 \
   -S 2 \
   -E 3 \
   -g 4 \
   -F 0 \
   -r 1 \
   -q 1 \
   -VS LENIENT \
   -k 1 \
   ${bed} \
  | ${vardict_path}/teststrandbias.R \
  | ${vardict_path}/var2vcf_valid.pl \
   -N ${sample_id} \
   -E -f 0.00001 \
   > ${sample_id}.${bam_type}.vardict.vcf
   """
 }

// process umivarcal_bam {
//   label 'bwa'

//   input:
//     tuple val(sample_id), val(bam_type), file(bam), file(bai) from umivarcal_bam

//   output:
//     tuple val(sample_id), val(bam_type), file('*.umivarcal.sam') into umivarcal_ready_sam

//   publishDir params.output, mode: 'copy', overwrite: true

//   // Copy the UMI from the RX:Z  into the read name in field 1
//   // Add 'chr' to the header and to each line, required by umi-varcal
//   // May require removing GL chromosomes due to samtools naming issue, still testing that
//   script:                                                                               
//   """
//   samtools view -h ${bam} \
//   | awk -F \$'\t' 'BEGIN {OFS = FS} { if(\$0 ~ "^@") {split(\$2,a,":");gsub(/SN:[^\t]*/, "SN:chr"a[2]); print} else if (\$3 != "*") {;split(\$23,a,":"); gsub(\$1,\$1"_"a[3]); sub("^", "chr",\$3); print} }' > ${sample_id}.${bam_type}.umivarcal.sam
  
//   """
// }

// process umivarcal{
//   label "umivarcal" 

//   input:
//     file(bed) from bed_targets_chr
//     file(reference_fasta_chr) from reference_fasta_chr
//     file("*") from umivarcal_ref_index.filter{ it.toString() =~ /fai$/ }.collect()
//     tuple val(sample_id), val(bam_type), file(sam) from umivarcal_ready_sam

//   output:
//     file("${sample_id}.${bam_type}.umivarcal.vcf")
//     file("${sample_id}.${bam_type}.umivarcal.gvcf")
//     file("${sample_id}.${bam_type}.umivarcal.variants")


//   publishDir params.output, mode: 'copy', overwrite: true

//   //This program is extremely unstable 
//   errorStrategy 'ignore'

//   script:                                                                               
//   """
//   python3 /home/umi-varcal/umi-varcal.py call \
//   --fasta ${reference_fasta_chr} \
//   --bed ${bed} \
//   --input ${sam} \
//   --output ${sample_id}.${bam_type}.umivarcal \
//   --cores ${task.cpus} \
//   --min_base_quality 1 \
//   --min_read_quality 1 \
//   --min_mapping_quality 1 \
//   --min_variant_umi 1 \
//   --gvcf True \
//   --alpha 0.9 
//   """
// }


process smc2{
  label "smc2" 

  input:
    file(bed) from bed_targets
    file(reference_fasta) from reference_fasta
    file("*") from smc2_ref_index.filter{ it.toString() =~ /fai$/ }.collect()
    tuple val(sample_id), val(bam_type), file(bam), file(bai) from smc2_bam

  output:
    file("${sample_id}.${bam_type}.smCounter.cut.vcf")

  publishDir params.output, mode: 'copy', overwrite: true

  //Some input files are in the docker image, packaged with program
  script:  
  """
  python /srv/qgen/code/smc2-af0bd7/run.py \
  --bedTarget ${bed} \
  --bamFile ${bam} \
  --outPrefix ${sample_id}.${bam_type} \
  --nCPU ${task.cpus} \
  --minBQ 10 \
  --minMQ 10 \
  --refGenome ${reference_fasta} \
  --bamType consensus \
  --repBed /srv/qgen/data/annotation/SR_LC_SL.full.bed \
  --srBed /srv/qgen/data/annotation/simpleRepeat.full.bed \
  --umiTag RX \
  --hpLen 8 --mismatchThr 6.0 --primerDist 2 --primerSide 1 --minAltUMI 3 --maxAltAllele 2 
  """
}
