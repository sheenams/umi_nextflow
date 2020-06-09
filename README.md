# umi_nextflow

Developmental pipeline for high-sensitivity variant calling using UMI-tagged short
sequence reads.

## Dependencies

  * Nextflow: `wget -qO- https://get.nextflow.io | bash`
  * Docker

## Running the pipeline

The pipeline is initiated from a directory containing paired-end fastqs. These
fastqs should be of the format `<sample_id>.{1,2}.fastq.gz`

`nextflow run main.nf -profile local --input_dir 302R_fastqs/ --run_id 302R 
--downsample_reads 25000000`

Nextflow will create the directory `work` for scratch work and analysis. By
default, the output from the pipeline will be saved to the `output` directory.

If the pipeline were to encounter an error, one that error is addressed, the
pipeline can be resumed from a cached state by appending the `-resume` flag to
the original command. Resuming the above example from a cached state would look
like:

`nextflow run main.nf -profile local --input_dir 302R_fastqs/ --run_id 302R 
--downsample_reads 25000000 -resume`

## Parameters

**Required**
*  `-profile` The config profile for the run environment. Options {local, uw_batch}.
*  `--input_folder` Path to a directory containing paired-end fastqs. Trailing / is required
*  `--run_id` An ID for the sequencing run

**Optional**
*  `--downsample_reads` Number of reads to sample from each fastq. Currently using 25000000 for research runs
*  `--output` Diretory in which to save pipeline output. Default './output'
*  `--save_intermediate_output` Flag to save intermediate files to output directory
*  `-resume` Flag to resume pipeline from a cached state, if possible

