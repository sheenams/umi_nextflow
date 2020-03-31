

## Generating References:

1. Download UCSC-RefSeq exon table ([Link](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=816668497_kLaRN7zpFAU6zyFw55iYbS4r4D4Q&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_refGene&hgta_ctDesc=table+browser+query+on+refGene&hgta_ctVis=pack&hgta_ctUrl=&fbUpBases=200&fbExonBases=10&fbIntronBases=0&fbQual=cds&fbDownBases=200&hgta_doGetBed=get+BED)) and save as: `reference/ucsc_refseq_coding_exons_0bp_padding.bed`

2. Add 10bp extra padding to all exons, sort and merge any overlaps:
```
bedtools slop -i ucsc_refseq_coding_exons_0bp_padding.bed -g hs37.fa.dict -b 10 \
 | bedtools sort -i stdin \
 | bedtools merge -i stdin \
 > ucsc_refseq_coding_exons_10bp_padding.merged.bed

sed 's/^chr//g' ucsc_refseq_coding_exons_10bp_padding.merged.bed > ucsc_refseq_coding_exons_10bp_padding.merged.nochr.bed
```

3. Get intersection of covered regions on MONC assay with expanded exon targets:
```
bedtools intersect \
    -a MONC_ctDNA1.3_designed_probe_coords_180314_no_chr.bed \
    -b ucsc_refseq_coding_exons_10bp_padding.merged.nochr.bed \
> exonic_targets.nochr.bed
```

4. Get all other covered regions which are NOT in the intersection set:
```
bedtools intersect \
    -v \
    -a MONC_ctDNA1.3_designed_probe_coords_180314_no_chr.bed \
    -b ucsc_refseq_coding_exons_10bp_padding.merged.nochr.bed \
    > nonexonic_targets.nochr.bed
```
5. Merge non-exonic and exonic together:

```
cat exonic_targets.nochr.bed nonexonic_targets.nochr.bed \
    | bedtools sort -i stdin \
    | bedtools merge -i stdin \
    > MONC_estimated_targets.nochr.bed


## Generating Picard interval lists:

```
docker run -v `pwd`:/data quay.io/biocontainers/picard:2.9.2--py35_1 \
picard -Xmx4g BedToIntervalList \
I=/data/MONC_estimated_targets.nochr.bed \
SD=/data/hs37.fa.dict \
O=/data/MONC_estimated_targets.nochr.interval_list
```

And then save to s3://uwlm-personal/nkrumm/umi/reference/MONC_estimated_targets.nochr.interval_list

