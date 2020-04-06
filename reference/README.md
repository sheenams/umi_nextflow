

## UMI/MONC V1.3 references

| file | description |
| ---- | ------------|
| MONC_ctDNA1.3_designed_probe_coords_180314_no_chr.probes.bed | 120bp IDT probes/baits |
| MONC_ctDNA1.3_designed_probe_coords_180314_no_chr.bed | Merged set of probes/baits |
| ctDNA_miniPanel_v1.3.designed_targets.bed | Bed file of designed targets |



## Generating Picard interval lists:

```
docker run -v `pwd`:/data quay.io/biocontainers/picard:2.9.2--py35_1 \
picard -Xmx4g BedToIntervalList \
I=/data/MONC_estimated_targets.nochr.bed \
SD=/data/hs37.fa.dict \
O=/data/MONC_estimated_targets.nochr.interval_list
```

And then save to s3://uwlm-personal/nkrumm/umi/reference/MONC_estimated_targets.nochr.interval_list

