#!/usr/bin/env bash
INPUT=test1
gmap --min-intronlength=15000 -t 32 -D data/gmap_genomes/C_NM_000546.5 -d C_NM_000546.5 -f samse --read-group-id=EORTC10994 --read-group-name=GemSim_test --read-group-library=MWG1 --read-group-platform=PACBIO data/synthetic/${INPUT}.fastq >  data/synthetic/${INPUT}_C_model_GMAPno40.sam
samtools view -b -S data/synthetic/${INPUT}_C_model_GMAPno40.sam > data/synthetic/${INPUT}_C_model_GMAPno40.bam
samtools sort data/synthetic/${INPUT}_C_model_GMAPno40.bam data/synthetic/${INPUT}_C_model_GMAPno40.sorted
samtools index data/synthetic/${INPUT}_C_model_GMAPno40.sorted.bam

gmap --min-intronlength=15000 -t 32 -D data/gmap_genomes/C_NM_000546.5 -d C_NM_000546.5 -f samse --read-group-id=EORTC10994 --read-group-name=GemSim_test --read-group-library=MWG1 --read-group-platform=PACBIO data/synthetic/${INPUT}_non_alt.fastq >  data/synthetic/${INPUT}_non_alt_C_model_GMAPno40.sam
samtools view -b -S data/synthetic/${INPUT}_non_alt_C_model_GMAPno40.sam > data/synthetic/${INPUT}_non_alt_C_model_GMAPno40.bam
samtools sort data/synthetic/${INPUT}_non_alt_C_model_GMAPno40.bam data/synthetic/${INPUT}_non_alt_C_model_GMAPno40.sorted
samtools index data/synthetic/${INPUT}_non_alt_C_model_GMAPno40.sorted.bam