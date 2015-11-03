#!/usr/bin/env bash
TGTGENOME=C_NM_000546.5
TGTGENOME=NM_000546.5


INPUT=test1
INPUT=C_test2
#INPUT=bug_1_three_missed_alterations/C_test2
#INPUT=bug_2_very_large_deletion/C_test2
gmap --min-intronlength=15000 -t 32 -D data/gmap_genomes/${TGTGENOME} -d ${TGTGENOME} -f samse --read-group-id=EORTC10994 --read-group-name=GemSim_test --read-group-library=MWG1 --read-group-platform=PACBIO data/synthetic/${INPUT}.fastq >  data/synthetic/${INPUT}_C_model_GMAPno40.sam
samtools view -b -S data/synthetic/${INPUT}_C_model_GMAPno40.sam > data/synthetic/${INPUT}_C_model_GMAPno40.bam
samtools sort data/synthetic/${INPUT}_C_model_GMAPno40.bam data/synthetic/${INPUT}_C_model_GMAPno40.sorted
samtools index data/synthetic/${INPUT}_C_model_GMAPno40.sorted.bam

gmap --min-intronlength=15000 -t 32 -D data/gmap_genomes/${TGTGENOME} -d ${TGTGENOME} -f samse --read-group-id=EORTC10994 --read-group-name=GemSim_test --read-group-library=MWG1 --read-group-platform=PACBIO data/synthetic/${INPUT}_non_alt.fastq >  data/synthetic/${INPUT}_non_alt_C_model_GMAPno40.sam
samtools view -b -S data/synthetic/${INPUT}_non_alt_C_model_GMAPno40.sam > data/synthetic/${INPUT}_non_alt_C_model_GMAPno40.bam
samtools sort data/synthetic/${INPUT}_non_alt_C_model_GMAPno40.bam data/synthetic/${INPUT}_non_alt_C_model_GMAPno40.sorted
samtools index data/synthetic/${INPUT}_non_alt_C_model_GMAPno40.sorted.bam