#!/usr/bin/env bash
./bin/STAR --genomeDir STAR_idx --readFilesIn ./data/synthetic/C_FOOFOO_6127_700_10_1_1-1-1_non_alt.fastq --outFileNamePrefix alignments/
samtools view -b -S alignments/Aligned.out.sam > alignments/Aligned.out.bam
samtools sort alignments/Aligned.out.bam alignments/Aligned.out.sorted
samtools index alignments/Aligned.out.sorted.bam