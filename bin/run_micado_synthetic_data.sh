#!/usr/bin/env bash
source ~/.virtualenvs/micado/bin/activate
export PYTHONPATH=`pwd`/src
# build a sample
#python src/read_sampler/altered_reads_sampler.py --input_sam "data/alignments/C_model_GMAPno40_NM_000546.5.sam" \
#            --output_file_prefix "data/synthetic/C_test2" \
#            --n_reads 500 --fraction_altered 0.3 --n_alterations 1 --alt_weight 1,1,1 --seed 1446413523
# analyse it

#python src/principal.py --fastq data/synthetic/C_test2.fastq --experiment TP53 --fasta data/reference/reference_TP53.fasta --samplekey synth2 --snp data/reference/snp_TP53.tab --npermutations 10
python src/principal.py --fastq data/synthetic/bug_0_large_deletion/C_test2.fastq --experiment TP53 --fasta data/reference/reference_TP53.fasta --samplekey synth2 --snp data/reference/snp_TP53.tab --npermutations 5