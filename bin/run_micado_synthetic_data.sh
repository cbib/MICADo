#!/usr/bin/env bash
source ~/.virtualenvs/micado/bin/activate
export PYTHONPATH=`pwd`/src
# build a sample
./bin/generate_synthetic_dataset.sh
# analyse it

python src/principal.py --fastq data/synthetic/C_test2.fastq --experiment TP53 --fasta data/reference/reference_TP53.fasta --samplekey synth2 --snp data/reference/snp_TP53.tab --npermutations 100