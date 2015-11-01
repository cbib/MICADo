#!/usr/bin/env bash
source ~/.virtualenvs/micado/bin/activate
export PYTHONPATH=`pwd`/src

python src/read_sampler/altered_reads_sampler.py --input_sam "data/alignments/C_model_GMAPno40_NM_000546.5.sam" \
            --output_file_prefix "data/synthetic/C_test2" \
            --n_reads 500 --fraction_altered 0.3 --n_alterations 2 --alt_weight 1,1,0