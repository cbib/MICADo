#!/usr/bin/env bash
source ~/.virtualenvs/micado/bin/activate
export PYTHONPATH=`pwd`/src

SAMPLENAME=data/synthetic/C_test2
#SAMPLENAME=data/synthetic/bug_2_very_large_deletion/C_test2

# build a sample
python src/read_sampler/altered_reads_sampler.py --input_sam "data/alignments/C_model_GMAPno40_NM_000546.5.sam" \
            --output_file_prefix "${SAMPLENAME}" \
            --n_reads 500 --fraction_altered 0.1 --n_alterations 5 --alt_weight 0,1,0 \
            --systematic_offset -202
# analyse it
# perform the gmap alignment (just in case )
./bin/gmap_align_synthetic_reads.sh

# run micado
python src/MICADo.py --fastq ${SAMPLENAME}.fastq --experiment TP53 --fasta data/reference/reference_TP53.fasta --samplekey synth2 --snp data/reference/snp_TP53.tab --npermutations 20 --pvalue 0.1 --results "${SAMPLENAME}.significant_alterations.json"

# merge known alterations and identified alterations
bin/merge_json_objects.py ${SAMPLENAME}.alterations.json ${SAMPLENAME}.significant_alterations.json > ${SAMPLENAME}.combined.alterations.json

jq < ${SAMPLENAME}.combined.alterations.json


