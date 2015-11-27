import random
REFFASTA="data/reference/reference_FLT3.fasta"
INPUTSAM="data/experimental_results/TP53/alignments/C_model_GMAPno40_NM_000546.5.sam"
REFFASTADICT="data/reference/reference_FLT3.dict"
IDXGENOMENAME="reference_FLT3"
XPCODE="FLT3"
MICADO_FLAGS="--disable_cycle_breaking"

SNPDATA="data/reference/snp_FLT3.tab"
SYNTHLABEL="SYNTHPFLT3"


XPDIR="data/flt3_analysis/"
SAMPLEKEY="flt3_hayssam"

VARSCAN="java -jar bin/VarScan.v2.4.0.jar"
GATK="java -jar bin/GenomeAnalysisTK.jar"
PICARD_DICT="java -jar bin/picard-1.140.jar CreateSequenceDictionary"
PICARD_RG="java -jar bin/picard-1.140.jar AddOrReplaceReadGroups"

MICADO_N_PERMUTATIONS=1000
include: "Snakefile_tools"



AFLT3SAMPLE="SRR413284_Normal3"
AFLT3SAMPLE="SRR413284_Normal3_oriented"

rule orient_sample:
    input:fastq=XPDIR+"reads/{sample}.fastq",code="bin/orient_reads_in_forward_direction.py"
    output:XPDIR+"reads/{sample}_oriented.fastq"
    shell:"""
    source ~/.virtualenvs/micado/bin/activate
    export PYTHONPATH=`pwd`/src

    python {input.code} \
             --fastq {input.fastq} \
            --forward_primers TGCTGTGCATACAATTCCCTTGGC,GAGAGGCACTCATGTCAGAACTCA\
            --reverse_primers TCTCTGCTGAAAGGTCGCCTGTTT,AGTCCTCCTCTTCTTCCAGCCTTT \
            --output_suffix oriented
    """

rule samtools_stats :
    input: bam=XPDIR+"alignments/GMAP/"+AFLT3SAMPLE+"_on_"+IDXGENOMENAME+".sorted.bam"
    shell:"""
        samtools view {input.bam}| awk '{{print $2}}' | sort | uniq -c
    """

rule view_GMAP_align_flt3_sample:
    input: bam=XPDIR+"alignments/GMAP/"+AFLT3SAMPLE+"_on_"+IDXGENOMENAME+".sorted.bam",ref_fasta=REFFASTA
    shell:"""
    samtools tview {input.bam} {input.ref_fasta}
    """
rule test_varscan:
    input : XPDIR+"results/varscan/"+AFLT3SAMPLE+"_on_"+IDXGENOMENAME+".vcf"
    shell: 'grep -v "^#" {input} | cut -f 2,4,5'

rule test_gatk:
    input : XPDIR+"results/gatk/"+AFLT3SAMPLE+"_on_"+IDXGENOMENAME+"_raw.vcf"
    shell: 'grep -v "^#" {input} | cut -f 2,4,5'


rule test_micado:
    input: XPDIR+"results/micado/"+AFLT3SAMPLE+".significant_alterations.json"
