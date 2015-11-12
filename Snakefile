import random


REFFASTA="data/reference/NM_000546.5.fasta"
REFFASTADICT="data/reference/NM_000546.5.dict"
IDXGENOMENAME="NM_000546_5"

#REFFASTA="data/reference/reference_TP53_var_alpha.fasta"
#REFFASTADICT="data/reference/reference_TP53_var_alpha.dict"
#IDXGENOMENAME="reference_TP53_var_alpha"

#REFFASTA="data/reference/reference_TP53.fasta"
#REFFASTADICT="data/reference/reference_TP53.dict"
#IDXGENOMENAME="reference_TP53"

SNPDATA="data/reference/snp_TP53.tab"
SYNTHLABEL="SYNTHP53"
INPUTSAM="data/experimental_results/TP53/alignments/C_model_GMAPno40_NM_000546.5.sam"

XPDIR="data/synthetic/"
SAMPLEKEY="synth2"

#STAR="bin/STAR_2.5.0a"
#STAR="bin/STARlong_2.5.0a" #cf https://groups.google.com/forum/#!topic/rna-star/vcNIAQScQt8
#STAR="bin/STAR" #cf https://groups.google.com/forum/#!topic/rna-star/vcNIAQScQt8
VARSCAN="java -jar bin/VarScan.v2.4.0.jar"
GATK="java -jar bin/GenomeAnalysisTK.jar"
PICARD_DICT="java -jar bin/picard-1.140.jar CreateSequenceDictionary"
PICARD_RG="java -jar bin/picard-1.140.jar AddOrReplaceReadGroups"

include: "Snakefile_tools"


rule view_GMAP_align_synth_1:
    input: bam=XPDIR+"alignments/GMAP/C_SYNTHP53_111_1000_05_1_1-1-1_on_NM_000546_5.sorted.bam",ref_fasta=REFFASTA
    shell:"""
    samtools tview {input.bam} {input.ref_fasta}
    """

TESTSAMPLEPARAMS="C_SYNTHP53_115_150_10_3_1-1-1"

rule test_varscan:
    input : XPDIR+"results/varscan/"+TESTSAMPLEPARAMS+"_on_NM_000546_5.vcf"
    shell: 'grep -v "^#" {input} | cut -f 2,4,5'

rule test_gatk:
    input : XPDIR+"results/gatk/"+TESTSAMPLEPARAMS+"_on_NM_000546_5_raw.vcf"
    shell: 'grep -v "^#" {input} | cut -f 2,4,5'


rule test_micado:
    input: XPDIR+"results/micado/"+TESTSAMPLEPARAMS+".combined_alterations.json"

# generate multi
AVAIL_SYNTH_READS_FILE, = glob_wildcards(XPDIR+"reads/{id}.fastq")

rule build_multi_sample:
    input : \
            expand(XPDIR+"reads/C_SYNTHP53_{seed}_500_05_3_1-1-1.fastq",seed=random.sample(range(50000),k=100))
            # expand(XPDIR+"reads/C_SYNTHP53_{seed}_500_10_3_1-1-1.fastq",seed=random.sample(range(50000),k=100))

rule eval_varscan_multi:
    input:expand(XPDIR+"results/varscan/{sample}_on_NM_000546_5.vcf",sample=AVAIL_SYNTH_READS_FILE)
    shell : "wc -l data/synthetic/results/varscan/C_SYNTHP53_*500_10_3_* | sort -n"

rule eval_gatk_multi:
    input:expand(XPDIR+"results/gatk/{sample}_on_NM_000546_5_raw.vcf",sample=AVAIL_SYNTH_READS_FILE)
    shell : "wc -l data/synthetic/results/gatk/C_SYNTHP53_*500_10_3_* | sort -n"

rule eval_micado_multi:
    input : expand(XPDIR+"results/micado/{sample}.combined_alterations.json",sample=AVAIL_SYNTH_READS_FILE)
    shell:"""
        cat data/synthetic/results/micado/C_SYNTHP53_*_500_10_3_1-1-1.combined_alterations.json | jq  -c \
        '{{"seed":.sampler.parameters.seed,"n":.significant_alterations|length,"inj":.sampler.injected_alterations|length}}'
    """

rule map_avail_reads:
    input: expand(XPDIR+"results/varscan/{sample}_on_NM_000546_5.vcf",sample=AVAIL_SYNTH_READS_FILE),\
           expand(XPDIR+"results/gatk/{sample}_on_NM_000546_5_raw.vcf",sample=AVAIL_SYNTH_READS_FILE),\
           expand(XPDIR+"results/micado/{sample}.combined_alterations.json",sample=AVAIL_SYNTH_READS_FILE),