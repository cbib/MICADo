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

MICADO_N_PERMUTATIONS=25

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



# Synthetic Experimental results
#rule specific:
#    input:expand("{dir}/results/C_{lbl}_{seed}_{nreads}_{frac}_{nalt}_1-1-1.combined_alterations.json data/synthetic/alignments/C_{lbl}_{seed}_{nreads}_{frac}_{nalt}_1-1-1_C_model_GMAPno40.sorted.bam".split(),\
#                dir=XPDIR, \
#                lbl=SYNTHLABEL,\
#                seed=[7004],\
#                nreads=500,\
#                nalt=[1],\
#                alt_type=['1-1-1'],\
#                frac=['50'])
#
#rule all_synthetic:
#    input:expand("{dir}/results/C_{lbl}_{seed}_{nreads}_{frac}_{nalt}_1-1-1.combined_alterations.json",\
#                dir=XPDIR, \
#                lbl=SYNTHLABEL,\
#                nreads=[150,500,700,1000],\
#                nalt=[1,2,3],\
#                alt_type=['1-1-1'],\
#                frac=['035','040','045',"05","10","50"],\
#                seed=random.sample(range(100000),k=4))
#
#
#
#rule clean:
#    shell:"""
#        rm data/synthetic/resuls/*.json
#        rm data/synthetic/reads/*.fastq
#        rm data/synthetic/alignments/*.bam
#        rm data/synthetic/alignments/*.sam
#        rm data/synthetic/alignments/*.bai
#        rm -rf data/synthetic/results/*
#        rm -fr data/gmap_genomes/*
#        rm -fr data/STAR_genomes/*
#
#
#    """


# generate multi synthetic reads, then call with the three tools
AVAIL_SYNTH_READS_FILE, = glob_wildcards(XPDIR+"reads/{id}.fastq")
print("Found %d read set"%(len(AVAIL_SYNTH_READS_FILE)))

#rule build_multi_sample:
#    input : \
#             expand(XPDIR+"reads/C_SYNTHP53_{seed}_500_05_3_1-1-1.fastq",seed=random.sample(range(50000),k=100)),
#             expand(XPDIR+"reads/C_SYNTHP53_{seed}_500_10_3_1-1-1.fastq",seed=random.sample(range(50000),k=100)),
#             expand(XPDIR+"reads/C_SYNTHP53_{seed}_500_50_3_1-1-1.fastq",seed=random.sample(range(50000),k=100)),


rule generate_synthetic_reads:
    input:expand(XPDIR+"reads/C_SYNTHP53_{seed}_{nreads}_{frac}_{nalt}_1-1-1.fastq",\
                dir=XPDIR, \
                lbl=SYNTHLABEL,\
                seed=random.sample(range(100000),k=20),\
                nreads=[150,500,700,1000],\
                nalt=[1,2,3],\
                alt_type=['1-1-1'],\
                frac=['035','040','045',"05","10","50"],\
                )


rule build_multi_sample:
    input :  expand(XPDIR+"reads/C_SYNTHP53_{seed}_{nreads}_{fraction}_3_1-1-1.fastq",seed=random.sample(range(50000),k=100),nreads=[150,800],fraction=["05","10","50","80"])

rule build_multi_sample_0075:
    input :  expand(XPDIR+"reads/C_SYNTHP53_{seed}_{nreads}_{fraction}_3_1-1-1.fastq",seed=random.sample(range(50000),k=250),nreads=[150,500,800],fraction=["075"])


rule eval_varscan_multi:
    input:expand(XPDIR+"results/varscan/{sample}_on_NM_000546_5.vcf",sample=AVAIL_SYNTH_READS_FILE)
    shell : "wc -l data/synthetic/results/varscan/C_SYNTHP53_*500_10_3_*.vcf | sort -n"


rule eval_gatk_multi:
    input:expand(XPDIR+"results/gatk/{sample}_on_NM_000546_5_raw.vcf",sample=AVAIL_SYNTH_READS_FILE)

rule eval_micado_multi:
    input : expand(XPDIR+"results/micado/{sample}.significant_alterations.json",sample=AVAIL_SYNTH_READS_FILE)
    shell:"""
        cat data/synthetic/results/micado/C_SYNTHP53_*_500_10_3_1-1-1.combined_alterations.json | jq  -c \
        '{{"seed":.sampler.parameters.seed,"n":.significant_alterations|length,"inj":.sampler.injected_alterations|length}}'
    """

rule map_avail_reads_micado:
    input: expand(XPDIR+"results/micado/{sample}.significant_alterations.json",sample=AVAIL_SYNTH_READS_FILE)

rule map_avail_reads_varscan:
    input: expand(XPDIR+"results/varscan/{sample}_on_NM_000546_5.vcf",sample=AVAIL_SYNTH_READS_FILE)


rule map_avail_reads:
#    input: expand(XPDIR+"results/micado/{sample}.combined_alterations.json",sample=AVAIL_SYNTH_READS_FILE)
    input: expand(XPDIR+"results/varscan/{sample}_on_NM_000546_5.vcf",sample=AVAIL_SYNTH_READS_FILE),\
           expand(XPDIR+"results/gatk/{sample}_on_NM_000546_5_raw.vcf",sample=AVAIL_SYNTH_READS_FILE),\
           expand(XPDIR+"results/micado/{sample}.combined_alterations.json",sample=AVAIL_SYNTH_READS_FILE),

rule specific_micado_bug0:
    input: XPDIR+"results/micado/C_SYNTHP53_22640_500_05_3_1-1-1.significant_alterations.json"