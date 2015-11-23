import random
REFFASTA="data/reference/NM_000546.5.fasta"
REFFASTADICT="data/reference/NM_000546.5.dict"
IDXGENOMENAME="NM_000546_5"

SNPDATA="data/reference/snp_TP53.tab"
SYNTHLABEL="SYNTHP53"
INPUTSAM="data/experimental_results/TP53/alignments/C_model_GMAPno40_NM_000546.5.sam"

XPDIR="data/tp53_analysis/"
SAMPLEKEY="tp53_hayssam"

VARSCAN="java -jar bin/VarScan.v2.4.0.jar"
GATK="java -jar bin/GenomeAnalysisTK.jar"
PICARD_DICT="java -jar bin/picard-1.140.jar CreateSequenceDictionary"
PICARD_RG="java -jar bin/picard-1.140.jar AddOrReplaceReadGroups"

MICADO_N_PERMUTATIONS=1000
include: "Snakefile_tools"



AP53SAMPLE="N_534_1"
AP53SAMPLE="C_158_1"
AP53SAMPLE="N_158_1"
AP53SAMPLE="N_193_1"
AP53SAMPLE="N_183_1"

rule view_GMAP_align_p53_sample:
    input: bam=XPDIR+"alignments/GMAP/"+AP53SAMPLE+"_on_NM_000546_5.sorted.bam",ref_fasta=REFFASTA
    shell:"""
    samtools tview {input.bam} {input.ref_fasta}
    """
rule test_varscan:
    input : XPDIR+"results/varscan/"+AP53SAMPLE+"_on_NM_000546_5.vcf"
    shell: 'grep -v "^#" {input} | cut -f 2,4,5'

rule test_gatk:
    input : XPDIR+"results/gatk/"+AP53SAMPLE+"_on_NM_000546_5_raw.vcf"
    shell: 'grep -v "^#" {input} | cut -f 2,4,5'


rule test_micado:
    input: XPDIR+"results/micado/"+AP53SAMPLE+".significant_alterations.json"

pool_0_samples= [x.split(",")[0] for x in open("data/experimental_results/TP53/pool_0_groups.tsv","r").readlines()[1:]]

rule test_micado_pool_0:
    input: expand(XPDIR+"results/micado/{sample}.significant_alterations.json",sample=pool_0_samples)


rule test_varscan_pool_0:
    input: expand(XPDIR+"results/varscan/{sample}_on_NM_000546_5.vcf",sample=pool_0_samples)

rule test_gatk_pool_0:
    input: expand(XPDIR+"results/gatk/{sample}_on_NM_000546_5_raw.vcf",sample=pool_0_samples)

rule micado_large_deletion_bug:
    input : XPDIR+"results/micado/N_215_1.significant_alterations.json",XPDIR+"results/micado/C_215_1.significant_alterations.json"

rule test_micadp_on_tips:
    input : XPDIR+"results/micado/N_272_2.significant_alterations.json",\
            XPDIR+"results/micado/C_272_2.significant_alterations.json",\
            XPDIR+"results/micado/N_183_1.significant_alterations.json",\
            XPDIR+"results/micado/C_183_1.significant_alterations.json",\
            XPDIR+"results/micado/N_183_2.significant_alterations.json",\
            XPDIR+"results/micado/C_183_2.significant_alterations.json"

rule test_micadp_on_neg_control:
    input : XPDIR+"results/micado/N_158_1.significant_alterations.json",\
            XPDIR+"results/micado/C_158_1.significant_alterations.json",\
            XPDIR+"results/micado/N_193_1.significant_alterations.json",\
            XPDIR+"results/micado/C_193_1.significant_alterations.json",\
            XPDIR+"results/micado/N_319_1.significant_alterations.json",\
            XPDIR+"results/micado/C_319_1.significant_alterations.json"

rule test_micadp_on_large_deletions:
    input : XPDIR+"results/micado/C_83_1.significant_alterations.json",\
            XPDIR+"results/micado/C_83_2.significant_alterations.json",\
            XPDIR+"results/micado/N_183_1.significant_alterations.json",\
            XPDIR+"results/micado/N_183_2.significant_alterations.json",\
            XPDIR+"results/micado/N_215_1.significant_alterations.json",\
            XPDIR+"results/micado/N_272_1.significant_alterations.json",\
            XPDIR+"results/micado/N_272_2.significant_alterations.json",\
            XPDIR+"results/micado/N_276_1.significant_alterations.json"



some_tp53_samples = random.sample([os.path.splitext(x)[0] for x in os.listdir(XPDIR+"/reads/") if x.endswith(".fastq") and x[0]=="N"],50)
# some_samples=["C_872_1", "C_1375_1", "C_1464_2", "C_1083_1", "C_1063_2", "C_522_1", "C_967_1", "C_579_2", "C_890_1", "C_412_1", "C_38_2", "C_986_2", "C_576_2", "C_49_1", "C_1212_1", "C_72_2", "C_565_1", "C_222_1", "C_1442_1"]
rule test_varscan_all:
    input : expand("{dir}results/varscan/{sample}_on_NM_000546_5.vcf",dir=XPDIR,sample=some_tp53_samples)