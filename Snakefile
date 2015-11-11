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

include: "Snakemake_tools"

#rule gmap_tgt_genome:
#    input : expand("data/gmap_genomes/{tgt}/{tgt}/{tgt}.version",tgt=IDXGENOMENAME)


#rule gen_STAR_align_synth_1:
#    input: bam=XPDIR+"alignments/STAR/C_SYNTHP53_111_1000_05_1_1-1-1_on_NM_000546_5.sorted.bam",ref_fasta=REFFASTA



#rule view_STAR_align_synth_1:
#    input: bam=XPDIR+"alignments/STAR/C_SYNTHP53_111_1000_05_1_1-1-1_on_NM_000546_5.sorted.bam",ref_fasta=REFFASTA
#    shell:"""
#    samtools tview {input.bam} {input.ref_fasta}
#    """

rule view_GMAP_align_synth_1:
    input: bam=XPDIR+"alignments/GMAP/C_SYNTHP53_111_1000_05_1_1-1-1_on_NM_000546_5.sorted.bam",ref_fasta=REFFASTA
    shell:"""
    samtools tview {input.bam} {input.ref_fasta}
    """

TESTSAMPLEPARAMS="C_SYNTHP53_115_1500_10_1_1-1-1"

rule test_varscan:
    input : XPDIR+"results/varscan/"+TESTSAMPLEPARAMS+"_on_NM_000546_5.vcf"
    shell: 'grep -v "^#" {input} | cut -f 2,4,5'

rule test_gatk:
    input : XPDIR+"results/gatk/"+TESTSAMPLEPARAMS+"_on_NM_000546_5_raw.vcf"
    shell: 'grep -v "^#" {input} | cut -f 2,4,5'


rule test_micado:
    input: XPDIR+"results/micado/"+TESTSAMPLEPARAMS+".combined_alterations.json"
