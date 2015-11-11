
import random

REFFASTA="data/reference/NM_000546.5.fasta"
IDXGENOMENAME="NM_000546_5"

#REFFASTA="data/reference/reference_TP53.fasta"
#IDXGENOMENAME="reference_TP53"

SNPDATA="data/reference/snp_TP53.tab"
SYNTHLABEL="SYNTHP53"
INPUTSAM="data/experimental_results/TP53/alignments/C_model_GMAPno40_NM_000546.5.sam"



#rule gmap_tgt_genome:
#    input : expand("data/gmap_genomes/{tgt}/{tgt}/{tgt}.version",tgt=IDXGENOMENAME)

rule view_STAR_align_synth_1:
    input: bam="data/synthetic/alignments/STAR/C_SYNTHP53_111_500_50_8_1-1-1_on_NM_000546_5.sorted.bam",ref_fasta=REFFASTA
    shell:"""
    samtools tview {input.bam} {input.ref_fasta}
    """

rule view_GMAP_align_synth_1:
    input: bam="data/synthetic/alignments/GMAP/C_SYNTHP53_1112_500_50_8_1-1-1_on_NM_000546_5.sorted.bam",ref_fasta=REFFASTA
    shell:"""
    samtools tview {input.bam} {input.ref_fasta}
    """



# Generic rules

## Alignments (for debugging and other pipeline)

### GMAP

rule gmap_build_genome_index:
    output: gmap_genome="data/gmap_genomes/{IDXGENOMENAME}/{IDXGENOMENAME}.version"
    shell:"""
       gmap_build -d {IDXGENOMENAME} {REFFASTA} -D data/gmap_genomes/
    """

rule gmap_align_synthetic_data:
    input: fastq="data/synthetic/reads/{sample}.fastq",\
            gmap_genome="data/gmap_genomes/{IDXGENOMENAME}/{IDXGENOMENAME}.version"
    params: sample="data/synthetic/reads/{sample}",\
            sorted_bam_prefix = "data/synthetic/alignments/GMAP/{sample}_on_{IDXGENOMENAME}.sorted"

    output: temp_sam = "data/synthetic/alignments/GMAP/{sample}_on_{IDXGENOMENAME}.sam",\
            sorted_bam = "data/synthetic/alignments/GMAP/{sample}_on_{IDXGENOMENAME}.sorted.bam",\
            sorted_bam_index = "data/synthetic/alignments/GMAP/{sample}_on_{IDXGENOMENAME}.sorted.bam.bai",\
            temp_bam = "data/synthetic/alignments/GMAP/{sample}_on_{IDXGENOMENAME}.bam"
    shell:"""
        gmap --min-intronlength=15000 -t 32 -D data/gmap_genomes/{IDXGENOMENAME} \
             -d {IDXGENOMENAME} -f samse --read-group-id=EORTC10994 \
             --read-group-name=GemSim_test --read-group-library=MWG1 \
             --read-group-platform=PACBIO {input.fastq} >  {output.temp_sam}

        samtools view -b -S {output.temp_sam} > {output.temp_bam}
        samtools sort {output.temp_bam} {params.sorted_bam_prefix}
        samtools index {output.sorted_bam}
    """

### STAR

rule STAR_build_genome_index:
    input:
        ref_fasta=REFFASTA
    output:
        star_index_dir='STAR_genomes/{IDXGENOMENAME}',star_index_file='STAR_genomes/{IDXGENOMENAME}/SAindex'
    shell:"""
        ./bin/STAR --runMode genomeGenerate --genomeFastaFiles {input.ref_fasta} --genomeDir {output.star_index_dir}
    """

rule STAR_align_synth:
    input:  star_index_file='STAR_genomes/{IDXGENOMENAME}/SAindex',\
            star_genome_dir='STAR_genomes/{IDXGENOMENAME}',\
            reads="data/synthetic/reads/{sample}.fastq"
    params:
        sample            = "data/synthetic/reads/{sample}",\
        alignment_ouput   = "data/synthetic/alignments/STAR/",\
        sorted_bam_prefix = "data/synthetic/alignments/STAR/{sample}_on_{IDXGENOMENAME}.sorted"
    output: sam="data/synthetic/alignments/STAR/{sample}_on_{IDXGENOMENAME}.sam",\
            bam="data/synthetic/alignments/STAR/{sample}_on_{IDXGENOMENAME}.bam", \
            sorted_bam="data/synthetic/alignments/STAR/{sample}_on_{IDXGENOMENAME}.sorted.bam",\
            indexed_bam="data/synthetic/alignments/STAR/{sample}_on_{IDXGENOMENAME}.sorted.bam.bai"

    shell:"""
        mkdir -p data/synthetic/alignments/STAR/
        ./bin/STAR --genomeDir {input.star_genome_dir} --readFilesIn {input.reads} --outFileNamePrefix {params.alignment_ouput}
        mv {params.alignment_ouput}/Aligned.out.sam {output.sam}
        samtools view -b -S {output.sam} > {output.bam}
        samtools sort {output.bam} {params.sorted_bam_prefix}
        samtools index {output.sorted_bam}
     """


# Synthetic Experimental results
rule specific:
    input:expand("data/synthetic/results/C_{lbl}_{seed}_{nreads}_{frac}_{nalt}_1-1-1.combined_alterations.json data/synthetic/alignments/C_{lbl}_{seed}_{nreads}_{frac}_{nalt}_1-1-1_C_model_GMAPno40.sorted.bam".split(),\
                lbl=SYNTHLABEL,\
                seed=[7004],\
                nreads=500,\
                nalt=[1],\
                alt_type=['1-1-1'],\
                frac=['50'])

rule all_synthetic:
    input:expand("data/synthetic/results/C_{lbl}_{seed}_{nreads}_{frac}_{nalt}_1-1-1.combined_alterations.json",\
                lbl=SYNTHLABEL,\
                nreads=[150,500,700,1000],\
                nalt=[1,2,3],\
                alt_type=['1-1-1'],\
                frac=['035','040','045',"05","10","50"],\
                seed=random.sample(range(100000),k=20))



rule clean:
    shell:"""
        rm data/synthetic/resuls/*.json
        rm data/synthetic/reads/*.fastq
        rm data/synthetic/alignments/*.bam
        rm data/synthetic/alignments/*.sam
        rm data/synthetic/alignments/*.bai
    """


rule run_micado:
    priority :2
    input : fasta_ref=REFFASTA,\
            random_sample="data/synthetic/reads/{sample}.fastq",\
            snp_data=SNPDATA
    params : sample_name= "data/synthetic/reads/{sample}"
    log : "exec_logs/micado_log_{sample}.txt"
    output:
            micado_results=temp("data/synthetic/results/{sample}.significant_alterations.json"),\

    shell:"""
        source ~/.virtualenvs/micado/bin/activate
        export PYTHONPATH=`pwd`/src

        # run micado
        /usr/bin/time -l python src/principal.py --fastq {input.random_sample} --experiment TP53 \
                                --fasta {input.fasta_ref} \
                                --samplekey synth2 \
                                --snp {input.snp_data} \
                                --npermutations 20 --pvalue 0.1 \
                                --results {output.micado_results} 2> {log}


"""


rule generate_sample:
    input:input_sam=INPUTSAM
    params:
            output_reads="data/synthetic/reads/{sample}_{seed}_{nreads}_{frac}_{nalt}_{altw}",\
            output_results="data/synthetic/results/{sample}_{seed}_{nreads}_{frac}_{nalt}_{altw}"
    log : "exec_logs/sampler_log_{sample}_{seed}_{nreads}_{frac}_{nalt}_{altw}.txt"
    output:random_alt=temp("data/synthetic/reads/{sample}_{seed,\d+}_{nreads,\d+}_{frac,\d+}_{nalt,\d+}_{altw}.fastq"),
           non_alt=temp("data/synthetic/reads/{sample}_{seed,\d+}_{nreads,\d+}_{frac,\d+}_{nalt,\d+}_{altw}_non_alt.fastq"),\
           sampler_results=temp("data/synthetic/results/{sample}_{seed,\d+}_{nreads,\d+}_{frac,\d+}_{nalt,\d+}_{altw}.alterations.json")
    shell:"""
        source ~/.virtualenvs/micado/bin/activate
        export PYTHONPATH=`pwd`/src

        # build a sample
        /usr/bin/time -l python src/read_sampler/altered_reads_sampler.py --input_sam {input.input_sam}  \
                    --output_reads_prefix "{params.output_reads}" \
                    --output_result_prefix "{params.output_results}" \
                    --n_reads {wildcards.nreads} --fraction_altered 0.{wildcards.frac} --n_alterations {wildcards.nalt} --alt_weight {wildcards.altw} \
                    --seed {wildcards.seed} \
                     --systematic_offset -202 2> {log} >> alterations.txt
                    # --output_lowercase \

    """



# other variant caller



# helpers

rule combine_json:
    priority : 50
    input : micado_results="data/synthetic/results/{sample}.significant_alterations.json",\
            sampler_results="data/synthetic/results/{sample}.alterations.json"

    output:combined_json=temp("data/synthetic/results/{sample}.combined_alterations.temp.json"),
           cleaned_json="data/synthetic/results/{sample}.combined_alterations.json"
    shell:"""
        # merge known alterations and identified alterations
        source ~/.virtualenvs/micado/bin/activate
        export PYTHONPATH=`pwd`/src
        bin/merge_json_objects.py {input.sampler_results} {input.micado_results} > {output.combined_json}
        jq "." < {output.combined_json} > {output.cleaned_json}

    """
