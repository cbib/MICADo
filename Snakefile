
import random

rule STAR_align_synth_1:
    input: "alignments/data/synthetic/C_FOOFOO_111_500_50_1_1-1-1.bam"

rule STAR_align_synth:
    input: star_index="STAR_idx/SAindex",reads="data/synthetic/C_FOOFOO_{seed}_{nreads}_{frac}_{nalt}_1-1-1.fastq"
    params:
        sample_name= "C_FOOFOO_{seed}_{nreads}_{frac}_{nalt}_1-1-1",\
        sorted_bam_prefix = "alignments/data/synthetic/C_FOOFOO_{seed}_{nreads}_{frac}_{nalt}_1-1-1.sorted"
    output: sam="alignments/data/synthetic/C_FOOFOO_{seed}_{nreads}_{frac}_{nalt}_1-1-1.sam",\
            bam="alignments/data/synthetic/C_FOOFOO_{seed}_{nreads}_{frac}_{nalt}_1-1-1.bam", \
            sorted_bam="alignments/data/synthetic/C_FOOFOO_{seed}_{nreads}_{frac}_{nalt}_1-1-1.sorted.bam",\
            indexed_bam="alignments/data/synthetic/C_FOOFOO_{seed}_{nreads}_{frac}_{nalt}_1-1-1.sorted.bam.bai"

    shell:"""
        ./bin/STAR --genomeDir STAR_idx --readFilesIn {input.reads} --outFileNamePrefix alignments/
        mv alignments/Aligned.out.sam {output.sam}
        samtools view -b -S {output.sam} > {output.bam}
        samtools sort {output.bam} {params.sorted_bam_prefix}
        samtools index {output.sorted_bam}
     """

rule specific:
    input:expand("data/synthetic/C_FOOFOO_{seed}_{nreads}_{frac}_{nalt}_1-1-1.combined_alterations.json data/synthetic/C_FOOFOO_{seed}_{nreads}_{frac}_{nalt}_1-1-1_C_model_GMAPno40.sorted.bam".split(),\
                seed=[7004],\
                nreads=500,\
                nalt=[1],\
                alt_type=['1-1-1'],\
                frac=['50'])

rule all:
    input:expand("data/synthetic/C_FOOFOO_{seed}_{nreads}_{frac}_{nalt}_1-1-1.combined_alterations.json",\
                nreads=[150,500,700,1000],\
                nalt=[1,2,3],\
                alt_type=['1-1-1'],\
                frac=['035','040','045',"05","10","50"],\
                seed=random.sample(range(100000),k=20))



rule clean:
    shell:"""
        rm data/synthetic/*.json
        rm data/synthetic/*.fastq
        rm data/synthetic/*.bam
        rm data/synthetic/*.sam
        rm data/synthetic/*.bai
    """

rule run_micado:
    priority :2
    input : fasta_ref="data/reference/reference_TP53.fasta",\
            random_sample="data/synthetic/{sample}.fastq",\
            snp_data="data/reference/snp_TP53.tab"
    params : sample_name= "data/synthetic/{sample}"
    log : "exec_logs/micado_log_{sample}.txt"
    output:
            micado_results=temp("data/synthetic/{sample}.significant_alterations.json"),\

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
    input:input_sam="data/alignments/C_model_GMAPno40_NM_000546.5.sam"
    params:
            sample_name="data/synthetic/{sample}_{seed}_{nreads}_{frac}_{nalt}_{altw}"
    log : "exec_logs/sampler_log_{sample}_{seed}_{nreads}_{frac}_{nalt}_{altw}.txt"
    output:random_alt=temp("data/synthetic/{sample}_{seed,\d+}_{nreads,\d+}_{frac,\d+}_{nalt,\d+}_{altw}.fastq"),
           non_alt=temp("data/synthetic/{sample}_{seed,\d+}_{nreads,\d+}_{frac,\d+}_{nalt,\d+}_{altw}_non_alt.fastq"),\
           sampler_results=temp("data/synthetic/{sample}_{seed,\d+}_{nreads,\d+}_{frac,\d+}_{nalt,\d+}_{altw}.alterations.json")
    shell:"""
        source ~/.virtualenvs/micado/bin/activate
        export PYTHONPATH=`pwd`/src

        # build a sample
        /usr/bin/time -l python src/read_sampler/altered_reads_sampler.py --input_sam {input.input_sam}  \
                    --output_file_prefix "{params.sample_name}" \
                    --n_reads {wildcards.nreads} --fraction_altered 0.{wildcards.frac} --n_alterations {wildcards.nalt} --alt_weight {wildcards.altw} \
                    --seed {wildcards.seed} \
                     --systematic_offset -202 2> {log} >> alterations.txt
                    # --output_lowercase \

    """



rule combine_json:
    priority : 50
    input : micado_results="data/synthetic/{sample}.significant_alterations.json",\
            sampler_results="data/synthetic/{sample}.alterations.json"

    output:combined_json=temp("data/synthetic/{sample}.combined_alterations.temp.json"),
           cleaned_json="data/synthetic/{sample}.combined_alterations.json" # we move them outside
    shell:"""
        # merge known alterations and identified alterations
        source ~/.virtualenvs/micado/bin/activate
        export PYTHONPATH=`pwd`/src
        bin/merge_json_objects.py {input.sampler_results} {input.micado_results} > {output.combined_json}
        jq "." < {output.combined_json} > {output.cleaned_json}

    """

rule align_synthetic_data:
    input: fastq="data/synthetic/{sample}.fastq"
    params: TGTGENOME="NM_000546.5",sample="data/synthetic/{sample}"
    output: "data/synthetic/{sample}_C_model_GMAPno40.sorted.bam"
    shell:"""
        gmap --min-intronlength=15000 -t 32 -D data/gmap_genomes/{params.TGTGENOME} \
             -d {params.TGTGENOME} -f samse --read-group-id=EORTC10994 \
             --read-group-name=GemSim_test --read-group-library=MWG1 \
             --read-group-platform=PACBIO {input.fastq} >  {params.sample}_C_model_GMAPno40.sam

        samtools view -b -S {params.sample}_C_model_GMAPno40.sam > {params.sample}_C_model_GMAPno40.bam
        samtools sort {params.sample}_C_model_GMAPno40.bam {params.sample}_C_model_GMAPno40.sorted
        samtools index {params.sample}_C_model_GMAPno40.sorted.bam
    """
# other variant caller

