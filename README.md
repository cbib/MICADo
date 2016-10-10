MICADo
======

**Looking for mutations in PacBio cancer data: an alignment-free method**

MICADo is a tool to perform variant calling on targeted (third) next-generation sequencing data. It's algorithm is based on colored de Bruijn graphs. It has been designed for high-precision variant calling for each sample in a cohort. The intuition behind the method is the substraction of systematic biases or errors present in a cohort in order to signle out real mutations. We evaluated its precision on a sample study of TP53 RNA-Seq targeted sequencing data generated by a PacBio sequencer. 

# Getting help

For any information or help running MICADo, you can get in touch with: 
* [Justine Rudewicz](mailto:justinerudewicz[AT]gmail.com)
* [Hayssam Soueidan](mailto:massyah[AT]gmail.com)
* [Macha Nikolski](mailto:macha[AT]labri.fr)

## Version history 

* v.1.0, 2015-12-01, first release, version used in the accompanying paper

## Installation

### System Requirements

MICADo was implemented in Python (python 2.7 ; http://www.python.org/) and tested under Linux and Mac OS environments. 
Python modules required are listed in requierement.txt. Use the following comand line to install them:

```{bash}
sudo pip install -r requirements.txt
```

## Usage

### Data Input

* Fastq file of the sample of interest, targeted capture
* Fastq files of the samples of the cohort, generated using the same library preparation and the same sequencer 
* Reference sequence of the targeted region (FASTA). This can be a multi-fasta file describing e.g. multiple isoforms. 
* (Optionnal) A TSV file describing known SNPs and indels to ignore 

### Minimal example 

```{bash}
python src/principal.py 
	--samplekey 158_1  # Sample label for the results 
	--fastq data/fastq/TP53/C_158_1.fastq,data/fastq/TP53/N_158_1.fastq  # Sample fastq file
	--fasta data/reference/reference_TP53.fasta  # Reference sequence in fasta format 
	--snp data/reference/snp_TP53.tab # snp file 
	--kmer_length 20 # K-mer used for the de Bruijn graph construction 
	--experiment TP53 # Experiment label to build a cohort sequence library for resampling
	--npermutations 100 # Number of resampling to perform for computation of significance scores
```

# LICENSE

    Copyright (c) 2015 Justine Rudewicz (1) (justinerudewicz@gmail.com) 
                Hayssam Soueidan (1) (massyah@gmail.com)
                Richard Iggo (2) (R.Iggo@bordeaux.unicancer.fr)
                Macha Nikolski (1,3) (macha@labri.fr)
    (1) CBiB - Universite Victor Segalen Bordeaux,
    146, rue Leo Saignat, 33076 Bordeaux, France
    (2) Institut Bergonié, unité VINCO INSERM 916, Bordeaux, France
    (3) CNRS / LaBRI, Universite Bordeaux 1, 351 cours de la Liberation,
    33405 Talence Cedex, France 

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
