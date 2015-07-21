#!/usr/bin/env python
# coding=utf8

# python src/principal.py --samplekey 83_1 --fastq /Users/rudewicz/didac/DiDaC/data/fastq/all_pool_trimmed0.1/C_83_1.fastq,/Users/rudewicz/didac/DiDaC/data/fastq/all_pool_trimmed0.1/N_83_1.fastq --fasta p53var1.fasta --kmer_length 20 --npermutations 10 --min_support_percentage 2

import networkx as nx
from argparse import ArgumentParser
from helpers.helpers import time_iterator, get_or_create_dir
from helpers.logger import init_logger
import sys

logger = init_logger('MICADo')

## imports
logger.info("Will import")
# import reference_graph as RG
# import visualization as VISU
# from individugraph import IndividuGraph as IG
# from randomreadsgraph import RandomReadsGraph as RRG
# import forannotation as ANNO

logger.info("Import finished")

def process_sample(kmer_length, min_support_percentage,  n_permutations, sample_key=None, fastq_files=None, fasta_file=None, destination_directory=".", export_gml=False):

	# split of all fastq files
	fastq_files = fastq_files.split(",")
	print fastq_files
	# g_ref construction
	# logger.info("Will build reference graph with k==%d", kmer_length)
	# g_ref = RG.ref_constructor(kmer_length)

if __name__ == "__main__":
	parser = ArgumentParser()
	parser.add_argument('--kmer_length', help='Size of k-mer words', default=20, type=int, required=False)
	parser.add_argument('--fastq', help='FASTQ files to analyse (sep="," ; with all the path)', required=False, type=str)
	parser.add_argument('--fasta', help='FASTA file of reference sequences (need to be in /data/reference)', required=False, type=str)
	parser.add_argument('--min_support_percentage', help='Minimum of read support percentage for node filter', default=3, type=int, required=False)
	parser.add_argument('--samplekey', help='Unique sample key', default="", type=str, required=True)
	parser.add_argument('--npermutations', help="number of permutations / random samples to perform", default=1000, type=int, required=False)
	parser.add_argument("--destdir", help="Output directory", default="output/gml", type=str, required=False)
	parser.add_argument("--export", help="Whether to export graphs to GML", action='store_true')

	args = parser.parse_args()

	process_sample(kmer_length=args.kmer_length, min_support_percentage=args.min_support_percentage, fastq_files=args.fastq, fasta_file=args.fasta, n_permutations=args.npermutations, sample_key=args.samplekey,
					destination_directory=args.destdir, export_gml=args.export)
