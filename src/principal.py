#!/usr/bin/env python
# coding=utf8

# python src/principal.py --samplekey 83_1 --fastq /Users/rudewicz/didac/DiDaC/data/fastq/all_pool_trimmed0.1/C_83_1.fastq,/Users/rudewicz/didac/DiDaC/data/fastq/all_pool_trimmed0.1/N_83_1.fastq --fasta /Users/rudewicz/didac/MICADo/data/reference_TP53test.fasta --kmer_length 20 --npermutations 10 --min_support_percentage 2

import networkx as nx
from argparse import ArgumentParser
from helpers.helpers import time_iterator, get_or_create_dir
from helpers.logger import init_logger
import sys

logger = init_logger('MICADo')

## imports
logger.info("Will import")
import reference_graph as RG
import visualization as VISU
from patient_graph import PatientGraph as PG
# from randomreadsgraph import RandomReadsGraph as RRG
# import forannotation as ANNO

logger.info("Import finished")

def process_sample(kmer_length, min_support_percentage,  n_permutations, sample_key=None, fastq_files=None, fasta_file=None, destination_directory=".", export_gml=False):

	# g_ref construction
	logger.info("Will build reference graph with k==%d and fasta=%s", kmer_length,fasta_file)
	g_reference = RG.ref_constructor(kmer_length,fasta_file)

	# Is there cycles in reference graph?
	if list(nx.simple_cycles(g_reference)):
		if kmer_length > 50:
			logger.info("There are always cycle(s) with k==50...exiting")
			sys.exit(0)
		# Check non depassement valeur limite de k 
		return process_sample(kmer_length=kmer_length+1,sample_key=sample_key,fastq_files=fastq_files,fasta_file=fasta_file, min_support_percentage=min_support_percentage, n_permutations=n_permutations, destination_directory=destination_directory, export_gml=export_gml)

	# g_patient construction
	logger.info("Will build patient graph for %s with k==%d and minimum support (percentage) = %d", fastq_files, kmer_length, min_support_percentage)
	fastq_files = fastq_files.split(",")
	g_patient = PG(fastq_files, kmer_length)
	g_patient.graph_cleaned_init(min_support_percentage) 

	# Is there cycles in patient graph?
	if list(nx.simple_cycles(g_patient.dbgclean)):
		if kmer_length > 50:
			logger.info("There are always cycle(s) with k==50...exiting")
			sys.exit(0)
		# Check non depassement valeur limite de k 
		return process_sample(kmer_length=kmer_length+1,sample_key=sample_key,c_fastq_file=c_fastq_file,n_fastq_file=n_fastq_file, min_support_percentage=min_support_percentage, n_permutations=n_permutations, destination_directory=destination_directory, export_gml=export_gml)

	# Some prints for stats 
	dir_stat = get_or_create_dir("output/statistics") 
	# graph stat
	graph_stat_file = open(dir_stat+"/graph_stat_file"+sample_key+".tsv", 'w')
	graph_stat_file.write(
		"%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d"%(
		kmer_length,
		g_reference.size(),
		sample_key,
		g_patient.coverage['total'],
		g_patient.dbg.size(),
		g_patient.dbgclean.size(),
		g_patient.dbg.in_degree().values().count(0),
		g_patient.dbg.out_degree().values().count(0),
		g_patient.dbgclean.in_degree().values().count(0),
		g_patient.dbgclean.out_degree().values().count(0)
		))
	# kmer stat
	kmer_stat_file = open(dir_stat+"/kmer_stat_file"+sample_key+".tsv", 'w')
	for node_print in g_patient.dbg.nodes():
		fragment_print = ",".join(g_patient.dbg.node[node_print]['fastq_id'])
		reads_print = len(g_patient.dbg.node[node_print]['read_list_n'])
		kmer_stat_file.write(
			"%s\t%s\t%s\t%d\n"%(
			sample_key,
			node_print,
			fragment_print,
			reads_print,
			))

	g_patient.graph_rmRefEdges_init(g_patient.dbgclean, g_reference)  # .dbg_refrm creation



# For visualisation
	graph_name = "G_%s_" % sample_key
	if export_gml:
		logger.info("Will save viz graph for %s with k==%d", sample_key, kmer_length)
		get_or_create_dir(destination_directory)
		# for the refrence graph
		g_reference_merge = VISU.merge_reference_graph(g_reference.copy())
		g_reference_visu = VISU.reference_graph_visualization_formatting(g_reference.copy())
		g_reference_merge_visu = VISU.reference_graph_merged_visualization_formatting(g_reference_merge.copy())
		nx.write_gml(g_reference_visu,destination_directory+"/g_reference_visu"+str(kmer_length)+".gml")
		nx.write_gml(g_reference_merge_visu,destination_directory+"/g_reference_merge_visu"+str(kmer_length)+".gml")
		

if __name__ == "__main__":
	parser = ArgumentParser()
	parser.add_argument('--kmer_length', help='Size of k-mer words', default=20, type=int, required=False)
	parser.add_argument('--fastq', help='FASTQ files to analyse (sep="," ; with all the path)', required=False, type=str)
	parser.add_argument('--fasta', help='FASTA file of reference sequences (with all the path)', required=False, type=str)
	parser.add_argument('--min_support_percentage', help='Minimum of read support percentage for node filter', default=3, type=int, required=False)
	parser.add_argument('--samplekey', help='Unique sample key', default="", type=str, required=True)
	parser.add_argument('--npermutations', help="number of permutations / random samples to perform", default=1000, type=int, required=False)
	parser.add_argument("--destdir", help="Output directory", default="output/gml", type=str, required=False)
	parser.add_argument("--export", help="Whether to export graphs to GML", action='store_true')

	args = parser.parse_args()

	process_sample(kmer_length=args.kmer_length, min_support_percentage=args.min_support_percentage, fastq_files=args.fastq, fasta_file=args.fasta, n_permutations=args.npermutations, sample_key=args.samplekey,
					destination_directory=args.destdir, export_gml=args.export)
