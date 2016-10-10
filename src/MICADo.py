#!/usr/bin/env python
# coding=utf8

import json
from Bio import pairwise2

import networkx as nx
from argparse import ArgumentParser
import time
from forannotation import find_edit_operations
from helpers.helpers import time_iterator, get_or_create_dir, get_timestamp, get_git_revision_hash
from helpers.logger import init_logger
import sys

logger = init_logger('MICADo')

## imports
logger.info("Will import")
from reference_graph import ReferenceGraph as RG
from patient_graph import PatientGraph as PG
from randomreadsgraph import RandomReadsGraph as RRG

logger.info("Import finished")


def process_sample(kmer_length, min_support_percentage, n_permutations, p_value_threshold, max_len, sample_key=None, fastq_files=None,
				   fasta_file=None, snp_file=None, experiment_name=None, output_results=None, disable_cycle_breaking=False):
	import seq_lib as seq_lib_module
	seq_lib_module.library_itit(experiment_name)


	# g_reference construction
	logger.info("Will build reference graph with k==%d and fasta=%s & snp=%s", kmer_length, fasta_file, snp_file)
	g_reference = RG(kmer_length, fasta_file, snp_file)

	# Is there cycles in reference graph?
	if list(nx.simple_cycles(g_reference.dbg)):
		if kmer_length >= 70:
			logger.info("There are always cycle(s) with k==70...exiting")
			sys.exit(0)
		# Check non depassement valeur limite de k
		logger.info("[Reference graph] Increasing k to %d to remove cycles", kmer_length+1)
		return process_sample(kmer_length=kmer_length + 1, min_support_percentage=min_support_percentage, n_permutations=n_permutations,
							p_value_threshold=p_value_threshold, max_len=max_len, sample_key=sample_key, fastq_files=fastq_files, 
							fasta_file=fasta_file, snp_file=snp_file, experiment_name=experiment_name, output_results=output_results, 
							disable_cycle_breaking=disable_cycle_breaking)

	# g_patient construction
	logger.info("Will build patient graph for %s with k==%d and minimum support = %dpct", fastq_files, kmer_length, min_support_percentage)
	fastq_files = fastq_files.split(",")
	g_patient = PG(fastq_files, kmer_length)
	logger.info("Before cleaning: %d nodes", len(g_patient.dbg))
	g_patient.graph_cleaned_init(min_support_percentage)
	logger.info("After cleaning: %d nodes", len(g_patient.dbgclean))

	# Is there cycles in patient graph?
	if not disable_cycle_breaking and list(nx.simple_cycles(g_patient.dbgclean)):
		if kmer_length >= 70:
			logger.info("There are still cycle(s) with k==70...exiting")
			sys.exit(0)
		# Check non depassement valeur limite de k
		logger.info("[Sample graph] Increasing k to %d to remove cycles", kmer_length+1)
		return process_sample(kmer_length=kmer_length + 1, min_support_percentage=min_support_percentage, n_permutations=n_permutations, 
							 p_value_threshold=p_value_threshold, max_len=max_len, sample_key=sample_key, fastq_files=",".join(fastq_files), 
							 fasta_file=fasta_file, snp_file=snp_file, experiment_name=experiment_name, output_results=output_results)

	# copy g_patient cleaned and remove reference edges on it (.dbg_refrm creation)
	g_patient.graph_rmRefEdges_init(g_patient.dbgclean, g_reference.dbg)

	# search for alternative paths in dbg_refrm (.alteration_list creation)
	g_patient.alteration_list_init(g_reference.dbg, kmer_length, min_support_percentage, max_len)

	### Permutation test ###
	logger.info("Will create random graphs")
	all_possible_kmers = set()
	for an_alt in g_patient.alteration_list:
		all_possible_kmers.update(an_alt.reference_path)
		all_possible_kmers.update(an_alt.alternative_path)

	for _, _ in time_iterator(range(0, n_permutations), logger, msg_prefix="permuting"):
		g_random = RRG(g_patient.coverage, kmer_length, restrict_to=all_possible_kmers, seq_lib_module=seq_lib_module)
		for i in range(0, len(g_patient.alteration_list)):
			i_alteration = g_patient.alteration_list[i]
			ref_path = i_alteration.reference_path
			alt_path = i_alteration.alternative_path
			g_random_data = g_random.check_path(ref_path,
												alt_path,
												i_alteration.min_coverage)
			i_alteration.random_ratio_list.append(g_random_data[0])
			i_alteration.random_reference_count_list.append(g_random_data[1])
			i_alteration.random_alternative_count_list.append(g_random_data[2])

	logger.info("Will generate p-values for %d possible alterations", len(g_patient.alteration_list))
	for i in range(0, len(g_patient.alteration_list)):
		g_patient.alteration_list[i].pvalue_init()

	g_patient.significant_alteration_list_init(p_value_threshold=p_value_threshold)

	# Annotation
	annotate_and_output_results(g_patient, g_reference, output_results)
	# SNP
	dir_stat = get_or_create_dir("output/snp")
	graph_snp = open(dir_stat + "/snp_" + sample_key + ".tsv", 'w')
	for snp_id in g_reference.snp.keys():
		if g_reference.snp[snp_id][1] in g_patient.dbgclean:
			if g_reference.snp[snp_id][0] in g_patient.dbgclean:
				graph_snp.write("%s\t%s\t%d\t%d\n" % (
					sample_key, snp_id, len(g_patient.dbg.node[g_reference.snp[snp_id][0]]['read_list_n']),
					len(g_patient.dbg.node[g_reference.snp[snp_id][1]]['read_list_n'])))
			else:
				graph_snp.write(
					"%s\t%s\t0\t%d\n" % (sample_key, snp_id, len(g_patient.dbg.node[g_reference.snp[snp_id][1]]['read_list_n'])))


def annotate_and_output_results(g_patient, g_reference, output_results):
	import forannotation as ANNO
	annotated_alterations = ANNO.alteration_list_to_transcrit_mutation(g_patient, g_reference)
	# add experiment arguments
	PROGRAMEND = time.time()
	experiment_description = {}
	this_timestamp = get_timestamp()
	experiment_description['timestamp'] = this_timestamp
	experiment_description['exec_time'] = PROGRAMEND - PROGRAMSTART
	experiment_description['parameters'] = vars(args)
	experiment_description['n_reads'] = g_patient.n_reads
	experiment_description['git_revision_hash'] = get_git_revision_hash()
	# experiment_description['memory_usage'] = process.memory_info().rss

	experiment_description['significant_alterations'] = annotated_alterations
	experiment_description['graphs'] = {
		"coverage_total": g_patient.coverage['total'],
		"before_cleaning": len(g_patient.dbg),
		"after_clearning": len(g_patient.dbgclean)
	}
	experiment_description['all_alterations'] = []
	for x in g_patient.alteration_list:
		alteration_description = x.__dict__
		del alteration_description['reference_path']
		del alteration_description['alternative_path']
		del alteration_description['random_alternative_count_list']
		del alteration_description['random_reference_count_list']
		del alteration_description['random_ratio_list']
		alteration_description['edit_operations'] = find_edit_operations(x.reference_sequence, x.alternative_sequence)
		alteration_description['alignment'] = pairwise2.align.globalms(x.reference_sequence, x.alternative_sequence, 2, -3, -5, -2)[0]
		experiment_description['all_alterations'].append(alteration_description)
	# print json.dumps(experiment_description)
	if output_results:
		with open(output_results, "w") as f:
			json.dump(experiment_description, f)


if __name__ == "__main__":
	global PROGRAMSTART
	PROGRAMSTART = time.time()
	parser = ArgumentParser()
	parser.add_argument('--samplekey', help='Unique sample key', default="", type=str, required=True)
	parser.add_argument('--fastq', help='FASTQ files to analyse (sep="," ; with all the path)', type=str, required=True)
	parser.add_argument('--experiment', help='Experiment name, unique for one study (used for library construction)',
						type=str, required=True)
	parser.add_argument('--fasta', help='FASTA file of reference sequences (with all the path)', type=str, required=True)
	parser.add_argument('--kmer_length', help='Size of k-mer words', default=20, type=int, required=False)
	parser.add_argument('--snp', help='SNP file for reference sequence (with all the path)', type=str, required=False)
	parser.add_argument('--min_support_percentage', help='Minimum of read support percentage for node filter', default=3.0, type=float,
						required=False)
	parser.add_argument('--npermutations', help="number of permutations / random samples to perform", default=1000, type=int,
						required=False)
	parser.add_argument("--max_len", help="Maximum allowed indel length", default=250, type=int, required=False)
	parser.add_argument("--pvalue", help="P value threshold for significance", type=float, default=0.001, required=False)
	parser.add_argument("--results", help="Output (as JSON) results file  ", type=str, default=None, required=False)
	parser.add_argument("--disable_cycle_breaking", help="Do not search for k-mer values yielding a DAG", action="store_true", required=False)

	args = parser.parse_args()

	process_sample(
		kmer_length=args.kmer_length,
		min_support_percentage=args.min_support_percentage,
		fastq_files=args.fastq,
		fasta_file=args.fasta,
		snp_file=args.snp,
		experiment_name=args.experiment,
		n_permutations=args.npermutations,
		sample_key=args.samplekey,
		p_value_threshold=args.pvalue,
		output_results=args.results,
		max_len=args.max_len,
		disable_cycle_breaking=args.disable_cycle_breaking
	)
