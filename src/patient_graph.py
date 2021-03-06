#!/usr/bin/env python
# coding=utf8
import pprint as pp
from itertools import ifilter
from operator import itemgetter
import re
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import Levenshtein
import networkx as nx
from Bio import pairwise2
from scipy.spatial.distance import euclidean, cdist
import numpy as np
import collections
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.neighbors import NearestNeighbors
from helpers.logger import init_logger
from alteration import alteration as ALT

logger = init_logger("Patient graph")

pattern = re.compile('.+\/(\w+)')


def identify_anchor_kmer_in_reference_graph_by_composition(reference_graph, kmer_to_anchor):
	bimer_vectorizer = CountVectorizer(ngram_range=(1, 1), analyzer='char')
	anchor_vec = bimer_vectorizer.fit_transform([kmer_to_anchor]).todense().A
	print anchor_vec
	ref_nodes = reference_graph.nodes()
	ref_vec = bimer_vectorizer.transform(ref_nodes).todense().A
	print ref_vec
	# find closest
	distances = cdist(anchor_vec, ref_vec)
	min_index = np.argmin(distances)
	min_dist = distances[0, min_index]
	other_min_index = [(i, x) for i, x in enumerate(distances[0]) if x == min_dist]
	print [ref_nodes[i] for i, x in other_min_index], min_dist, kmer_to_anchor


def identify_anchor_kmer_in_reference_graph(reference_graph, kmer_to_anchor, leftmost=None, rightmost=None, path_length=None):
	"""

	:type reference_graph: nx.DiGraph
	"""
	toposort = {v: k for k, v in enumerate(nx.topological_sort(reference_graph))}
	# print "Righmost is ",rightmost,toposort[rightmost]
	nodes_to_consider = reference_graph.nodes()
	if rightmost:
		idx = toposort[rightmost]
		nodes_to_consider = ifilter(lambda x: toposort[x] <= idx, nodes_to_consider)
	# print "Max is ", idx
	if leftmost:
		idx = toposort[leftmost]
		nodes_to_consider = ifilter(lambda x: toposort[x] >= idx, nodes_to_consider)
	# print "Min is ", idx
	nodes_to_consider = list(nodes_to_consider)

	node_dists = [(node, Levenshtein.distance(node, kmer_to_anchor), Levenshtein.editops(node, kmer_to_anchor)) for node in
				  nodes_to_consider]
	# print "Will search anchor in ",list(node_dists)
	min_dist = min(node_dists, key=itemgetter(1))[1]
	node_dists = [x for x in node_dists if x[1] == min_dist]
	print "Min possible dist is", min_dist
	if rightmost:
		score_func = lambda x: (x[1] - min_dist) + abs(toposort[x[0]] - (toposort[rightmost] - path_length))
	elif leftmost:
		score_func = lambda x: (x[1] - min_dist) + abs(toposort[x[0]] - (toposort[leftmost] + path_length))
	dist_sorted = sorted(node_dists, key=score_func)
	# identify the rightmost node with minimal distance
	return dist_sorted[0][0]


class PatientGraph:
	def __init__(self, fastq_files, kmer_length):
		self.coverage = {}
		self.coverage['total'] = 0
		self.alteration_list = []
		self.dbg = nx.DiGraph()
		self.dbgclean = None
		self.n_reads = 0
		self.kmer_start_set = set()
		self.kmer_end_set = set()

		for f in fastq_files:
			fastq_id = pattern.search(f).group(1)
			logger.info("Considering file %s for fastq %s", f, fastq_id)
			comp = 0
			for record_s in SeqIO.parse(f, "fastq", generic_dna):
				self.n_reads += 1
				sequence = str(record_s.seq)
				comp += 1
				# For tips search
				self.kmer_start_set.add(sequence[0:kmer_length])
				self.kmer_end_set.add(sequence[len(sequence) - kmer_length:len(sequence)])
				for i2 in range(0, len(sequence) - kmer_length):
					curr_kmer = sequence[(i2):(i2 + kmer_length)]
					next_kmer = sequence[(i2 + 1):(i2 + 1 + kmer_length)]
					if next_kmer not in self.dbg:
						self.dbg.add_node(next_kmer, read_list_n={fastq_id + "_" + str(comp)}, fastq_id={fastq_id})
					if curr_kmer in self.dbg:
						self.dbg.node[curr_kmer]['read_list_n'].add(fastq_id + "_" + str(comp))
						self.dbg.node[curr_kmer]['fastq_id'].add(fastq_id)
						if next_kmer not in self.dbg[curr_kmer]:
							self.dbg.add_edge(curr_kmer, next_kmer)
					else:
						self.dbg.add_node(curr_kmer, read_list_n={fastq_id + "_" + str(comp)}, fastq_id={fastq_id})
						self.dbg.add_edge(curr_kmer, next_kmer)
			self.coverage[fastq_id] = comp
			self.coverage['total'] += comp

	# Compute coverage from a node for the graph .dbg 
	def total_coverage_node(self, node):
		coverage_node = 0
		if node not in self.dbg:
			return 0
		for fastq_id in self.dbg.node[node]['fastq_id']:
			coverage_node += self.coverage[fastq_id]
		return coverage_node

	# Delete nodes of a graph G with count < coverage * min_support %
	def graph_cleaned_init(self, min_support):
		self.dbgclean = self.dbg.copy()
		nodes_count_inf_seuil = []
		for n in self.dbg:
			total_coverage = self.total_coverage_node(n)
			if len(self.dbg.node[n]['read_list_n']) <= total_coverage * min_support / 100:
				nodes_count_inf_seuil.append(n)
		self.dbgclean.remove_nodes_from(nodes_count_inf_seuil)

	# Removes edges in G_sample_test which are present in g_reference
	def graph_rmRefEdges_init(self, G2analyse, g_reference):
		self.dbg_refrm = G2analyse.copy()
		self.dbg_refrm.remove_edges_from(g_reference.edges())

	# Creation of the altertion list 
	def alteration_list_init(self, G_ref, kmer_length, min_support, max_len):
		self.alteration_list = []
		# Only nodes in dbg_refrm & G_ref and with in degree > 0 for end nodes and out degree > 0 for start nodes  
		G_ref_nodes_set = set(G_ref.nodes())
		shared_nodes = list(set(self.dbg_refrm.nodes()) & G_ref_nodes_set)
		out_d = self.dbg_refrm.out_degree()
		in_d = self.dbg_refrm.in_degree()
		shared_nodes_start = [x for x in shared_nodes if out_d[x] > 0]
		shared_nodes_end = [x for x in shared_nodes if in_d[x] > 0]
		# Add tips end & start in shared_nodes_end & start
		out_degree_g_testclean_dict = self.dbgclean.out_degree()
		in_degree_g_testclean_dict = self.dbgclean.in_degree()
		out_degree_g_ref_dict = G_ref.out_degree()
		in_degree_g_ref_dict = G_ref.in_degree()
		end_tips_list = [key for key, v in self.dbgclean.out_degree().items() if
						 out_degree_g_testclean_dict[key] == 0 and key not in G_ref and key in self.kmer_end_set]
		start_tips_list = [key for key, v in self.dbgclean.in_degree().items() if
						   in_degree_g_testclean_dict[key] == 0 and key not in G_ref and key in self.kmer_start_set]
		shared_nodes_start.extend(start_tips_list)
		shared_nodes_end.extend(end_tips_list)
		# Search for alternative paths
		for node_start in shared_nodes_start:
			start_node = node_start
			for node_end in shared_nodes_end:
				end_node = node_end
				for alternative_path in nx.all_simple_paths(self.dbg_refrm, node_start, node_end):
					if len(set(alternative_path) & G_ref_nodes_set) > 2:
						continue
					# Compute coverage of the altenative path
					total_coverage = max([self.total_coverage_node(alt_nodes) for alt_nodes in alternative_path])
					# Read intersection of all nodes in the alt path for G_sample 
					read_set_pathAlt_G_sample = []
					for node in alternative_path:
						read_set_pathAlt_G_sample.append(set(self.dbg_refrm.node[node]['read_list_n']))
					intersect_allnodes_pathAlt_G_sample = set.intersection(*read_set_pathAlt_G_sample)
					if len(intersect_allnodes_pathAlt_G_sample) <= total_coverage * min_support / 100:
						continue
					# Reference path choice
					# Replace start/end if it's a tips
					if node_start not in G_ref:
						logger.critical("The node %s (read support : %d) is a tip (start)", node_start,
										len(self.dbg_refrm.node[alternative_path[1]]['read_list_n']))
						anchor = identify_anchor_kmer_in_reference_graph(G_ref, node_start, rightmost=node_end,
																		 path_length=len(alternative_path))
						logger.critical("Node %s anchored to %s", node_start, anchor)
						node_start = anchor

					if node_end not in G_ref:
						logger.critical("The node %s (read support : %d) is a tip (end)", node_end,
										len(self.dbg_refrm.node[alternative_path[1]]['read_list_n']))
						anchor = identify_anchor_kmer_in_reference_graph(G_ref, node_start, leftmost=node_start,
																		 path_length=len(alternative_path))
						logger.critical("Node %s anchored to %s", node_end, anchor)
						node_end = anchor

					reference_path_list = []
					reference_path = ""
					for i_path in nx.all_simple_paths(G_ref, node_start, node_end):
						reference_path_list.append(i_path)

					if len(reference_path_list) == 0:
						logger.critical("No reference path between %s and %s", node_start, node_end)
						logger.critical("Alternative path : %s", alternative_path)
						continue

					# if there is multiple references paths, check the largest read intersection 
					# if read intersection are equal, the reference path is the one with the smaller delta size accordind to the alternative path
					if len(reference_path_list) > 1:
						logger.debug("Trying to identify actual reference")
						reference_path = reference_path_list[0]
						size_biggest_intersection = len(list(set(alternative_path) & set(reference_path)))
						logger.debug("Selected ref path num 0 with size %d", size_biggest_intersection)
						for i_reference_path in range(1, len(reference_path_list)):
							curr_reference_path = reference_path_list[i_reference_path]
							size_intersection = len(list(set(alternative_path) & set(curr_reference_path)))
							if size_intersection > size_biggest_intersection:
								size_biggest_intersection = size_intersection
								reference_path = curr_reference_path
								logger.debug("Switching to ref path num %d with size %d", i_reference_path, size_biggest_intersection)
							elif size_intersection == size_biggest_intersection:
								size_reference_path = len(reference_path)
								size_curr_reference_path = len(curr_reference_path)
								size_alternative_path = len(alternative_path)
								delta_1 = abs(size_reference_path - size_alternative_path)
								delta_2 = abs(size_curr_reference_path - size_alternative_path)
								if delta_2 < delta_1:
									size_biggest_intersection = size_intersection
									reference_path = curr_reference_path
									logger.debug("Switching to ref path num %d with size %d and deltas: %d--%d ", i_reference_path,
												 size_biggest_intersection, delta_2, delta_1)
						assert reference_path
						assert size_biggest_intersection
					else:
						reference_path = reference_path_list[0]
					# Read intersection of all nodes in the reference path for g_patient 
					condition = 0
					read_set_pathRef_G_sample = []
					for node in reference_path:
						if node not in self.dbg:
							condition = 1
							logger.critical("Identified node %s absent from the input DBG", node)
							intersect_allnodes_pathRef_G_sample = "0"  # Weird smoothing, TODO check with justine if required
							# intersect_allnodes_pathRef_G_sample = []
							break
						read_set_pathRef_G_sample.append(set(self.dbg.node[node]['read_list_n']))
					if condition == 0:
						intersect_allnodes_pathRef_G_sample = set.intersection(*read_set_pathRef_G_sample)
					if abs(len(reference_path) - len(alternative_path)) > max_len:
						logger.critical("Disregarding large alteration %s vs %s", reference_path, alternative_path)
						continue

					reference_sequence = ALT.kmerpathToSeq(reference_path, kmer_length)
					# Decompose path if it is multiple
					for atomic_sequence, atomic_path in decompose_multiple_alterations(reference_path, alternative_path, kmer_length):
						self.alteration_list.append(ALT(reference_path, atomic_path, reference_sequence, atomic_sequence,
														len(intersect_allnodes_pathRef_G_sample),
														len(intersect_allnodes_pathAlt_G_sample), kmer_length,
														max(self.total_coverage_node(node_start),
															self.total_coverage_node(node_end)) * min_support / 100))

				# Replace start/end if it was a tips
				node_end = end_node
				node_start = start_node

	def significant_alteration_list_init(self, p_value_threshold=0.001):
		self.significant_alteration_list = []
		for alteration in self.alteration_list:
			if alteration.pvalue_ratio <= p_value_threshold:
				self.significant_alteration_list.append(alteration)

def decompose_multiple_alterations(reference_path, alternative_path, kmer_length):
	reference_sequence = ALT.kmerpathToSeq(reference_path, kmer_length)
	multi_alternative_sequence = ALT.kmerpathToSeq(alternative_path, kmer_length)

	edit_ops = Levenshtein.editops(reference_sequence, multi_alternative_sequence)
	if len(edit_ops) > 2:
		logger.info("Multiple alt when considering ref %s vs alt %s", reference_sequence, multi_alternative_sequence)
		logger.info("Globally apply %s", edit_ops)
	start, end = 0, 0
	while start < len(edit_ops):
		if edit_ops[start] == 'replace':
			atomic_sequence = Levenshtein.apply_edit([edit_ops[start]], reference_sequence, multi_alternative_sequence)
			# print atomic_sequence
			atomic_path = ALT.kmerize(atomic_sequence, kmer_length)
			start += 1
		else:
			start_e = edit_ops[start]
			end = start + 1
			while (end < len(edit_ops)
				   and edit_ops[end][0] == start_e[0]
				   and (start_e[1] == edit_ops[end][1] or start_e[2] == edit_ops[end][2])):
				end += 1
			edit_op_to_apply = edit_ops[start:end]
			start = end
			logger.info("Will apply %s", edit_op_to_apply)
			atomic_sequence = Levenshtein.apply_edit(edit_op_to_apply, reference_sequence, multi_alternative_sequence)
			atomic_path = ALT.kmerize(atomic_sequence, kmer_length)
		# record each atomic alteration
		logger.info("Adding atomic alteration for ref %s vs alt %s", reference_sequence, atomic_sequence)
		yield atomic_sequence, atomic_path
