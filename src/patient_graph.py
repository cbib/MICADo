#!/usr/bin/env python
# coding=utf8
import re
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import networkx as nx
from Bio import pairwise2
import collections
from helpers.logger import init_logger

from alteration import alteration as ALT

# from alteration import alteration as ALT
logger = init_logger("Patient graph")

pattern = re.compile('.+\/(\w+)')


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
	def alteration_list_init(self, G_ref, k, min_support, max_len):
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
		start_g_ref = [key for key, v in G_ref.in_degree().items() if in_degree_g_ref_dict[key] == 0][0]  # only one in TP53
		end_g_ref = [key for key, v in G_ref.out_degree().items() if out_degree_g_ref_dict[key] == 0][0]  # only one in TP53
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
					# Read intersection of all nodes in the alt path for G_sample 
					read_set_pathAlt_G_sample = []
					for node in alternative_path:
						read_set_pathAlt_G_sample.append(set(self.dbg_refrm.node[node]['read_list_n']))
					intersect_allnodes_pathAlt_G_sample = set.intersection(*read_set_pathAlt_G_sample)
					if len(intersect_allnodes_pathAlt_G_sample) == 0:
						continue
					# Reference path choice
					# Replace start/end if it's a tips
					if node_start not in G_ref:
						logger.critical("The node %s (read support : %d) is a tips(start)", node_start,
										len(self.dbg_refrm.node[alternative_path[1]]['read_list_n']))
						node_start = start_g_ref
					if node_end not in G_ref:
						logger.critical("The node %s (read support : %d) is a tips(end)", node_end,
										len(self.dbg_refrm.node[alternative_path[1]]['read_list_n']))
						node_end = end_g_ref
					reference_path_list = []
					reference_path = ""
					for i_path in nx.all_simple_paths(G_ref, node_start, node_end):
						reference_path_list.append(i_path)

					if len(reference_path_list) == 0:
						logger.critical("No reference path between %s and %s", node_start, node_end)
						logger.critical("Alternative path : %s", alternative_path)
						continue

					# if there is multiple references paths, check the largest read intersection or the smallest reference tags
					# if no clear criteria for choice is found we keep the first reference path
					if len(reference_path_list) > 1:
						logger.debug("Trying to identify actual reference")
						reference_path = reference_path_list[0]
						size_biggest_intersection = len(list(set(alternative_path) & set(reference_path)))
						logger.debug("Selected ref path num 0 with size %d", size_biggest_intersection)
						# reference_path = None
						# size_biggest_intersection = 0
						# for i_reference_path in range(0, len(reference_path_list)):
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

					self.alteration_list.append(ALT(reference_path, alternative_path, len(intersect_allnodes_pathRef_G_sample),
													len(intersect_allnodes_pathAlt_G_sample), k,
													max(self.total_coverage_node(node_start),
														self.total_coverage_node(node_end)) * min_support / 100))
				# Replace start/end if it was a tips
				node_end = end_node
				node_start = start_node

	def significant_alteration_list_init(self, p_value_threshold=0.001):
		self.significant_alteration_list = []
		for alteration in self.alteration_list:
			# Pour avoir l'ensemble des paths dans signif alt list
			if alteration.pvalue_ratio <= p_value_threshold:
				# if alteration.pvalue_ratio <= 1:
				self.significant_alteration_list.append(alteration)

	def multiple_alternative_path_filter(self):
		to_remove = []
		node_dict = {"end": collections.defaultdict(list), "start": collections.defaultdict(list)}
		for i_alteration in range(0, len(self.significant_alteration_list)):
			node_start = self.significant_alteration_list[i_alteration].reference_path[0]
			node_end = self.significant_alteration_list[i_alteration].reference_path[
				len(self.significant_alteration_list[i_alteration].reference_path) - 1]
			node_dict["start"][node_start].append(i_alteration)
			node_dict["end"][node_end].append(i_alteration)
		for extremity in node_dict.keys():
			for node, liste_of_alterations in node_dict[extremity].items():
				if len(node_dict[extremity][node]) > 1:
					ratio_max = 0
					for i_alteration in node_dict[extremity][node]:
						if self.significant_alteration_list[i_alteration].ratio_read_count > ratio_max:
							ratio_max = self.significant_alteration_list[i_alteration].ratio_read_count
					for i_alteration in node_dict[extremity][node]:
						if self.significant_alteration_list[i_alteration].ratio_read_count != ratio_max:
							to_remove.append(self.significant_alteration_list[i_alteration])
		for alteration in set(to_remove):
			self.significant_alteration_list.remove(alteration)
