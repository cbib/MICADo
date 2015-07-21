#!/usr/bin/env python
# coding=utf8
import re

from Bio import SeqIO
from Bio.Alphabet import generic_dna
import networkx as nx
from Bio import pairwise2
import collections
from helpers.logger import init_logger

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
		self.kmer_start_set = set()
		self.kmer_end_set = set()

		for f in fastq_files:
			fastq_id = pattern.search(f).group(1)
			logger.info("Considering file %s for fastq %s", f, fastq_id)
			comp = 0
			for record_s in SeqIO.parse(f, "fastq", generic_dna):
				sequence = str(record_s.seq)
				comp += 1
				# For tips search
				self.kmer_start_set.add(sequence[0:kmer_length])
				self.kmer_end_set.add(sequence[len(sequence)-kmer_length:len(sequence)])
				for i2 in range(0, len(sequence) - kmer_length):
					curr_kmer = sequence[(i2):(i2 + kmer_length)]
					next_kmer = sequence[(i2 + 1):(i2 + 1 + kmer_length)]
					if next_kmer not in self.dbg:
						self.dbg.add_node(next_kmer, read_list_n=set([fastq_id + "_" + str(comp)]), fastq_id=set([fastq_id]))
					if curr_kmer in self.dbg:
						self.dbg.node[curr_kmer]['read_list_n'].add(fastq_id + "_" + str(comp))
						self.dbg.node[curr_kmer]['fastq_id'].add(fastq_id)
						if next_kmer not in self.dbg[curr_kmer]:
							self.dbg.add_edge(curr_kmer, next_kmer)
					else:
						self.dbg.add_node(curr_kmer, read_list_n=set([fastq_id + "_" + str(comp)]), fastq_id=set([fastq_id]))
						self.dbg.add_edge(curr_kmer, next_kmer)
			self.coverage[fastq_id] = comp
			self.coverage['total'] += comp

	# Compute coverage from a node for the graph .dbg 
	def total_coverage_node(self, node):
		coverage_node = 0
		for fastq_id in self.dbg.node[node]['fastq_id']:
			coverage_node += self.coverage[fastq_id]
		return coverage_node

	## Delete nodes of a graph G with count < coverage * alpha %
	def graph_cleaned_init(self, alpha):
		self.dbgclean = self.dbg.copy()
		nodes_count_inf_seuil = []
		for n in self.dbg:
			total_coverage = self.total_coverage_node(n)
			if len(self.dbg.node[n]['read_list_n']) <= total_coverage * alpha / 100:
				nodes_count_inf_seuil.append(n)
		self.dbgclean.remove_nodes_from(nodes_count_inf_seuil)

	# Removes edges in G_sample_test which are present in g_reference
	def graph_rmRefEdges_init(self, G2analyse, g_reference):
		self.dbg_refrm = G2analyse.copy()
		self.dbg_refrm.remove_edges_from(g_reference.edges())



