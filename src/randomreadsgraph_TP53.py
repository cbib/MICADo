#!/usr/bin/env python
# coding=utf8
import collections
import networkx as nx
from helpers.logger import init_logger
import seq_lib_TP53 as SL

logger = init_logger("RANDGRAPH")

cached_kmers = {}
cached_pairs = {}
KMER_UID = {}
curr_uid = 0
last_sample = None


class RandomReadsGraph:
	def __init__(self, coverage_dict, k, restrict_to={}):
		global cached_kmers, curr_uid, last_sample
		self.coverage = coverage_dict
		self.kmer_map = collections.defaultdict(set)
		self.restrict_to = set(restrict_to)
		self.possible_pairs = set()
		read_list = SL.sampling(self.coverage)

		self.dbg = nx.DiGraph()
		# logger.info("Will process %d reads",len(read_list))
		for i_read in range(0, len(read_list)):
			# for i_read in range(0, len(read_list[0:1000])):

			this_read = read_list[i_read]
			# if this_read not in cached_kmers:
			# 	these_kmers=[this_read[i:i+k] for i in xrange(len(this_read)-k)]
			# 	for km in these_kmers:
			# 		if km not in KMER_UID:
			# 			KMER_UID[km]=curr_uid
			# 			curr_uid+=1
			# 	cached_kmers[this_read]=map(lambda x:KMER_UID[x],these_kmers)
			# 	cached_pairs[this_read]=[(f1,f2) for f1,f2 in zip(cached_kmers[this_read],cached_kmers[this_read][1:])]
			#
			# kkmers=cached_kmers[this_read]
			# kmers_pairs=cached_pairs[this_read]

			kkmers = [this_read[i:i + k] for i in xrange(len(this_read) - k) if this_read[i:i + k] in self.restrict_to]
			kmers_pairs = [(f1, f2) for f1, f2 in zip(kkmers, kkmers[1:])]

			for kmer in kkmers:
				self.kmer_map[kmer].add(i_read)

			self.possible_pairs.update(kmers_pairs)

		# self.dbg.add_path(kkmers)
		# for kmer in kkmers:
		# 	if 'read_list_n' not in self.dbg.node[kmer]:
		# 		self.dbg.node[kmer]['read_list_n']=set()
		# 	self.dbg.node[kmer]['read_list_n'].add(i_read)

		# for i_kmer in range(0, len(read_list[i_read]) - k):
		#
		# 	curr_kmer = read_list[i_read][(i_kmer):(i_kmer + k)]
		# 	next_kmer = read_list[i_read][(i_kmer + 1):(i_kmer + 1 + k)]
		#
		# 	if curr_kmer not in self.dbg:
		# 		self.dbg.add_node(curr_kmer, read_list_n={i_read})
		# 	else:
		# 		self.dbg.node[curr_kmer]['read_list_n'].add(i_read)
		#
		# 	if next_kmer not in self.dbg:
		# 		self.dbg.add_node(next_kmer, read_list_n={i_read})
		# 	else:
		# 		self.dbg.node[next_kmer]['read_list_n'].add(i_read)
		#
		# 	self.dbg.add_edge(curr_kmer, next_kmer)

	def build_read_set_for_path(self, a_path, verbose=False):
		# a_path=map(lambda x:KMER_UID[x],a_path)
		missing_kmers = set(a_path).difference(self.kmer_map)
		if len(missing_kmers):
			# logger.critical("Completely missing kmer (%d): %s", len(missing_kmers), missing_kmers)
			return set()

		# current_set = set(self.kmer_map[a_path[0]]['read_list_n'])
		current_set = set(self.kmer_map[a_path[0]])
		if verbose:
			print len(current_set)

		assert isinstance(current_set, set)

		for i, a_node in enumerate(a_path[1:]):
			# if not self.dbg.has_edge(a_path[i], a_node):
			# 	logger.critical("Missing edge between %s -> %s",a_path[i],a_node)

			if not (a_path[i], a_node) in self.possible_pairs:
				logger.critical("Missing edge between %s -> %s. Actual read set is of size %d", a_path[i], a_node, len(current_set))


			# current_set.intersection_update(self.dbg.node[a_node]['read_list_n'])
			current_set.intersection_update(self.kmer_map[a_node])

			if not (a_path[i], a_node) in self.possible_pairs:
				logger.critical("Missing edge between %s -> %s. After update, read set is of size %d", a_path[i], a_node, len(current_set))

			if verbose:
				print len(current_set)
			if len(current_set) < 1:
				# premature exit, we reached an empty set 	
				return current_set
		return current_set

	def check_path(self, reference_path, alternative_path, min_cov):
		ref_path_read_set = self.build_read_set_for_path(reference_path)
		alt_path_read_set = self.build_read_set_for_path(alternative_path)
		# if len(alt_path_read_set) == 0:
		# 	logger.critical("Empty ALT read set for path %s", alternative_path)
		# if len(ref_path_read_set) == 0:
		# 	logger.critical("Empty REF read set for path %s", reference_path)
		if len(alt_path_read_set) + len(ref_path_read_set) <= min_cov:
			ratio = 0.0
		else:
			ratio = float(len(alt_path_read_set)) / (len(alt_path_read_set) + len(ref_path_read_set))
		return ratio, len(ref_path_read_set), len(alt_path_read_set)
