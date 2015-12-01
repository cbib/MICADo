#!/usr/bin/env python
# coding=utf8
import collections
import networkx as nx
from helpers.logger import init_logger

logger = init_logger("RANDGRAPH")

cached_kmers = {}
cached_pairs = {}
KMER_UID = {}
curr_uid = 0
last_sample = None


def kmerize_iter(s, k):
	for i in range(0, len(s) - k):
		yield s[i:i + k]


class RandomReadsGraph:
	def __init__(self, coverage_dict, k, seq_lib_module, restrict_to=None):
		global cached_kmers, curr_uid, last_sample
		self.coverage_dict = coverage_dict
		# TODO optimize indexation
		self.kmer_map = collections.defaultdict(set)
		self.restrict_to = set(restrict_to) if restrict_to else None
		self.possible_pairs = set()
		read_list = seq_lib_module.sampling(self.coverage_dict)
		self.dbg = nx.DiGraph()

		for i_read in range(0, len(read_list)):
			this_read = read_list[i_read]
			if self.restrict_to:
				kkmers = [this_read[i:i + k] for i in xrange(len(this_read) - k) if this_read[i:i + k] in self.restrict_to]
			else:
				kkmers = [this_read[i:i + k] for i in xrange(len(this_read) - k) if this_read[i:i + k]]
			kmers_pairs = [(f1, f2) for f1, f2 in zip(kkmers, kkmers[1:])]

			for kmer in kkmers:
				self.kmer_map[kmer].add(i_read)

			self.possible_pairs.update(kmers_pairs)
		# print len(self.kmer_map), np.mean(map(len, self.kmer_map.values()))

	def build_read_set_for_path(self, a_path, verbose=False):
		missing_kmers = set(a_path).difference(self.kmer_map)
		if len(missing_kmers):
			# logger.critical("Completely missing kmer (%d): %s", len(missing_kmers), missing_kmers)
			return set()
		current_set = set(self.kmer_map[a_path[0]])
		if verbose:
			print len(current_set)

		assert isinstance(current_set, set)

		for i, a_node in enumerate(a_path[1:]):
			if not (a_path[i], a_node) in self.possible_pairs:
				logger.critical("Missing edge between %s -> %s. Actual read set is of size %d", a_path[i], a_node, len(current_set))

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
			logger.critical("Total coverage %d below min_cov %d", len(alt_path_read_set) + len(ref_path_read_set), min_cov)
			ratio = 0.0
		else:
			ratio = float(len(alt_path_read_set)) / (len(alt_path_read_set) + len(ref_path_read_set))
		return ratio, len(ref_path_read_set), len(alt_path_read_set)
