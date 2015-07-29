# For one PCR amplicon
#!/usr/bin/env python
# coding=utf8
import collections
import networkx as nx
from helpers.logger import init_logger
import seq_lib as SL

logger = init_logger("RANDGRAPH")

class RandomReadsGraph:
	def __init__(self, coverage_dict, k,restrict_to={}):
		self.coverage = sum(coverage_dict.values())
		self.kmer_map=collections.defaultdict(set)
		self.restrict_to=set(restrict_to)
		self.possible_pairs=set()
		read_list = SL.sampling(self.coverage)
		self.dbg = nx.DiGraph()
		for i_read in range(0, len(read_list)):
			this_read=read_list[i_read]
			kkmers=[this_read[i:i+k] for i in xrange(len(this_read)-k) if this_read[i:i+k] in self.restrict_to]
			kmers_pairs=[(f1,f2) for f1,f2 in zip(kkmers,kkmers[1:])]
			for kmer in kkmers:
				self.kmer_map[kmer].add(i_read)
			self.possible_pairs.update(kmers_pairs)

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

			if not (a_path[i],a_node) in self.possible_pairs:
				logger.critical("Missing edge between %s -> %s. Actual read set is of size %d",a_path[i],a_node,len(current_set))


			# current_set.intersection_update(self.dbg.node[a_node]['read_list_n'])
			current_set.intersection_update(self.kmer_map[a_node])

			if not (a_path[i],a_node) in self.possible_pairs:
				logger.critical("Missing edge between %s -> %s. After update, read set is of size %d",a_path[i],a_node,len(current_set))

			if verbose:
				print len(current_set)
			if len(current_set)<1:
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


