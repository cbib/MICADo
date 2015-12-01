#!/usr/bin/env python
# coding=utf8

import numpy as np


class alteration:
	def __init__(self, reference_path, alternative_path, reference_sequence, alternative_sequence, reference_read_count, alternative_read_count, k, min_cov):
		self.reference_path = reference_path
		self.alternative_path = alternative_path
		self.reference_sequence = reference_sequence
		self.alternative_sequence = alternative_sequence
		self.reference_read_count = reference_read_count
		self.alternative_read_count = alternative_read_count
		self.min_coverage = min_cov

		self.ratio_read_count = float(alternative_read_count) / (reference_read_count + alternative_read_count)
		self.random_ratio_list = []
		self.random_reference_count_list = []
		self.random_alternative_count_list = []


		print self.reference_sequence,reference_read_count
		print self.alternative_sequence,alternative_read_count

	def pvalue_init(self):
		# self.pvalue_count = float(len([i for i in self.random_alternative_count_list if i >= self.alternative_read_count])) / len(self.random_alternative_count_list)
		self.pvalue_ratio = float(len([i for i in self.random_ratio_list if i >= self.ratio_read_count])) / len(self.random_ratio_list)
		standard_deviation = np.std(self.random_ratio_list)
		if standard_deviation == 0.0:
			self.zscore = "inf"
		else:
			self.zscore = float((self.ratio_read_count - np.mean(self.random_ratio_list)) / np.std(self.random_ratio_list))

	@staticmethod
	# Kmerize a sequence, return a kmer list
	def kmerpathToSeq(kmer_list, k):
		sequence = kmer_list[0]
		for i_kmer in range(1, len(kmer_list)):
			sequence += kmer_list[i_kmer][k - 1]
		return sequence
	@staticmethod
	# kmerize a sequence, return a kmer list
	def kmerize(sequence, kmer_length):
		return [sequence[i:i + kmer_length] for i in range(0, len(sequence) - kmer_length)]
