#!/usr/bin/env python
# coding=utf8
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import networkx as nx


class ReferenceGraph:
	def __init__(self, kmer_length, fasta_file, snp_file):
		self.dbg = nx.DiGraph()
		self.nt_ref = {}
		self.snp = {}
		self.kmer_length = kmer_length
		# Reference
		for record in SeqIO.parse(fasta_file, "fasta", generic_dna):
			try:
				startposition = int(record.description.split("\t")[1]) + 1
			except IndexError:
				startposition = 0
			seq_s = str(record.seq)
			for i2 in range(0, len(seq_s) - kmer_length):
				curr_kmer = seq_s[(i2):(i2 + kmer_length)]
				next_kmer = seq_s[(i2 + 1):(i2 + 1 + kmer_length)]
				if next_kmer not in self.dbg:
					self.dbg.add_node(next_kmer, ref_list={record.id: startposition + i2})
				if curr_kmer in self.dbg:
					if self.dbg[curr_kmer].get(next_kmer, 0) == 0:
						self.dbg.add_edge(curr_kmer, next_kmer)
				else:
					# for the first kmer
					self.dbg.add_node(curr_kmer, ref_list={record.id: startposition + i2})
					self.dbg.add_edge(curr_kmer, next_kmer)
		# if snp_file:
		# 	self.inject_known_snps(snp_file)

	def inject_known_snps(self, snp_file):
		k = self.kmer_length
		# SNP
		IN_SNP = open(snp_file, 'r')
		lines = IN_SNP.readlines()
		lines = map(str.strip, lines)
		for l in range(1, len(lines), +2):
			line_split = lines[l].split("\t")
			line_before_split = lines[l - 1].split("\t")
			left = line_split[0]
			left_len = len(left)
			start_position = int(line_before_split[1]) + left_len - k + 1
			snp_code = line_split[1]
			right = line_split[2]
			reference_base = snp_code.split("|")[0]
			kmer_around_ref = left[left_len - k:left_len] + reference_base + right[0:k]
			snp_base = snp_code.split("|")[1]

			kmer_around_snp = left[left_len - k:left_len] + snp_base + right[0:k]
			for i2 in range(0, len(kmer_around_snp) - k):
				self.nt_ref[line_before_split[0]] = {start_position + i2: kmer_around_snp[i2]}
				curr_kmer = kmer_around_snp[i2:(i2 + k)]
				next_kmer = kmer_around_snp[(i2 + 1):(i2 + 1 + k)]
				if next_kmer not in self.dbg:
					self.dbg.add_node(next_kmer, ref_list={line_before_split[0]: start_position + i2})
				# elif line_before_split[0] not in self.dbg.node[next_kmer]['ref_list']:
				# 	self.dbg.node[next_kmer]['ref_list'][line_before_split[0]] = start_position+i2
				if curr_kmer in self.dbg:
					if self.dbg[curr_kmer].get(next_kmer, 0) == 0:
						self.dbg.add_edge(curr_kmer, next_kmer)
				else:
					self.dbg.add_node(curr_kmer, ref_list={line_before_split[0]: start_position + i2})
					self.dbg.add_edge(curr_kmer, next_kmer)
			self.nt_ref[line_before_split[0]] = {}
			for i2 in range(0, len(kmer_around_snp)):
				self.nt_ref[line_before_split[0]][start_position + i2] = kmer_around_ref[i2]
			self.snp[line_before_split[0]] = [
				left[left_len - k + 1:left_len] + reference_base,
				left[left_len - k + 1:left_len] + snp_base]
