#!/usr/bin/env python
# coding=utf8
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import networkx as nx

class ReferenceGraph:
	def __init__(self,kmer_length,fasta_file,snp_file):
		self.dbg = nx.DiGraph()
		self.nt_ref = {}
		# Reference
		for record in SeqIO.parse(fasta_file, "fasta", generic_dna):
			startposition = int(record.description.split("\t")[1]) + 1
			seq_s = str(record.seq)
			for i2 in range(0, len(seq_s) - kmer_length):
				curr_kmer = seq_s[(i2):(i2 + kmer_length)]
				next_kmer = seq_s[(i2 + 1):(i2 + 1 + kmer_length)]
				if next_kmer not in self.dbg:
					self.dbg.add_node(next_kmer, ref_list={record.id:startposition+i2})
				if curr_kmer in self.dbg:
					if self.dbg[curr_kmer].get(next_kmer, 0) == 0:
						self.dbg.add_edge(curr_kmer, next_kmer)
				else:
					# for the first kmer
					self.dbg.add_node(curr_kmer, ref_list={record.id:startposition+i2})
					self.dbg.add_edge(curr_kmer, next_kmer)
		# SNP
		IN_SNP = open(snp_file, 'r')
		lines = IN_SNP.readlines()
		lines = map(str.strip, lines)
		for l in range(1,len(lines),+2):
			line_split = lines[l].split("\t")
			line_before_split = lines[l-1].split("\t")
			start_position = int(line_before_split[1]) + len(line_split[0])-kmer_length + 1
			kmer_around_ref = line_split[0][len(line_split[0])-kmer_length:len(line_split[0])]+line_split[1].split("|")[0]+line_split[2][0:kmer_length]
			kmer_around_snp = line_split[0][len(line_split[0])-kmer_length:len(line_split[0])]+line_split[1].split("|")[1]+line_split[2][0:kmer_length]
			for i2 in range(0, len(kmer_around_snp) - kmer_length):
				self.nt_ref[line_before_split[0]]={start_position+i2:kmer_around_snp[i2]}
				curr_kmer = kmer_around_snp[(i2):(i2 + kmer_length)]
				next_kmer = kmer_around_snp[(i2 + 1):(i2 + 1 + kmer_length)]
				if next_kmer not in self.dbg:
					self.dbg.add_node(next_kmer, ref_list={line_before_split[0]:start_position+i2})
				# elif line_before_split[0] not in self.dbg.node[next_kmer]['ref_list']:
				# 	self.dbg.node[next_kmer]['ref_list'][line_before_split[0]] = start_position+i2
				if curr_kmer in self.dbg:
					if self.dbg[curr_kmer].get(next_kmer, 0) == 0:
						self.dbg.add_edge(curr_kmer, next_kmer)
				else:
					self.dbg.add_node(curr_kmer, ref_list={line_before_split[0]:start_position+i2})
					self.dbg.add_edge(curr_kmer, next_kmer)
			self.nt_ref[line_before_split[0]]={}
			for i2 in range(0, len(kmer_around_snp)):
				self.nt_ref[line_before_split[0]][start_position+i2]=kmer_around_ref[i2]
