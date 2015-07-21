#!/usr/bin/env python
# coding=utf8
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import networkx as nx

def ref_constructor(kmer_length,fasta_file):
	g_reference = nx.DiGraph()
	for record in SeqIO.parse(fasta_file, "fasta", generic_dna):
		seq_s = str(record.seq)
		for i2 in range(0, len(seq_s) - kmer_length):
			curr_kmer = seq_s[(i2):(i2 + kmer_length)]
			next_kmer = seq_s[(i2 + 1):(i2 + 1 + kmer_length)]
			if next_kmer not in g_reference:
				g_reference.add_node(next_kmer, ref_list=[record.id])
			else:
				g_reference.node[next_kmer]['ref_list'].append(record.id)	
			if curr_kmer in g_reference:
				g_reference.node[curr_kmer]['ref_list'].append(record.id)
				if g_reference[curr_kmer].get(next_kmer, 0) == 0:
					g_reference.add_edge(curr_kmer, next_kmer)
			else:
				g_reference.add_node(curr_kmer, ref_list=[record.id])
				g_reference.add_edge(curr_kmer, next_kmer)
	return g_reference