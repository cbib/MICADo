#! /usr/bin/env python
# encoding: utf-8
# FORMAT = '%(asctime)-15s %(message)s'
# logging.basicConfig(format=FORMAT)
# logger = logging.getLogger('mix')

# This script is used to construct snp list for MICADo 
# Example for the FLT3 transcript: python SNPfile_makor.py NM_004119.2.fasta refseq_22221.txt snp_FLT3.tab NM_004119 200

import os
import sys
import logging
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

FASTA = sys.argv[1]
REFSEQ = sys.argv[2]
OUT_TAB = sys.argv[3]
SPLICING_VARIANT = sys.argv[4]
NT_AROUND_SNP = int(sys.argv[5])

print "SNPfile_makor.py: %s %s %s %d"%(FASTA,REFSEQ,OUT_TAB,NT_AROUND_SNP)

OUT_TAB = open(OUT_TAB,'w')
lines = open(REFSEQ, 'r').readlines()

for record in SeqIO.parse(FASTA, "fasta", generic_dna):
	record.name = record.id.split("|")[3]
	for i in range(1,len(lines)):
		line_ = lines[i].rstrip()
		line_split = line_.split("\t")
		if len(line_split) != 14:
			break
		if line_split[4] == SPLICING_VARIANT:
			split_nt_changes = line_split[1].split("|")
			print "Write info for SNP %s on %s"%(line_split[0],SPLICING_VARIANT)
			for i in range(1,len(split_nt_changes)):
				alt_pos = int(line_split[7]) - 1 # index of nucleotide (1st nt of have index 0)
				OUT_TAB.write("%s\t%d\n%s\t%s\t%s\n"%(line_split[0],alt_pos-NT_AROUND_SNP,str(record.seq[alt_pos-NT_AROUND_SNP:alt_pos]),str(Seq(record.seq[alt_pos],generic_dna))+"|"+str(Seq(split_nt_changes[i],generic_dna).reverse_complement()),str(record.seq[alt_pos+1:alt_pos+1+NT_AROUND_SNP])))

OUT_TAB.close()