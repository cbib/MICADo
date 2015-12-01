from unittest import TestCase
from Bio import SeqIO
from Bio.Alphabet import generic_dna

__author__ = 'hayssam'


class TestPatientGraph(TestCase):
	def test_fastq_parser(self):
		f = "data/unit_test_data/C_SYNTHP53_22640_500_05_3_1-1-1.fastq"
		for seq in SeqIO.parse(f, "fastq", generic_dna):
			print seq
