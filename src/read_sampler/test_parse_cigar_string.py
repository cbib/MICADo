from unittest import TestCase
from helpers.helpers import data_dir
from read_sampler.altered_reads_sampler import parse_sam_file
from read_sampler.cigar_parser import parse_cigar_string

__author__ = 'hayssam'


class TestParse_cigar_string(TestCase):
	def test_parse_cigar_string(self):
		starting_file = data_dir + "/alignments/C_model_GMAPno40_NM_000546.5.sam"
		reads = parse_sam_file(starting_file)
		for i,r in enumerate(reads.CIGAR):
			print i,parse_cigar_string(r)



