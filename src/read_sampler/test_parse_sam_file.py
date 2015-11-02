from unittest import TestCase
from helpers.helpers import data_dir
from read_sampler.altered_reads_sampler import parse_sam_file

__author__ = 'hayssam'


class TestParse_sam_file(TestCase):
	def test_parse_sam_file(self):
		starting_file = data_dir + "/alignments/C_model_GMAPno40_NM_000546.5.sam"
		reads = parse_sam_file(starting_file)
		# print reads
		print reads.sample(10)