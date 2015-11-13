from unittest import TestCase
from read_sampler.altered_reads_sampler import parse_sam_file,build_a_sample

__author__ = 'hayssam'


class TestBuild_a_sample(TestCase):
	def test_build_a_sample(self):
		# tgt is C_SYNTHP53_22640_500_05_3_1-1-1.fastq
		source_read = parse_sam_file("data/experimental_results/TP53/alignments/C_model_GMAPno40_NM_000546.5.sam")

