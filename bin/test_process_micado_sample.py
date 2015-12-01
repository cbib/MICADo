from unittest import TestCase
from bin.tabulate_and_aggregate_xp_results import process_micado_sample

__author__ = 'hayssam'


class TestProcess_micado_sample(TestCase):
	def test_process_micado_sample(self):
		sample_tag = "C_SYNTHP53_10010_500_05_3_1-1-1"
		process_micado_sample(sample_tag)
