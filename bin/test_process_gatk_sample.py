from unittest import TestCase
from bin.aggregate_caller_results import process_gatk_sample

__author__ = 'hayssam'


class TestProcess_gatk_sample(TestCase):
	def test_process_gatk_sample(self):
		a_sample='C_SYNTHP53_23406_500_10_3_1-1-1'
		process_gatk_sample(a_sample)

