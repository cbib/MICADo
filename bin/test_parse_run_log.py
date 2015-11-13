from unittest import TestCase
from bin.aggregate_caller_results import parse_run_log

__author__ = 'hayssam'


class TestParse_run_log(TestCase):
	def test_parse_run_log(self):
		a_log = "exec_logs/micado_log_C_SYNTHP53_10010_500_05_3_1-1-1.txt"
		print parse_run_log(a_log)
