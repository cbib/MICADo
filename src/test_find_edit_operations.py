from unittest import TestCase
from Bio import pairwise2
from forannotation import find_edit_operations

__author__ = 'hayssam'


class TestFind_edit_operations(TestCase):
	def test_find_edit_operations(self):
		ref="CTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCAC"
		alt="CTGGGAGAGACCGGCGCACAccGAGGAAGAGAATCTCcGCAAGAAAGGGGAGCCTCAC"
		print find_edit_operations(ref,alt)
		# compare to
		# alignments = pairwise2.align.globalms(ref, alt, 2, -3, -5, -2)
		# %timeit pairwise2.align.globalms(ref, alt, 2, -3, -5, -2)
