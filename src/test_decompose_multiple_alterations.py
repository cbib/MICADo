from unittest import TestCase
from patient_graph import decompose_multiple_alterations
from test_randomReadsGraph import kmerize

__author__ = 'hayssam'


class TestDecompose_multiple_alterations(TestCase):
	def test_decompose_multiple_alterations(self):
		ref_seq = "CAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		ref_path = kmerize(ref_seq, 18)
		alt_seq = "CAGGTCCAGATGAAGCTCGCTCCCCGCGTGGCCCCTGCACCAGC"
		alt_path = kmerize(alt_seq, 18)
		decomposed_alt = list(decompose_multiple_alterations(ref_path, alt_path, 18))
		for alt_seq, at_path in decomposed_alt:
			print ref_seq, alt_seq
		self.assertEqual(len(decomposed_alt), 2)

	def test_decompose_single_alterations_with_deletions(self):
		ref_seq =  "CAGGTCCAGATGAAGCTCCXXXXXCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		ref_path = kmerize(ref_seq, 18)
		alt_seq = "CAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		alt_path = kmerize(alt_seq, 18)
		decomposed_alt = list(decompose_multiple_alterations(ref_path, alt_path, 18))
		for alt_seq, at_path in decomposed_alt:
			print ref_seq, alt_seq
			self.assertEqual(len([x for x in alt_seq if x == "X"]), 0)  # at most 5 alt
		self.assertEqual(len(decomposed_alt), 1)

	def test_decompose_multiple_alterations_with_deletions(self):
		ref_seq =  "CAGGTCCAGATGAAGCTCCXXXXXCAGAATGCCAGAXXXXXGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		ref_path = kmerize(ref_seq, 18)
		alt_seq = "CAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		alt_path = kmerize(alt_seq, 18)
		decomposed_alt = list(decompose_multiple_alterations(ref_path, alt_path, 18))
		for alt_seq, at_path in decomposed_alt:
			print ref_seq, alt_seq
			self.assertEqual(len([x for x in alt_seq if x == "X"]), 5)  # at most 0 X in the alt
		self.assertEqual(len(decomposed_alt), 2)

	def test_decompose_single_alterations_with_insertions(self):
		ref_seq =  "CAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		ref_path = kmerize(ref_seq, 18)
		alt_seq = "CAGGTCCAGATGAAGCTCCXXXXXCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		alt_path = kmerize(alt_seq, 18)
		decomposed_alt = list(decompose_multiple_alterations(ref_path, alt_path, 18))
		for alt_seq, at_path in decomposed_alt:
			print ref_seq, alt_seq
			self.assertEqual(len([x for x in alt_seq if x == "X"]), 5)  #at most 5 alt
		self.assertEqual(len(decomposed_alt), 1)

	def test_decompose_multiple_alterations_with_insertions(self):
		ref_seq =  "CAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		ref_path = kmerize(ref_seq, 18)
		alt_seq = "CAGGTCCAGATGAAGCTCCXXXXXCAGAATGCCAGAGGXXXXXCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		alt_path = kmerize(alt_seq, 18)
		decomposed_alt = list(decompose_multiple_alterations(ref_path, alt_path, 18))
		for alt_seq, at_path in decomposed_alt:
			print ref_seq, alt_seq
			self.assertEqual(len([x for x in alt_seq if x == "X"]), 5)  #at most 5 alt
		self.assertEqual(len(decomposed_alt), 2)

	def test_decompose_multiple_alterations_with_end_insertions(self):
		ref_seq =  "CAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		ref_path = kmerize(ref_seq, 18)
		alt_seq = "CAGGTCCAGATGAAGCTCCXXXXXCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCXXXXX"
		alt_path = kmerize(alt_seq, 18)
		decomposed_alt = list(decompose_multiple_alterations(ref_path, alt_path, 18))
		for alt_seq, at_path in decomposed_alt:
			print ref_seq, alt_seq
			self.assertEqual(len([x for x in alt_seq if x == "X"]), 5)  #at most 5 alt
		self.assertEqual(len(decomposed_alt), 2)

	def test_decompose_single_alterations_with_end_insertions(self):
		ref_seq =  "AAAA"
		ref_path = kmerize(ref_seq, 18)
		alt_seq = "XXXXXAAXXXXXAAXXXXX"
		alt_path = kmerize(alt_seq, 18)
		decomposed_alt = list(decompose_multiple_alterations(ref_path, alt_path, 18))
		for alt_seq, at_path in decomposed_alt:
			print ref_seq, alt_seq
			self.assertEqual(len([x for x in alt_seq if x == "X"]), 5)  #at most 5 alt
		self.assertEqual(len(decomposed_alt), 3)

	def test_decompose_single_alterations(self):
		ref_seq = "CAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		ref_path = kmerize(ref_seq, 18)
		alt_seq = "CAGGTCCAGATGAAGCTCCCXGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		alt_path = kmerize(alt_seq, 18)
		decomposed_alt = list(decompose_multiple_alterations(ref_path, alt_path, 18))
		for alt_seq, at_path in decomposed_alt:
			print ref_seq, alt_seq
		self.assertEqual(len(decomposed_alt), 1)

	def test_decompose_multiple_single_alterations(self):
		ref_seq = "CAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		ref_path = kmerize(ref_seq, 18)
		alt_seq = "CAGGTCCAGATGAAGCTCCCXGAATGCCAGXGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		alt_path = kmerize(alt_seq, 18)
		decomposed_alt = list(decompose_multiple_alterations(ref_path, alt_path, 18))
		for alt_seq, at_path in decomposed_alt:
			print ref_seq, alt_seq
			self.assertEqual(len([x for x in alt_seq if x == "X"]), 1)  # at most 1 alt
		self.assertEqual(len(decomposed_alt), 2)

	def test_decompose_multiple_mismatches(self):
		ref_seq = "CAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGC"
		ref_path = kmerize(ref_seq, 18)
		alt_seq = "CAGGTCCAGATGAAGCTCXXAGAATGCCAGXGGCTGCTCCCCXCGTGGCCCCTGCACCAGC"
		alt_path = kmerize(alt_seq, 18)
		decomposed_alt = list(decompose_multiple_alterations(ref_path, alt_path, 18))
		for alt_seq, at_path in decomposed_alt:
			print ref_seq, alt_seq
			self.assertEqual(len([x for x in alt_seq if x == "X"]), 1)  # at most 1 alt
			self.assertEqual(len(ref_seq), len(alt_seq))  # at most 1 alt
		self.assertEqual(len(decomposed_alt), 4)
