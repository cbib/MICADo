import collections
import difflib
from unittest import TestCase
from Bio import pairwise2
import randomreadsgraph
import numpy as np
import reference_graph
import patient_graph
import Levenshtein

__author__ = 'hayssam'

negative_controls = ["C_158_1", "C_193_1", "C_256_1", "C_267_1", "C_284_1", "C_285_1", "C_288_1", "C_306_1", "C_312_1", "C_319_1",
					 "C_322_1", "N_158_1", "N_193_1", "N_256_1", "N_267_1", "N_284_1", "N_285_1", "N_288_1", "N_306_1", "N_312_1",
					 "N_319_1", "N_322_1"]
positive_controls = ["C_169_1", "C_169_2", "C_207_1", "C_207_2", "C_221_1", "C_221_2", "C_279_1", "C_279_2", "C_316_1", "C_316_2",
					 "C_318_1", "C_318_2", "C_320_1", "C_320_2", "C_340_1", "C_340_2", "C_341_1", "C_341_2", "N_169_1", "N_169_2",
					 "N_207_1", "N_207_2", "N_221_1", "N_221_2", "N_279_1", "N_279_2", "N_316_1", "N_316_2", "N_318_1", "N_318_2",
					 "N_320_1", "N_320_2", "N_340_1", "N_340_2", "N_341_1", "N_341_2"]


def micado_multi(sample_key, n_perm=25):
	kmer_length = 18
	max_len = 10
	# build reference graph
	g_reference = reference_graph.ReferenceGraph(kmer_length,
												 fasta_file='data/reference/NM_000546.5.fasta',
												 snp_file='data/reference/snp_TP53.tab')
	# build patient graph
	g_patient = patient_graph.SampleGraph(['data/tp53_analysis/reads/%s.fastq' % sample_key], kmer_length)
	g_patient.graph_cleaned_init(3.0)
	# copy g_patient cleaned and remove reference edges on it (.dbg_refrm creation)
	g_patient.graph_remove_reference_edges(g_patient.dbgclean, g_reference.dbg)
	# search for alternative paths in dbg_refrm (.alteration_list creation)
	g_patient.alteration_list_init(g_reference.dbg, kmer_length, 3.0, max_len)

	# TODO build real set of possible k-mers
	all_possible_kmers = set()
	for an_alt in g_patient.alteration_list:
		all_possible_kmers.update(an_alt.reference_path)
		all_possible_kmers.update(an_alt.alternative_path)

	# build a random read graph
	import seq_lib_TP53 as seq_lib
	random_ratio_dict = collections.defaultdict(list)
	lonely_ratio_dict = {}
	ref_seq_dict = {}
	alt_seq_dict = {}
	for n_perm in range(n_perm):
		print n_perm
		rg = randomreadsgraph.RandomReadsGraph({"N": 0, "C": 0}, k=kmer_length, seq_lib_module=seq_lib, restrict_to=None)
		for alt_i, putative_alt in enumerate(g_patient.alteration_list):
			# determine number of edit ops
			# There's at least one (since it's an alternative path)
			ref_seq = putative_alt.reference_sequence
			ref_seq_dict[alt_i] = ref_seq
			patient_seq = putative_alt.alternative_sequence
			edit_ops = Levenshtein.editops(ref_seq, patient_seq)
			lonely_ratio = putative_alt.ratio_read_count
			lonely_ratio_dict[alt_i] = lonely_ratio
			# print n_perm, alt_i, lonely_ratio, edit_ops
			for e in edit_ops:
				# print "Considering atomic edit op", e
				transformed = Levenshtein.apply_edit([e], ref_seq, patient_seq)
				ratio_random = rg.check_path(kmerize(ref_seq, kmer_length),
											 kmerize(transformed, kmer_length),
											 min_cov=putative_alt.min_coverage)
				random_ratio_dict[(alt_i, (e,))].append(ratio_random[0])
				alt_seq_dict[(alt_i, (e,))] = transformed
			# perform for all edit_ops
			ratio_random = rg.check_path(kmerize(ref_seq, kmer_length), kmerize(patient_seq, kmer_length),
										 min_cov=putative_alt.min_coverage)
			random_ratio_dict[(alt_i, tuple(edit_ops))].append(ratio_random[0])
			alt_seq_dict[(alt_i, tuple(edit_ops))] = patient_seq
	for (alt_i, edit_ops_i), ratios in sorted(random_ratio_dict.items(), key=lambda x: x[0][0]):
		this_patient_ratio = lonely_ratio_dict[alt_i]
		random_ratios = alt_seq_dict[(alt_i, edit_ops_i)]
		print "Alt %d with real ratio %f, Edit ops %s, random_ratios :%s" % (alt_i, this_patient_ratio, edit_ops_i, map(str, ratios))
		print "Ref seq %s" % (ref_seq_dict[alt_i])
		print "Alt seq %s" % (random_ratios)
		print "N higher: %d" % (len([x for x in ratios if x > this_patient_ratio]))
		standard_deviation = np.std(ratios)
		zscore = float((this_patient_ratio - np.mean(ratios)) / np.std(ratios))
		print "Z-Score", zscore


def kmerize(s, k):
	if (k + 2) <= len(s):
		return [s[i:(i + k)] for i in range(0, len(s) - (k - 1))]
	else:
		return [s]


def align(ref, transformed):
	return pairwise2.align.globalms(ref, transformed, 2, -3, -5, -2)[0]


class TestRandomReadsGraph(TestCase):
	def test_kmerize(self):
		self.assertEqual(len(kmerize("A" * 30, 20)), 10)

	def test_op_edits_for_N_193_1(self):
		ref = "ATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCC"
		alt = "ATGCCAGAGGCTGCTCCCGCGTGGCCCTGCACCAGCAGCTCC"
		# matcher = difflib.SequenceMatcher(a=ref, b=alt)
		# print matcher.get_opcodes()
		# op = [x[0] for x in matcher.get_opcodes() if x[0] != 'equal']
		# print op
		# alignments = pairwise2.align.globalms(ref, alt, 2, -3, -5, -2)
		# print alignments
		#
		# matcher2= difflib.SequenceMatcher(a="CCC",b="GC")
		# print matcher2.get_opcodes()

		editops = Levenshtein.editops(ref, alt)
		print editops
		# print opcodes
		# print Levenshtein.apply_edit(opcodes,ref,alt)
		for e in editops:
			print "applying", e
			try:
				transformed = Levenshtein.apply_edit([e], ref, alt)
				print align(ref, transformed)
			except Exception:
				print "Fail"

	def test_build_read_set_for_path_N_193_1(self):
		micado_multi("N_193_1")

	def test_build_read_set_for_path_N_183_1(self):
		micado_multi("N_183_1")  # tips

	def test_C_221_2(self):
		micado_multi("C_221_2")  # tips

	def test_build_read_set_for_path_N_192_2(self):
		micado_multi("N_192_2")

	def test_build_read_set_for_path_N_207_2(self):
		micado_multi("N_207_2")

	def test_build_read_set_for_path_N_215_1(self):
		micado_multi("N_215_1")

	def test_build_read_set_for_path_N_272_1(self):
		micado_multi("N_272_1")

	def test_build_read_set_for_all_negative_controls(self):
		for neg_sample in negative_controls:
			print "Processing sample", neg_sample
			micado_multi(neg_sample)

	def test_build_read_set_for_all_positive_controls(self):
		for pos_sample in positive_controls:
			print "Processing sample", pos_sample
			micado_multi(pos_sample)

	def test_rrg_creation_speed(self):
		import seq_lib_TP53 as seq_lib
		for i in range(10):
			rg = randomreadsgraph.RandomReadsGraph({"N": 0, "C": 0}, k=18, seq_lib_module=seq_lib, restrict_to=None)
			print i

	def test_sampling_speed(self):
		import seq_lib_TP53 as seq_lib
		for i in range(10):
			read_list = seq_lib.sampling({"N": 0, "C": 0})
			print i

# ref = "ATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCC"
# alt = "ATGCCAGAGGCTGCTCCCGCGTGGCCCTGCACCAGCAGCTCC"
#

# found_ratios =
# ratio_random = rg.check_path(kmerize(ref, 18), kmerize(alt, 18), min_cov=1)
# print ratio_random
