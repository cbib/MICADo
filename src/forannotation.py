import difflib
import re

from Bio import pairwise2

from helpers.logger import init_logger

logger = init_logger('Annotation')


def find_edit_operations(seq1, seq2):
	matcher = difflib.SequenceMatcher(a=seq1, b=seq2)
	return matcher.get_opcodes()


def alteration_list_to_transcrit_mutation(g_test, g_ref):
	annotated_alterations = []
	for i_alteration in range(0, len(g_test.significant_alteration_list)):
		# TODO substitute biopython global alignment by python difflib
		# TODO allow for multiple mutation to be returned as non distinguishable, e.g. they are reported as being in the same group.
		curr_alteration = g_test.significant_alteration_list[i_alteration]
		ref_seq = curr_alteration.reference_sequence
		alt_seq = curr_alteration.alternative_sequence
		# logger.info("Will perform alignment between \n %s \n %s", ref_seq, alt_seq)
		alignments = pairwise2.align.globalms(ref_seq, alt_seq, 2, -3, -5, -2)
		# if more than one alignment, choose the one which with the alteration at the left in the genome
		if len(alignments) > 1:
			# logger.critical("More than one alignment for %s vs %s", g_test.significant_alteration_list[i_alteration].reference_sequence,g_test.significant_alteration_list[i_alteration].alternative_sequence) 
			alignments = [alignments[0]]
		uncompact_cigar = ""
		compact_cigard = []
		for i_nucleotide in range(0, alignments[0][4]):
			if alignments[0][0][i_nucleotide] == alignments[0][1][i_nucleotide]:
				uncompact_cigar += "M"
			elif alignments[0][0][i_nucleotide] == "-":
				uncompact_cigar += "I"
			elif alignments[0][1][i_nucleotide] == "-":
				uncompact_cigar += "D"
			else:
				uncompact_cigar += "X"
		# print uncompact_cigar
		operation = uncompact_cigar[0]
		count = 0
		for i_nucleotide in range(0, len(uncompact_cigar)):
			if uncompact_cigar[i_nucleotide] != operation:
				compact_cigard += [count, operation]
				operation = uncompact_cigar[i_nucleotide]
				count = 0
			count += 1
		compact_cigard += [count, operation]
		# print compact_cigard

		if len(compact_cigard) != 6:
			logger.critical("More than one alteration: for %s", compact_cigard)
			continue

		alteration_type = compact_cigard[3]
		# print g_ref.dbg.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']
		ref_path_list = g_ref.dbg.node[curr_alteration.reference_path[0]]['ref_list']
		if len(ref_path_list) == 1:
			splicing_variant = ref_path_list.keys()[0]
		elif "NM_000546.5" in ref_path_list:
			splicing_variant = "NM_000546.5"
		elif "NM_001126114.2" in ref_path_list:
			splicing_variant = "NM_001126114.2"
		elif "NM_001126113.2" in ref_path_list:
			splicing_variant = "NM_001126113.2"
		position = ref_path_list[splicing_variant] + compact_cigard[0]
		if re.match("rs", splicing_variant):
			reference_sequence = ""
			# print splicing_variant
			# print compact_cigard
			# print g_ref.nt_ref[splicing_variant]
			for i_pos in range(0, compact_cigard[2]):
				if position + i_pos not in g_ref.nt_ref[splicing_variant]:
					logger.critical("%d not in nt_ref dict of %s", position + i_pos, splicing_variant)
					continue
				reference_sequence += g_ref.nt_ref[splicing_variant][position + i_pos]
			splicing_variant = "NM_000546.5"
		else:
			reference_sequence = ref_seq[compact_cigard[0]:compact_cigard[0] + compact_cigard[2]]


		alteration_description = {
			"splicing_variant": splicing_variant,
			"start": position,
			"end": None,
			'alt_sequence': None,
			"reference_sequence": reference_sequence,
			"alt_type": alteration_type,
			"pvalue_ratio": curr_alteration.pvalue_ratio,
			"alt_read_count": curr_alteration.alternative_read_count,
			"z_score": float(curr_alteration.zscore),

		}

		if alteration_type == "X":
			# c.76A>T
			if compact_cigard[2] != 1:
				logger.critical("More than one alteration: for %s", compact_cigard)
				continue
			alteration = alt_seq[compact_cigard[0]:compact_cigard[0] + 1]
			alteration_description['end'] = alteration_description['start']
			alteration_description['alt_sequence'] = alteration
			annotated_alterations.append(alteration_description)
			print "%s:c.%d%s>%s\t%d\t%d\t%f\t%f" % (splicing_variant, position, reference_sequence, alteration, curr_alteration.reference_read_count,
													curr_alteration.alternative_read_count,
													curr_alteration.pvalue_ratio,
													float(curr_alteration.zscore))
		elif alteration_type == "D":
			# c.76_78delACT
			alteration_description['end'] = position + len(reference_sequence) - 1
			annotated_alterations.append(alteration_description)
			print "%s:c.%d_%ddel%s\t%d\t%d\t%f\t%f" % (
				splicing_variant, position, position + len(reference_sequence) - 1, reference_sequence, curr_alteration.reference_read_count,
				curr_alteration.alternative_read_count, curr_alteration.pvalue_ratio,
				float(curr_alteration.zscore))
		else:
			# c.76_77insG
			alteration = alt_seq[compact_cigard[0]:compact_cigard[0] + compact_cigard[2]]
			alteration_description['alt_sequence'] = alteration
			alteration_description['end'] = position + 1
			annotated_alterations.append(alteration_description)

			print "%s:c.%d_%dins%s\t%d\t%d\t%f\t%f" % (
				splicing_variant, position, position + 1, alteration, curr_alteration.reference_read_count,
				curr_alteration.alternative_read_count, curr_alteration.pvalue_ratio,
				float(curr_alteration.zscore))

	return annotated_alterations

def splicing_variant_converter(spv):
	if spv not in ["NM_000546.5", "NM_001126114.2", "NM_001126113.2"]:
		spv = "NM_000546.5"
	return spv
