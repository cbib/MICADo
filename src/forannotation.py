import difflib
import re

from Bio import pairwise2
import itertools
from helpers import intset

from helpers.logger import init_logger

logger = init_logger('Annotation')


def find_edit_operations(seq1, seq2):
	matcher = difflib.SequenceMatcher(a=seq1, b=seq2)
	return matcher.get_opcodes()


def merge_identical_alterations(annotated_alterations):
	# identify groups having exactly the same start and end
	annotated_alterations.sort(key=lambda alt: (alt['start'], alt['end']))
	merged_alterations = []
	for g, alt_group in itertools.groupby(annotated_alterations, key=lambda alt: (alt['start'], alt['end'])):
		best_alt = max(list(alt_group), key=lambda alt: alt['alt_read_count'])
		merged_alterations.append(best_alt)
	# identify alterations belonging to the same interval
	all_ranges = [(x['start'], x['end']) for x in merged_alterations]
	logger.info("Will merge ranges: %s",all_ranges)
	int_sets = intset.IntSet(*all_ranges)
	normalized_ranges = [intset.IntSet(r) for r in int_sets._ranges]
	# assign each alteration to a range
	for alt in merged_alterations:
		this_range = (alt['start'], alt['end'])
		for i_range, a_range in enumerate(normalized_ranges):
			if a_range.overlaps(intset.IntSet(this_range)):
				alt['position_cluster'] = i_range
				break
	#
	one_alteration_per_cluster = []
	merged_alterations.sort(key=lambda x:x['position_cluster'])
	for k, alt_group in itertools.groupby(merged_alterations,key=lambda x:x['position_cluster']):
		representative_alt = max(alt_group,key=lambda x:x['alt_read_count'])
		one_alteration_per_cluster.append(representative_alt)
	#
	#
	# final_list = []
	# merged_alterations.sort(key=lambda alt: (alt['start']))
	# for g, alt_group in itertools.groupby(annotated_alterations,key=lambda alt: (alt['start'], alt['end'])):
	# 	best_alt = max(list(alt_group), key=lambda alt: alt['alt_read_count'])
	# 	merged_alterations.append(best_alt)

	return one_alteration_per_cluster


def alteration_list_to_transcrit_mutation(g_test, reference_graph):
	annotated_alterations = []
	for i_alteration in range(0, len(g_test.significant_alteration_list)):
		is_multi = False
		curr_alteration = g_test.significant_alteration_list[i_alteration]
		ref_seq = curr_alteration.reference_sequence
		alt_seq = curr_alteration.alternative_sequence
		alignments = pairwise2.align.globalms(ref_seq, alt_seq, 2, -3, -5, -2)
		if len(alignments) > 1:
			# logger.critical("More than one alignment for %s vs %s", g_test.significant_alteration_list[i_alteration].reference_sequence,g_test.significant_alteration_list[i_alteration].alternative_sequence)
			alignments = [alignments[0]]

		compact_cigard, uncompact_cigar = compute_cigar_string(alignments)

		if len(compact_cigard) != 6:
			logger.critical("More than one alteration: for %s", compact_cigard)
			is_multi = True
		# continue

		alteration_type = compact_cigard[3]
		# print reference_graph.dbg.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']
		ref_path_list = reference_graph.dbg.node[curr_alteration.reference_path[0]]['ref_list']
		if len(ref_path_list) == 1:
			splicing_variant = ref_path_list.keys()[0]
		else:
			splicing_variant = reference_graph.ref

		position = ref_path_list[splicing_variant] + compact_cigard[0]
		base_position = ref_path_list[splicing_variant]

		if re.match("rs", splicing_variant):
			reference_sequence = ""
			# print reference_graph.nt_ref[splicing_variant]
			for i_pos in range(0, compact_cigard[2]):
				if position + i_pos not in reference_graph.nt_ref[splicing_variant]:
					logger.critical("%d not in nt_ref dict of %s", position + i_pos, splicing_variant)
					continue
				reference_sequence += reference_graph.nt_ref[splicing_variant][position + i_pos]
			splicing_variant = reference_graph.ref
		else:
			reference_sequence = ref_seq[compact_cigard[0]:compact_cigard[0] + compact_cigard[2]]

		alteration_description = {
			"splicing_variant": splicing_variant,
			"compact_cigar": "".join(map(str, compact_cigard)),
			"is_multi": is_multi,
			"uncompact_cigar": uncompact_cigar,
			"base_position": base_position,
			"start": position,
			"end": None,
			'alt_sequence': None,
			"n_alterations": len([x for x in uncompact_cigar if x != "M"]),
			"reference_sequence": reference_sequence,
			"alt_type": alteration_type,
			"pvalue_ratio": curr_alteration.pvalue_ratio,
			"alt_read_count": curr_alteration.alternative_read_count,
			"ref_read_count": curr_alteration.reference_read_count,
			"z_score": float(curr_alteration.zscore),
			"alignment": alignments[0]
		}

		if alteration_type == "X":
			# c.76A>T
			if compact_cigard[2] != 1:
				logger.critical("More than one alteration: for %s", compact_cigard)
				alteration_description['is_multi'] = True
			# continue
			alteration_description['end'] = alteration_description['start']
			alteration_description['alt_sequence'] = alt_seq[compact_cigard[0]:compact_cigard[0] + compact_cigard[2]]
		elif alteration_type == "D":
			# c.76_78delACT
			alteration_description['end'] = position + len(reference_sequence) - 1
		else:
			# c.76_77insG
			alteration_description['alt_sequence'] = alt_seq[compact_cigard[0]:compact_cigard[0] + compact_cigard[2]]
			alteration_description['end'] = position + 1

		annotated_alterations.append(alteration_description)
	annotated_alterations = merge_identical_alterations(annotated_alterations)
	for alt in annotated_alterations:
		print_alteration(alt)
	return annotated_alterations


def print_alteration(alteration_description):
	alt_type = alteration_description['alt_type']
	alt_sequence = alteration_description['alt_sequence']
	reference_sequence = alteration_description['reference_sequence']
	ref_read_count = alteration_description['ref_read_count']
	splicing_variant = alteration_description['splicing_variant']
	start = alteration_description['start']
	alt_read_count = alteration_description['alt_read_count']
	p_value = alteration_description['pvalue_ratio']
	z_score = alteration_description['z_score']
	end = alteration_description['end']

	if alt_type == "X":
		print "%s:c.%d%s>%s\t%d\t%d\t%f\t%f" % (
			splicing_variant, start, reference_sequence, alt_sequence, ref_read_count,
			alt_read_count, p_value, z_score)
	elif alt_type == "D":
		print "%s:c.%d_%ddel%s\t%d\t%d\t%f\t%f" % (
			splicing_variant, start, end, reference_sequence, ref_read_count,
			alt_read_count, p_value, z_score)
	elif alt_type == "I":
		print "%s:c.%d_%dins%s\t%d\t%d\t%f\t%f" % (
			splicing_variant, start, end, alt_sequence, ref_read_count,
			alt_read_count, p_value, z_score)


def compute_cigar_string(alignments):
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
	return compact_cigard, uncompact_cigar
