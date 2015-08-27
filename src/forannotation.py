from helpers.helpers import time_iterator, get_or_create_dir
from helpers.logger import init_logger
from Bio import pairwise2
import re
logger = init_logger('Annotation')

def alteration_list_to_transcrit_mutation(g_test,g_ref):
	for i_alteration in range(0, len(g_test.significant_alteration_list)):
		alignments = pairwise2.align.globalms(g_test.significant_alteration_list[i_alteration].reference_sequence, g_test.significant_alteration_list[i_alteration].alternative_sequence, 2, -3, -5, -2)
		# if more than one alignment, choose the one which with the alteration at the left in the genome
		if len(alignments) > 1:
			# logger.critical("More than one alignment for %s vs %s", g_test.significant_alteration_list[i_alteration].reference_sequence,g_test.significant_alteration_list[i_alteration].alternative_sequence) 
			alignments = [alignments[0]]
		uncompact_cigar = ""
		compact_cigard = []
		for i_nucleotide in range(0,alignments[0][4]):
			if alignments[0][0][i_nucleotide] == alignments[0][1][i_nucleotide] :
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
		for i_nucleotide in range(0,len(uncompact_cigar)):
			if uncompact_cigar[i_nucleotide] != operation:
				compact_cigard += [count,operation]
				operation = uncompact_cigar[i_nucleotide]
				count = 0
			count += 1
		compact_cigard += [count,operation]
		# print compact_cigard
		if len(compact_cigard) == 6:
			alteration_type = compact_cigard[3]
			# print g_ref.dbg.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']
			if len(g_ref.dbg.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']) == 1:
				splicing_variant = g_ref.dbg.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list'].keys()[0]
			elif "NM_000546.5" in g_ref.dbg.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']:
				splicing_variant = "NM_000546.5"
			elif "NM_001126114.2" in g_ref.dbg.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']:
				splicing_variant = "NM_001126114.2"
			elif "NM_001126113.2" in g_ref.dbg.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']:
				splicing_variant = "NM_001126113.2"
			position = g_ref.dbg.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list'][splicing_variant]+compact_cigard[0]
			if re.match("rs",splicing_variant):
				reference = ""
				# print splicing_variant
				# print compact_cigard
				# print g_ref.nt_ref[splicing_variant]
				for i_pos in range(0,compact_cigard[2]):
					if position+i_pos not in g_ref.nt_ref[splicing_variant]:
						logger.critical("%d not in nt_ref dict of %s",position+i_pos,splicing_variant)
						continue
				 	reference += g_ref.nt_ref[splicing_variant][position+i_pos]
				splicing_variant = "NM_000546.5"
			else:
				reference = g_test.significant_alteration_list[i_alteration].reference_sequence[compact_cigard[0]:compact_cigard[0]+compact_cigard[2]]			
			if alteration_type == "X":
				# c.76A>T
				if compact_cigard[2] != 1:
					logger.critical("More than one alteration: for %s",compact_cigard)
					continue
				alteration = g_test.significant_alteration_list[i_alteration].alternative_sequence[compact_cigard[0]:compact_cigard[0]+1]
				print "%s:c.%d%s>%s\t%d\t%d\t%f\t%f"%(splicing_variant,position,reference,alteration,g_test.significant_alteration_list[i_alteration].reference_read_count,g_test.significant_alteration_list[i_alteration].alternative_read_count,g_test.significant_alteration_list[i_alteration].pvalue_ratio,float(g_test.significant_alteration_list[i_alteration].zscore))
			elif alteration_type == "D":
				# c.76_78delACT
				print "%s:c.%d_%ddel%s\t%d\t%d\t%f\t%f"%(splicing_variant,position,position+len(reference)-1,reference,g_test.significant_alteration_list[i_alteration].reference_read_count,g_test.significant_alteration_list[i_alteration].alternative_read_count,g_test.significant_alteration_list[i_alteration].pvalue_ratio,float(g_test.significant_alteration_list[i_alteration].zscore))
			else:
				# c.76_77insG
				alteration = g_test.significant_alteration_list[i_alteration].alternative_sequence[compact_cigard[0]:compact_cigard[0]+compact_cigard[2]]
				print "%s:c.%d_%dins%s\t%d\t%d\t%f\t%f"%(splicing_variant,position,position+1,alteration,g_test.significant_alteration_list[i_alteration].reference_read_count,g_test.significant_alteration_list[i_alteration].alternative_read_count,g_test.significant_alteration_list[i_alteration].pvalue_ratio,float(g_test.significant_alteration_list[i_alteration].zscore))
		else:
			logger.critical("More than one alteration: for %s",compact_cigard)

def splicing_variant_converter(spv):
	if spv not in ["NM_000546.5","NM_001126114.2","NM_001126113.2"]:
		spv = "NM_000546.5"
	return spv
