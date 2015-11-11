# coding=utf-8
from argparse import ArgumentParser
import collections
import pprint as pp
import json
import random
import itertools
import pandas as pd
import time
import sys
from pyparsing import ParseException
from helpers import helpers
from helpers.logger import init_logger
from helpers.helpers import time_iterator, get_git_revision_hash
from read_sampler.cigar_parser import parse_cigar_string

pd.set_option('display.width', 250)

__author__ = 'hayssam'

logger = init_logger(name="READSAMPLER")


# convert coordinates
def coordinate_map(an_alignment_row):
	global args
	range_accumulator = []
	cigar = an_alignment_row.CIGAR
	start_pos = an_alignment_row.POS + args.systematic_offset
	label = an_alignment_row.QNAME
	read_start_range = (0, 0)
	ref_start_range = (start_pos, start_pos)
	read_last_range = read_start_range
	ref_last_range = ref_start_range
	try:
		parse_results = parse_cigar_string(cigar)
	except ParseException:
		logger.critical("Parse error for cigar string %s", cigar)
		return []
	for length, type in parse_results:
		if type == "M":
			read_current_range = (read_last_range[1], read_last_range[1] + length)
			ref_current_range = (ref_last_range[1], ref_last_range[1] + length)
		elif type == "I":
			read_current_range = (read_last_range[1], read_last_range[1] + length)
			ref_current_range = (ref_last_range[1], ref_last_range[1])
		elif type == "D":
			read_current_range = (read_last_range[1], read_last_range[1])
			ref_current_range = (ref_last_range[1], ref_last_range[1] + length)
		elif type == "S":
			read_current_range = (read_last_range[1], read_last_range[1] + length)
			ref_current_range = (ref_last_range[1], ref_last_range[1])
		range_accumulator.append({"label": label, "type": type, "length": length, "ref_coord": ref_current_range, "read_coord": read_current_range})
		read_last_range = read_current_range
		ref_last_range = ref_current_range
	return range_accumulator


def random_alteration(start, end, weights, multi_mismatch=False):
	# sample a position, alt type, length and content
	a_pos = random.randint(start, end - MAX_LEN)
	print "POS:",a_pos
	a_type = helpers.weighted_choice(zip("IMD", weights))

	if a_type == "M" and not multi_mismatch:
		a_length = 1
	else:
		a_length = random.randint(1, 5)
	a_content = "".join([random.choice(random_nt) for _ in range(a_length)])
	a_qual = "q" * a_length
	# TODO parameter to force mismatches to be of length 1 or 2 (more realistic)
	if a_type == "D":
		return (a_pos, a_pos + a_length - 1), (a_type, None, None)
	elif a_type == "I":
		return (a_pos, a_pos), (a_type, a_content, a_qual)
	else:
		return (a_pos, a_pos + a_length - 1), (a_type, a_content, a_qual)


# find all regions containing the sampled pos
def region_overlap(x, y):
	return x[0] <= y[1] and y[0] <= x[1]


def region_overlap_right_strict(x, y):
	return x[0] <= y[1] and y[0] < x[1]


def transform_coordinate(ref_range, read_range, pos):
	return read_range[0] + (pos[0] - ref_range[0]), read_range[0] + (pos[1] - ref_range[0])


def mutating_sequence_iterator(read_label, alterations=None, output="seq"):
	"""
	Workhorse to alter a read. Given an alteration dict in genomic coordinates, translate it into reads coordinates.
	Then (lazily) iterate other all bases from the read sequence, and yield base after base. Yielded bases can be altered, subistitued or skipped.

	:param read_label: string: a read label from the index
	:param alterations: dict: an alteration dict, as generated from random_alteration
	:param output: string: kind of output
	:raise StopIteration: at the end of the read
	"""
	if not alterations:
		alterations = {}

	global all_ranges, aligned_reads

	# Switch output mode
	if output == "seq":
		read_sequence = aligned_reads.ix[read_label].SEQ
	else:
		read_sequence = aligned_reads.ix[read_label].QUAL

	# coordinate mapping
	read_alterations = map_alterations_to_read_coordinates(alterations, read_label)
	read_length = len(read_sequence)

	# iterator over read positions
	current_read_position = 0
	while current_read_position < read_length:
		base_is_covered = False

		for read_coord, (alt_type, alt, qual) in read_alterations.items():

			if (read_coord[0] <= current_read_position) and (current_read_position <= read_coord[1]):
				# Switch output mode
				# TODO should yield both quality and base in a single call => 2x speedup
				if output == "qual":
					alt_seq = qual
				else:
					alt_seq = alt
				if alt_type == "M":
					base_is_covered = True
					for alt_i, i in enumerate(range(read_coord[0], read_coord[1] + 1)):
						yield alt_seq[alt_i]
						current_read_position += 1
				elif alt_type == "I":
					for base in alt_seq:
						yield base
				elif alt_type == "D":
					base_is_covered = True
					for alt_i, i in enumerate(range(read_coord[0], read_coord[1] + 1)):
						current_read_position += 1
				break
		if base_is_covered:
			continue
		yield read_sequence[current_read_position]
		current_read_position += 1
	raise StopIteration


def map_alterations_to_read_coordinates(alterations, read_label):
	global all_ranges
	coord_map = all_ranges.ix[read_label]
	if len(coord_map.shape) == 1:  # we only had one row, make it a DF
		coord_map = pd.DataFrame([coord_map])
	read_alterations = {}
	for alt_ref_coord, alt in alterations.items():
		affected_regions = coord_map[coord_map.ref_coord.apply(lambda r: region_overlap_right_strict(r, alt_ref_coord))].query("type=='M'")
		# translate alteration map to read coordinate
		for i, region in affected_regions.iterrows():
			read_coordinates = transform_coordinate(region.ref_coord, region.read_coord, alt_ref_coord)
			read_alterations[read_coordinates] = alt
	return read_alterations


def clean_label(lbl):
	"""
	Stupid cleaner for STAR aligner who disallows "/"
	:type lbl: string
	"""
	return lbl.replace("/", "_")


def min_dist(int_list):
	_min_dist = int_list[0]
	for pair in itertools.izip(int_list, int_list[1:]):
		this_dist = abs(pair[0] - pair[1])
		if this_dist < _min_dist:
			_min_dist = abs(pair[0] - pair[1])
	return _min_dist


def build_a_sample(n_reads, fraction_altered, n_alterations, output_file_prefix, alterations_weight=None, multi_mismatch=False):
	if not alterations_weight:
		alterations_weight = [1.0, 1.0, 1.0]
	global all_ranges

	# sample some reads
	sub_reads = aligned_reads.sample(n=n_reads, random_state=args.seed, replace=False)

	# compute reference coordinates using the CIGAR
	all_ranges = []
	for i, an_alignment in sub_reads.iterrows():
		all_ranges.extend(coordinate_map(an_alignment))
	logger.info("Mapped coordinates to reference")
	all_ranges = pd.DataFrame.from_records(all_ranges)
	all_ranges.set_index("label", inplace=True, drop=False)

	# sample altered reads
	altered_reads_labels = sub_reads.QNAME.sample(int(len(sub_reads) * fraction_altered), random_state=args.seed, replace=False)
	altered_read_rows = all_ranges.ix[altered_reads_labels]
	non_altered_reads_labels = set(sub_reads.QNAME).difference(altered_reads_labels)
	assert set(altered_reads_labels).isdisjoint(set(non_altered_reads_labels))

	# pick a random label to test alterations
	a_label = random.choice(altered_reads_labels)

	some_alterations = generate_alterations(a_label, alterations_weight, altered_read_rows, multi_mismatch, n_alterations, sub_reads)

	if args.do_not_output_reads:
		return some_alterations

	# generate original reads
	with open(output_file_prefix + "_non_alt.fastq", "w") as f:
		for i, read_label in time_iterator(sub_reads.QNAME, logger, msg_prefix="Generating non altered fastq, non altered reads", delta_percent=0.1):
			print >> f, "@%s" % (clean_label(read_label)) + "_ORIG"
			print >> f, sub_reads.ix[read_label].SEQ
			print >> f, "+"
			print >> f, sub_reads.ix[read_label].QUAL
			# print >> f, "\n"

	# generate altered reads fastq files
	output_reads = set()
	with open(output_file_prefix + ".fastq", "w") as f:

		for i, read_label in time_iterator(altered_reads_labels, logger, msg_prefix="Generating altered fastq, altered reads", delta_percent=0.1):
			assert read_label not in output_reads
			output_reads.add(read_label)
			print >> f, "@%s" % (clean_label(read_label)) + "_ALT"
			print >> f, "".join(mutating_sequence_iterator(read_label=read_label, alterations=some_alterations))
			print >> f, "+"
			print >> f, "".join(mutating_sequence_iterator(read_label=read_label, alterations=some_alterations, output="qual"))
			# print >> f, "\n"

		for i, read_label in time_iterator(non_altered_reads_labels, logger, msg_prefix="Generating altered fastq, non altered reads", delta_percent=0.1):
			assert read_label not in output_reads
			output_reads.add(read_label)

			print >> f, "@%s" % (clean_label(read_label)) + "_ORIG"
			print >> f, sub_reads.ix[read_label].SEQ
			print >> f, "+"
			print >> f, sub_reads.ix[read_label].QUAL
			# print >> f, "\n"
	serialize_results(output_file_prefix, some_alterations)

	logger.info("finished generation for %d reads, %d alterations, output files are", n_reads, n_alterations)
	logger.info("%s: Original sampled reads", output_file_prefix + "_non_alt.fastq")
	logger.info("%s: Altered sampled reads", output_file_prefix + ".fastq")
	logger.info("%s: Alterations description", output_file_prefix + ".alterations.txt")
	logger.info("Alterations are %s", some_alterations)


def generate_alterations(a_label, alterations_weight, altered_reads_row, multi_mismatch, n_alterations, sub_reads):
	some_alterations = []
	# identify start and stop positions of reads that should be altered (with 10nt slack...)
	ref_start = min([min(x) for x in altered_reads_row.ref_coord]) + 10
	ref_end = max([max(x) for x in altered_reads_row.ref_coord]) - 10
	logger.info("Will sample between %d and %d ",ref_start,ref_end)
	# sample random alterations, reads that should be altered / kept as it
	alterations_modify_content = False
	max_try = 100
	i = 0
	while (not alterations_modify_content) and (i < max_try):
		some_alterations = dict([random_alteration(ref_start, ref_end, weights=alterations_weight, multi_mismatch=multi_mismatch) for _ in range(n_alterations)])
		# check that artificial alterations actually modify reads (case of generating a substitution corresponding to the actual content of the read)
		altered_sequence = "".join(mutating_sequence_iterator(read_label=a_label, alterations=some_alterations))
		non_altered_sequence = sub_reads.ix[a_label].SEQ
		if min_dist([x[0] for x in some_alterations]) <= 20:
			logger.info("Alterations %s are too close, iterating", some_alterations)
		elif altered_sequence != non_altered_sequence:
			alterations_modify_content = True
		else:
			logger.info("Alterations %s correspond to the real read content, iterating", some_alterations)
		i += 1
	logger.info("Generated alterations %s after %d trial", some_alterations, i)
	return some_alterations


def serialize_results(output_file_prefix, some_alterations):
	# prepare alterations for output
	PROGRAMEND = time.time()
	sampling_experiment = {}
	described_alterations = []
	for (start, end), (alt_type, alt_content, qual_content) in some_alterations.items():
		described_alterations.append({
			"start": start,
			'end': end,
			'alt_type': alt_type if alt_type != "M" else "X",  # for micado compatibility
			"alt_content": alt_content,
			"qual_string": qual_content,
			"alt_length": len(alt_content) if alt_content else 0
		})
	sampling_experiment['exec_time'] = PROGRAMEND - PROGRAMSTART
	sampling_experiment['parameters'] = vars(args)
	sampling_experiment['injected_alterations'] = described_alterations
	sampling_experiment['input_n_reads'] = len(aligned_reads)
	sampling_experiment['git_revision_hash'] = get_git_revision_hash()
	sampling_experiment['program_name'] = __file__
	with open(output_file_prefix + ".alterations.json", "w") as f:
		json.dump({"sampler": sampling_experiment}, f)


def parse_sam_file(sam_file):
	# sam_content = pd.DataFrame.from_csv(sam_file, sep="\t", header=5, index_col=None)
	sam_content = pd.read_csv(sam_file, sep="\t", header=5, index_col=None, quotechar=u"±")
	logger.info("Read %d lines SAM file", len(sam_content))
	sam_columns = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]
	this_columns = sam_columns + ["custom_%d" % i for i in range(len(sam_columns), sam_content.shape[1])]
	sam_content.columns = this_columns
	sam_content.set_index("QNAME", drop=False, inplace=True)
	return sam_content


if __name__ == '__main__':
	global args, random_nt
	parser = ArgumentParser()
	parser.add_argument('--n_reads', help='Number of reads for the generated sample', default=500, type=int, required=False)
	parser.add_argument('--fraction_altered', help='Fraction (between 0 and 1) of reads that should have recurrent alterations', required=False, type=float, default=0.1)
	parser.add_argument('--n_alterations', help='Number of alterations to insert', required=False, type=int, default=1)
	parser.add_argument('--output_file_prefix', help='output file prefix', required=True, type=str)
	parser.add_argument('--input_sam', help='Input SAM file', required=True, type=str)
	parser.add_argument('--alt_weight', help='Comma or "-" separated weights for alterations, in order:insertions, mismatches and deletions', required=False, type=str,
						default="1,1,1")
	parser.add_argument('--max_len', help='Max INDEL length', required=False, type=int, default=5)
	parser.add_argument('--seed', help='Random seed to use, by default: system time in seconds.', required=False, type=int, default=None)
	parser.add_argument('--multi_mismatch', help='Allow for mismatches over multiple nt', required=False, action='store_true')
	parser.add_argument('--output_lowercase', help='Output altered bases in lowercase', required=False, action='store_true')
	parser.add_argument('--systematic_offset', help='Systematically offset all reported reference coordinates by this value', required=False, type=int, default=0)
	parser.add_argument('--do_not_output_reads', help='Output sampled alteration information only', action="store_true")

	args = parser.parse_args()
	# starting_file = data_dir + "/alignments/C_model_GMAPno40_NM_000546.5.sam"
	PROGRAMSTART = time.time()
	if not args.seed:
		args.seed = int(time.time())
	else:
		args.seed = int(args.seed)

	logger.info("Random seed is %d[%s]", args.seed, type(args.seed))
	random.seed(args.seed)
	logger.info("Will parse input SAM")
	aligned_reads = parse_sam_file(args.input_sam)
	logger.info("Parsed %d reads from the input SAM file", len(aligned_reads))
	logger.info("Starting reads generation")
	MAX_LEN = args.max_len
	all_ranges = None

	if args.output_lowercase:
		random_nt = "actg"
	else:
		random_nt = "ACTG"
	# parse alteration weight
	try:
		if "-" in args.alt_weight:
			alt_weights = map(float, args.alt_weight.split("-"))
		else:
			alt_weights = map(float, args.alt_weight.split(","))
		assert len(alt_weights) == 3, "3 comma separated weights should be provided"
	except:
		raise SyntaxError
	sampled_alt = []
	alterations = build_a_sample(
		n_reads=args.n_reads, n_alterations=args.n_alterations,
		fraction_altered=args.fraction_altered,
		output_file_prefix=args.output_file_prefix,
		alterations_weight=alt_weights,
		multi_mismatch=args.multi_mismatch
	)
	# pp.pprint(sampled_alt)
