from argparse import ArgumentParser
import os
import pprint as pp
import collections
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import Levenshtein
import sys
from helpers.logger import init_logger

__author__ = 'hayssam'

logger = init_logger('Reads pre-processor')


def identify_primer_with_mismatches(seq, fw, rev, max_mismatch=8):
	for primer in fw:
		d = Levenshtein.distance(seq[:len(primer)], primer)
		if d < max_mismatch:
			return +1, primer
	for primer in rev:
		d = Levenshtein.distance(seq[:len(primer)], primer)
		if d < max_mismatch:
			return -1, primer
	return None, None


def identify_primer(seq, fw, rev):
	for primer in fw:
		if seq.startswith(primer):
			return +1, primer
	for primer in rev:
		if seq.startswith(primer):
			return -1, primer
	return None, None


def process_file(input_file, forward_primers, reverse_primers):
	# parse input FASTQ file
	with open(input_file) as f:
		sequence_iterator = SeqIO.parse(f, format="fastq", alphabet=generic_dna)
		# TODO perform lazy evaluation ?
		sequence_records = list(sequence_iterator)
	logger.info("Loaded %d sequences", len(sequence_records))
	# determine first nt
	primer_to_seq = collections.defaultdict(set)
	seq_without_primers = set()
	first_kmer_dict = collections.Counter()
	first_kmer_length = max(map(len, forward_primers + reverse_primers))
	for record in sequence_records:
		sequence = str(record.seq)
		orientation, primer = identify_primer(sequence, forward_primers, reverse_primers)
		if primer:
			primer_to_seq[(primer, orientation)].add(record)
			continue
		# Try the RC of the seq
		rc_sequence = str(record.reverse_complement().seq)
		orientation, primer = identify_primer(rc_sequence, forward_primers, reverse_primers)
		if primer:
			primer_to_seq[(primer, orientation * -1)].add(record)
			continue
		# try with mismatch
		orientation, primer = identify_primer_with_mismatches(seq=sequence, fw=forward_primers, rev=reverse_primers, max_mismatch=8)
		if primer:
			primer_to_seq[(primer + "_X", orientation)].add(record)
			continue
		orientation, primer = identify_primer_with_mismatches(seq=rc_sequence, fw=forward_primers, rev=reverse_primers, max_mismatch=8)
		if primer:
			primer_to_seq[(primer + "_X", orientation * -1)].add(record)
			continue
		# no primer found
		seq_without_primers.add(record)
		first_kmer_dict[str(record[:first_kmer_length].seq)] += 1

	logger.info("Primer mapping stats: %d without identified primers, %s", len(seq_without_primers),
				{k: len(v) for k, v in primer_to_seq.items()})
	logger.info("Not mappable most frequent:\n%s", pp.pformat(first_kmer_dict.most_common(10)))
	# perform the RC when needed and output
	output_records = []
	for (primer, orientation), record_list in primer_to_seq.items():
		if orientation == -1:
			output_records.extend(map(lambda x: reverse_record(x), record_list))
		else:
			output_records.extend(record_list)
	return output_records


def reverse_record(x):
	rec = x.reverse_complement()
	rec.id = x.id + "_rev"
	rec.description = "REV " + x.description
	return rec


def validate_primers(args):
	# global forward_primers, reverse_primers
	forward_primers = [x.strip() for x in args.forward_primers.split(",")]
	reverse_primers = [x.strip() for x in args.reverse_primers.split(",")]
	assert len([x for x in forward_primers if set(x) != {"A", "T", "C", "G"}]) == 0, "Incorrect nucleotide in primer list %s" % (
		forward_primers)
	assert len([x for x in reverse_primers if set(x) != {"A", "T", "C", "G"}]) == 0, "Incorrect nucleotide in primer list %s" % (
		reverse_primers)
	primer_length = map(len, forward_primers + reverse_primers)
	logger.info("Identifed %d primers of len between %d and %d", len(forward_primers + reverse_primers), min(primer_length),
				max(primer_length))
	return forward_primers, reverse_primers


if __name__ == '__main__':
	parser = ArgumentParser()
	parser.add_argument('--fastq', help='FASTQ file to process', required=True, type=str)
	parser.add_argument('--output_suffix', help='Suffix extension', required=True, type=str)
	parser.add_argument('--output_prefix', help='Prefix extension', required=True, type=str)
	parser.add_argument('--forward_primers', help='Comma separated list of forward primers', required=True, type=str)
	parser.add_argument('--reverse_primers', help='Comma separated list of reverse primers', required=True, type=str)
	args = parser.parse_args()

	forward_primers, reverse_primers = validate_primers(args)
	processed_reads = process_file(input_file=args.fastq, forward_primers=forward_primers, reverse_primers=reverse_primers)
	ext = os.path.splitext(args.fastq)[1]
	output_file_name = args.output_prefix + "/" + os.path.basename(args.fastq) + "_" + args.output_suffix + ext
	with open(output_file_name, "w") as f:
	 	SeqIO.write(processed_reads, f, "fastq")
	sys.exit(0)
