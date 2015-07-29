#!/usr/bin/env python
# coding=utf8
import collections

import os
import re
import random
import msgpack
import time
import glob
import sys
from helpers.helpers import time_iterator, get_or_create_dir
from helpers.logger import init_logger

logger = init_logger('SEQLIB')
logger.info("Setting up SEQLIB")

experiment_name = "TP53"

def build_read_library(FASTQFILE_PATH):
	pattern = re.compile('([NC])_(\d+)_(\d+)')
	read_library = {'N': collections.defaultdict(set), 'C': collections.defaultdict(set)}
	FASTQFILE_ALL = os.listdir(FASTQFILE_PATH)
	logger.info("Found %d fastq file to process", len(FASTQFILE_ALL))
	for j, a_fastq_file in time_iterator(FASTQFILE_ALL, logger, msg_prefix="Building read library"):
		if a_fastq_file == ".DS_Store":
			continue
		fragment = pattern.search(a_fastq_file).group(1)
		individu = pattern.search(a_fastq_file).group(2)
		fastq = open(FASTQFILE_PATH + "/" + a_fastq_file, 'r')
		lines = fastq.readlines()
		fastq.close()
		lines = map(str.strip, lines)
		# if individu not in read_library[fragment]:
		# read_library[fragment][individu] = []
		for i_line in range(1, len(lines), 4):
			read_library[fragment][individu].add(lines[i_line])
	#mutate everything back to lists
	read_library['N']={k:list(v) for k,v in read_library['N'].items()}
	read_library['C']={k:list(v) for k,v in read_library['C'].items()}
	return read_library


def build_serialize_library(FASTQFILE_PATH):
	logger.info("will rebuild library")
	read_library = build_read_library(FASTQFILE_PATH)
	logger.info("Will save %d items", len(read_library['N'])+len(read_library['C']))
	packed_docs = msgpack.packb(read_library, default=lambda x: x.__dict__)
	logger.info("Packed to %d chars", len(packed_docs))
	get_or_create_dir("data")
	get_or_create_dir("data/seq")
	tgt_file = "data/seq/"+experiment_name+"_%s_%d.packb" % ((int(time.time())), len(read_library))
	with open(tgt_file, "w") as f:
		f.write(packed_docs)

	logger.info("Serialized to file %s" % tgt_file)

# def library_itit(experiment_name):
FASTQFILE_PATH = "data/fastq/"+experiment_name
if 'rebuild_library' in sys.argv:
	build_serialize_library(FASTQFILE_PATH)
else:
	avail_files = {x: x.split("_") for x in glob.glob("data/seq/"+experiment_name+"_*_*.packb")}.items()
	if len(avail_files)<1:
		logger.info("Force rebuilding")
		build_serialize_library(FASTQFILE_PATH)
		avail_files = {x: x.split("_") for x in glob.glob("data/seq/"+experiment_name+"_*_*.packb")}.items()
	avail_files.sort(key=lambda x: float(x[1][1]), reverse=True)
	most_recent = avail_files[0][0]

	logger.info("will unpack read library ")
	with open(most_recent, "r") as f:
		read_library = msgpack.unpack(f)
	# un_packed_idx = [(ObjectId(x[0]), x[1]) for x in un_packed_idx]
	logger.info("De-Serialized %d read lib" % (len(read_library['N'])+len(read_library['C'])))  


# Sampling fonction from a coverage dict (with keys 'N' et 'C')
def sampling(coverage_dict):
	read_sampling = []
	for fragment, coverage in coverage_dict.items():
		if fragment == "total":
			continue
		fragment=fragment[0]
		read_sampling_for_remainder = []
		read_number_by_sample = coverage / len(read_library[fragment])
		read_number_remainder = coverage % len(read_library[fragment])
		for sample in read_library[fragment]:
			read_sampling += random.sample(read_library[fragment][sample], read_number_by_sample)
			read_sampling_for_remainder += random.sample(read_library[fragment][sample], 1)
		read_sampling += random.sample(read_sampling_for_remainder, read_number_remainder)
	return read_sampling
