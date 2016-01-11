# For one PCR amplicon
# !/usr/bin/env python
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

def build_read_library(FASTQFILE_PATH):
	read_library = collections.defaultdict(list)
	FASTQFILE_ALL = os.listdir(FASTQFILE_PATH)
	logger.info("Found %d fastq file to process", len(FASTQFILE_ALL))
	for j, a_fastq_file in time_iterator(FASTQFILE_ALL, logger, msg_prefix="Building read library"):
		if a_fastq_file == ".DS_Store":
			continue
		fastq = open(FASTQFILE_PATH + "/" + a_fastq_file, 'r')
		lines = fastq.readlines()
		fastq.close()
		lines = map(str.strip, lines)
		for i_line in range(1, len(lines), 4):
			read_library[a_fastq_file].append(lines[i_line])
	return read_library


def build_serialize_library(experiment_name,FASTQFILE_PATH):
	logger.info("will rebuild library")
	read_library = build_read_library(FASTQFILE_PATH)
	logger.info("Will save %d items", len(read_library))
	packed_docs = msgpack.packb(read_library, default=lambda x: x.__dict__)
	logger.info("Packed to %d chars", len(packed_docs))
	get_or_create_dir("data")
	get_or_create_dir("data/seq")
	tgt_file = "data/seq/" + experiment_name + "_%s_%d.packb" % ((int(time.time())), len(read_library))
	with open(tgt_file, "w") as f:
		f.write(packed_docs)

	logger.info("Serialized to file %s" % tgt_file)


def library_itit(experiment_name):
	global read_library
	avail_files = {x: x.split("_") for x in glob.glob("data/seq/" + experiment_name + "_*_*.packb")}.items()
	FASTQFILE_PATH = "data/fastq/"+experiment_name
	if len(avail_files) < 1:
		logger.info("Force rebuilding")
		build_serialize_library(experiment_name,FASTQFILE_PATH)
		avail_files = {x: x.split("_") for x in glob.glob("data/seq/" + experiment_name + "_*_*.packb")}.items()
	avail_files.sort(key=lambda x: float(x[1][1]), reverse=True)
	most_recent = avail_files[0][0]

	logger.info("will unpack read library ")
	with open(most_recent, "r") as f:
		read_library = msgpack.unpack(f)
	# un_packed_idx = [(ObjectId(x[0]), x[1]) for x in un_packed_idx]
	logger.info("De-Serialized %d read lib" % (len(read_library)))

# Sampling fonction from a coverage dict 
def sampling(coverage):
	coverage = 1000
	read_sampling = []
	read_sampling_for_remainder = []
	read_number_by_sample = coverage / len(read_library)
	read_number_remainder = coverage % len(read_library)
	for sample in read_library:
		read_sampling += random.sample(read_library[sample], read_number_by_sample)
		read_sampling_for_remainder += random.sample(read_library[sample], 1)
	read_sampling += random.sample(read_sampling_for_remainder, read_number_remainder)
	return read_sampling
