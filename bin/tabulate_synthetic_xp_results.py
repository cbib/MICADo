import itertools

__author__ = 'hayssam'
import difflib
import os
import re
import simplejson as json
from helpers.helpers import time_iterator
from helpers.logger import init_logger

__author__ = 'hayssam'
import pandas as pd

logger = init_logger("GATKPPROC", {})


def hash_dict(d):
	return hash(tuple(sorted(d.items())))


def enrich_caller_record(a_rec, caller):
	res_dict = {'alt_sequence': a_rec['ALT'], 'ref_sequence': a_rec['REF']}
	res_dict['alt_length'] = abs(len(res_dict['alt_sequence']) - len(res_dict['ref_sequence']))

	matcher = difflib.SequenceMatcher(a=res_dict['ref_sequence'], b=res_dict['alt_sequence'])
	op = [x[0] for x in matcher.get_opcodes() if x[0] != 'equal']
	if "insert" in op:
		res_dict['alt_type'] = "I"
	elif "delete" in op:
		res_dict['alt_type'] = "D"
	elif "replace" in op:
		res_dict['alt_type'] = "X"

	res_dict['start'] = int(a_rec['POS'])
	res_dict['origin'] = caller
	res_dict['hash'] = hash_dict(res_dict)
	return res_dict


def parse_run_log(a_run_log):
	if not os.path.isfile(a_run_log):
		logger.critical("No run log for %s", a_run_log)
		return [-1], [-1]
	log_content = open(a_run_log).readlines()
	run_time_re = re.compile(r"\s+([0-9\.]+) real")
	run_times = [run_time_re.findall(line) for line in log_content]
	run_times = [float(x[0]) for x in run_times if len(x) > 0]
	memory_usage_re = re.compile(r"([0-9]+)  maximum resident set size")
	memory_usage = [memory_usage_re.findall(line) for line in log_content]
	memory_usage = [float(x[0]) for x in memory_usage if len(x) > 0]

	if len(memory_usage) < 1:
		memory_usage = [-1]
	if len(run_times) < 1:
		run_times = [-1]

	return run_times, memory_usage


def process_varscan_sample(sample_name):
	# logger.info("Processing sample %s", sample_name)
	a_sample = "data/synthetic/results/varscan/%s_on_NM_000546_5.vcf" % sample_name
	keys = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
	content = pd.DataFrame.from_records([dict(zip(keys, x.split("\t"))) for x in open(a_sample).readlines() if not x.startswith("#")])

	altered = content
	result_dict = {}
	if len(altered) < 1:
		result_dict['varscan'] = []
	else:
		enriched_varscan = list(altered.apply(lambda x: enrich_caller_record(x, "varscan"), axis=1))
		result_dict['varscan'] = enriched_varscan
	result_dict['sample_name'] = sample_name
	# Add runtime info
	run_log = "exec_logs/varscan_%s_on_NM_000546_5.txt" % sample_name
	mapper_log = "exec_logs/gmap_log_%s.txt" % sample_name

	run_time, memory = parse_run_log(run_log)
	run_time_gmap, memory_gmap = parse_run_log(mapper_log)

	result_dict['total_run_time'] = sum(run_time) + sum(run_time_gmap)
	result_dict['max_memory'] = max(memory + memory_gmap)

	return result_dict


def process_gatk_sample(sample_name):
	# logger.info("Processing sample %s", sample_name)
	a_sample = "data/synthetic/results/gatk/%s_on_NM_000546_5_raw.vcf" % sample_name

	alteration_description = json.loads(open("data/synthetic/results/sampler/%s.alterations.json" % sample_name).read())
	keys = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
	content = pd.DataFrame.from_records([dict(zip(keys, x.split("\t"))) for x in open(a_sample).readlines() if not x.startswith("#")])

	if len(content) == 0:
		altered = []
	else:
		altered = content.query("ALT!='<NON_REF>'")
		altered['ALT'] = [a_rec['ALT'].split(",")[0] for _, a_rec in altered.iterrows()]
	result_dict = {}
	result_dict.update(alteration_description)
	if len(altered) < 1:
		result_dict['gatk'] = []
	else:
		enriched_gatk = list(altered.apply(lambda x: enrich_caller_record(x, "gatk"), axis=1))
		result_dict['gatk'] = enriched_gatk

	result_dict['sample_name'] = sample_name
	run_log = "exec_logs/gatk_%s_on_NM_000546_5.txt" % sample_name
	mapper_log = "exec_logs/gmap_log_%s.txt" % sample_name
	run_time, memory = parse_run_log(run_log)
	run_time_gmap, memory_gmap = parse_run_log(mapper_log)

	result_dict['total_run_time'] = sum(run_time) + sum(run_time_gmap)
	result_dict['max_memory'] = max(memory + memory_gmap)

	return result_dict


def process_micado_sample(sample_name):
	# micado_content = json.load(open("../micado_synthetic_results/synthetic/%s.combined_alterations.json" % sample_name))
	micado_content = json.load(open("data/synthetic/results/micado/%s.significant_alterations.json" % sample_name))
	alteration_description = json.loads(open("data/synthetic/results/sampler/%s.alterations.json" % sample_name).read())
	run_log = "exec_logs/micado_log_%s.txt" % sample_name
	altered = micado_content['significant_alterations']
	result_dict = {}
	result_dict.update(alteration_description)
	if len(altered) < 1:
		result_dict['micado'] = []
	else:
		for x in altered:
			x['origin'] = 'micado'

		result_dict['micado'] = altered

	result_dict['sample_name'] = sample_name
	run_time, memory = parse_run_log(run_log)
	result_dict['total_run_time'] = sum(run_time)
	result_dict['max_memory'] = max(memory)
	return result_dict


def flatten_sample(result_dict, caller_key):
	base_dict = {
		"git_revision_hash": result_dict['sampler']['git_revision_hash'],
		"sample": result_dict['sample_name'],
		"n_reads": result_dict['sampler']['parameters']['n_reads'],
		"n_caller_alterations": len(result_dict[caller_key]),
		'total_run_time': result_dict['total_run_time'],
		'max_memory': result_dict['max_memory'],
		"caller": caller_key,
	}
	for x in result_dict[caller_key]:
		alt_dict = x
		alt_dict.update(base_dict)
		yield alt_dict

pd.set_option('display.width', 250)


def process_gatk_samples():
	avail_samples = [x for x in os.listdir("data/synthetic/results/gatk/") if x.endswith(".vcf")]
	for _, samp in time_iterator(avail_samples, logger, msg_prefix="GATK results"):
		name = samp.split("_on_")[0]
		try:
			res = process_gatk_sample(sample_name=name)
		except Exception as e:
			print "failed on sample", samp
			# raise e
			continue
		yield res


def process_varscan_samples():
	avail_samples = [x for x in os.listdir("data/synthetic/results/varscan/") if x.endswith(".vcf")]
	for _, samp in time_iterator(avail_samples, logger, msg_prefix="VARSCAN results"):
		name = samp.split("_on_")[0]
		try:
			res = process_varscan_sample(sample_name=name)
		except Exception as e:
			print "failed on sample", samp
			continue
		# raise e
		yield res


def process_micado_samples():
	# avail_samples = [x for x in os.listdir("../micado_synthetic_results/synthetic/") if x.endswith(".significant_alterations.json")]
	avail_samples = [x for x in os.listdir("data/synthetic/results/micado/") if x.endswith(".significant_alterations.json")]
	for _, samp in time_iterator(avail_samples, logger, msg_prefix="MICADo results"):
		name = samp.split(".")[0]
		try:
			res = process_micado_sample(sample_name=name)
		except Exception as e:
			print "failed on sample", samp
			# raise e
			continue
		yield res


if __name__ == '__main__':
	gatk_aggregated_results = pd.DataFrame.from_records(itertools.chain(*[flatten_sample(x, "gatk") for x in process_gatk_samples()]))
	gatk_aggregated_results.to_csv("data/synthetic/summary/gatk_results_on_synthetic_data.csv")

	varscan_aggregated_results = pd.DataFrame.from_records(itertools.chain(*[flatten_sample(x, "varscan") for x in process_varscan_samples()]))
	varscan_aggregated_results.to_csv("data/synthetic/summary/varscan_results_on_synthetic_data.csv")

	micado_aggregated_results = pd.DataFrame.from_records(itertools.chain(*[flatten_sample(x, "micado") for x in process_micado_samples()]))
	micado_aggregated_results.to_csv("data/synthetic/summary/micado_results_on_synthetic_data.csv")
