import difflib
import itertools
import numpy as np
import os
import re
import simplejson as json
from helpers.helpers import time_iterator
from helpers.logger import init_logger

__author__ = 'hayssam'
import pandas as pd



# mode="SUPERVISED"
# XPDIR = "data/synthetic/"
# XPKEY = "synthetic"

mode = "UNSUPERVISED"
XPDIR = "data/tp53_analysis/"
XPKEY = "pool0"


logger = init_logger("RESULTPROCESSING[%s]"%(mode), {})
def hash_dict(d):
	return hash(tuple(sorted(d.items())))


def closest_alteration(alt, alt_list):
	if len(alt_list) < 1:
		return None
	closest = min(alt_list, key=lambda x: abs(x['start'] - alt['start']))

	return closest


def is_match(x, closest):
	if not closest:
		return False
	return abs(x['start'] - closest['start']) <= 8 and x['alt_type'] == closest['alt_type']


def enrich_caller_record_unsupervised(a_rec, caller):
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


## Unsupervised (no real mutations are known) mode, all identified alterations are returned

def process_varscan_sample_unsupervised(sample_name):
	# logger.info("Processing sample %s", sample_name)
	a_sample = XPDIR + "results/varscan/%s_on_NM_000546_5.vcf" % sample_name
	keys = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
	content = pd.DataFrame.from_records([dict(zip(keys, x.split("\t"))) for x in open(a_sample).readlines() if not x.startswith("#")])

	altered = content
	result_dict = {}
	if len(altered) < 1:
		result_dict['varscan'] = []
	else:
		enriched_varscan = list(altered.apply(lambda x: enrich_caller_record_unsupervised(x, "varscan"), axis=1))
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


def process_gatk_sample_unsupervised(sample_name):
	# logger.info("Processing sample %s", sample_name)
	a_sample = XPDIR + "results/gatk/%s_on_NM_000546_5_raw.vcf" % sample_name

	# alteration_description = json.loads(open(XPDIR + "results/sampler/%s.alterations.json" % sample_name).read())
	keys = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
	content = pd.DataFrame.from_records([dict(zip(keys, x.split("\t"))) for x in open(a_sample).readlines() if not x.startswith("#")])

	if len(content) == 0:
		altered = []
	else:
		altered = content.query("ALT!='<NON_REF>'")
		altered['ALT'] = [a_rec['ALT'].split(",")[0] for _, a_rec in altered.iterrows()]
	result_dict = {}
	if len(altered) < 1:
		result_dict['gatk'] = []
	else:
		enriched_gatk = list(altered.apply(lambda x: enrich_caller_record_unsupervised(x, "gatk"), axis=1))
		result_dict['gatk'] = enriched_gatk

	result_dict['sample_name'] = sample_name
	run_log = "exec_logs/gatk_%s_on_NM_000546_5.txt" % sample_name
	mapper_log = "exec_logs/gmap_log_%s.txt" % sample_name
	run_time, memory = parse_run_log(run_log)
	run_time_gmap, memory_gmap = parse_run_log(mapper_log)

	result_dict['total_run_time'] = sum(run_time) + sum(run_time_gmap)
	result_dict['max_memory'] = max(memory + memory_gmap)

	return result_dict


def process_micado_sample_unsupervised(sample_name):
	# micado_content = json.load(open("../micado_synthetic_results/synthetic/%s.combined_alterations.json" % sample_name))
	micado_content = json.load(open(XPDIR + "results/micado/%s.significant_alterations.json" % sample_name))
	# alteration_description = json.loads(open(XPDIR + "results/sampler/%s.alterations.json" % sample_name).read())
	run_log = "exec_logs/micado_log_%s.txt" % sample_name
	altered = micado_content['significant_alterations']
	result_dict = {}
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
	print sample_name
	return result_dict


def flatten_sample_unsupervised(result_dict, caller_key):
	base_dict = {
		"sample": result_dict['sample_name'],
		"sample_key": "_".join(result_dict['sample_name'].split("_")[1:]),
		"sample_fragment": result_dict['sample_name'].split("_")[0],
		"n_caller_alterations": len(result_dict[caller_key]),
		'total_run_time': result_dict['total_run_time'],
		'max_memory': result_dict['max_memory'],
		"caller": caller_key,
	}
	if len(result_dict[caller_key])<1:
		yield base_dict
	for x in result_dict[caller_key]:
		alt_dict = x
		alt_dict.update(base_dict)
		yield alt_dict


## Supervised mode

# enrich




def enrich_caller_record(a_rec, alteration_description, caller):
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

	closest = closest_alteration(res_dict, alteration_description['sampler']['injected_alterations'])
	res_dict['is_match'] = is_match(res_dict, closest)
	return res_dict


def classification_report(result_dict, tgt):
	result_dict['tp'] = len([x for x in result_dict[tgt] if x['is_match']])
	result_dict['fp'] = len([x for x in result_dict[tgt] if not x['is_match']])
	result_dict['fn'] = len([x for x in result_dict['alterations'] if not x['is_matched']])

	if len(result_dict[tgt]) < 1:
		result_dict['precision'] = np.NaN
	else:
		result_dict['precision'] = result_dict['tp'] * 1.0 / len(result_dict[tgt])
	result_dict['recall'] = len([x for x in result_dict['alterations'] if x['is_matched']]) * 1.0 / len(result_dict['alterations'])
	return result_dict


def process_varscan_sample(sample_name):
	# logger.info("Processing sample %s", sample_name)
	a_sample = XPDIR + "results/varscan/%s_on_NM_000546_5.vcf" % sample_name

	alteration_description = json.loads(open(XPDIR + "results/sampler/%s.alterations.json" % sample_name).read())
	keys = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
	content = pd.DataFrame.from_records([dict(zip(keys, x.split("\t"))) for x in open(a_sample).readlines() if not x.startswith("#")])

	altered = content
	result_dict = {}
	result_dict.update(alteration_description)
	if len(altered) < 1:
		result_dict['varscan'] = []
	else:
		enriched_varscan = list(altered.apply(lambda x: enrich_caller_record(x, alteration_description, "varscan"), axis=1))
		result_dict['varscan'] = enriched_varscan

	# enriched injected
	enriched_sampler = []
	for alt in alteration_description['sampler']['injected_alterations']:
		closest = closest_alteration(alt, result_dict['varscan'])
		if is_match(alt, closest):
			alt['is_matched'] = True
		else:
			alt['is_matched'] = False
		enriched_sampler.append(alt)
	result_dict['alterations'] = enriched_sampler
	result_dict['sample_name'] = sample_name
	# Add runtime info
	run_log = "exec_logs/varscan_%s_on_NM_000546_5.txt" % sample_name
	mapper_log = "exec_logs/gmap_log_%s.txt" % sample_name

	run_time, memory = parse_run_log(run_log)
	run_time_gmap, memory_gmap = parse_run_log(mapper_log)

	result_dict['total_run_time'] = sum(run_time) + sum(run_time_gmap)
	result_dict['max_memory'] = max(memory + memory_gmap)

	return classification_report(result_dict, "varscan")


def process_gatk_sample(sample_name):
	# logger.info("Processing sample %s", sample_name)
	a_sample = XPDIR + "results/gatk/%s_on_NM_000546_5_raw.vcf" % sample_name

	alteration_description = json.loads(open(XPDIR + "results/sampler/%s.alterations.json" % sample_name).read())
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
		enriched_gatk = list(altered.apply(lambda x: enrich_caller_record(x, alteration_description, "gatk"), axis=1))
		result_dict['gatk'] = enriched_gatk

	# enriched injected
	enriched_sampler = []
	for alt in alteration_description['sampler']['injected_alterations']:
		closest = closest_alteration(alt, result_dict['gatk'])
		if is_match(alt, closest):
			alt['is_matched'] = True
		else:
			alt['is_matched'] = False
		enriched_sampler.append(alt)
	result_dict['alterations'] = enriched_sampler
	result_dict['sample_name'] = sample_name
	run_log = "exec_logs/gatk_%s_on_NM_000546_5.txt" % sample_name
	mapper_log = "exec_logs/gmap_log_%s.txt" % sample_name
	run_time, memory = parse_run_log(run_log)
	run_time_gmap, memory_gmap = parse_run_log(mapper_log)

	result_dict['total_run_time'] = sum(run_time) + sum(run_time_gmap)
	result_dict['max_memory'] = max(memory + memory_gmap)

	return classification_report(result_dict, "gatk")


def process_micado_sample(sample_name):
	# micado_content = json.load(open("../micado_synthetic_results/synthetic/%s.combined_alterations.json" % sample_name))
	micado_content = json.load(open(XPDIR + "results/micado/%s.significant_alterations.json" % sample_name))
	alteration_description = json.loads(open(XPDIR + "results/sampler/%s.alterations.json" % sample_name).read())
	run_log = "exec_logs/micado_log_%s.txt" % sample_name
	altered = micado_content['significant_alterations']
	result_dict = {}
	result_dict.update(alteration_description)
	if len(altered) < 1:
		result_dict['micado'] = []
	else:
		for x in altered:
			x['origin'] = 'micado'
			closest = closest_alteration(x, alteration_description['sampler']['injected_alterations'])
			x['is_match'] = is_match(x, closest)

		result_dict['micado'] = altered
	# enriched injected
	enriched_sampler = []
	for alt in alteration_description['sampler']['injected_alterations']:
		closest = closest_alteration(alt, result_dict['micado'])
		if is_match(alt, closest):
			alt['is_matched'] = True
		else:
			alt['is_matched'] = False
		enriched_sampler.append(alt)
	result_dict['alterations'] = enriched_sampler
	result_dict['sample_name'] = sample_name
	run_time, memory = parse_run_log(run_log)
	result_dict['total_run_time'] = sum(run_time)
	result_dict['max_memory'] = max(memory)
	return classification_report(result_dict, "micado")


def flatten_sample(result_dict, caller_key):
	return {
		"git_revision_hash": result_dict['sampler']['git_revision_hash'],
		"sample": result_dict['sample_name'],
		"n_reads": result_dict['sampler']['parameters']['n_reads'],
		"fraction_altered": result_dict['sampler']['parameters']['fraction_altered'],
		"n_alterations": result_dict['sampler']['parameters']['n_alterations'],
		"n_caller_alterations": len(result_dict[caller_key]),
		"precision": result_dict['precision'],
		'recall': result_dict['recall'],
		'total_run_time': result_dict['total_run_time'],
		'max_memory': result_dict['max_memory'],
		"tp": result_dict['tp'],
		"fp": result_dict['fp'],
		'fn': result_dict['fn'],
		"caller": caller_key,
	}


pd.set_option('display.width', 250)


def process_gatk_samples():
	avail_samples = [x for x in os.listdir(XPDIR + "results/gatk/") if x.endswith(".vcf")]
	for _, samp in time_iterator(avail_samples, logger, msg_prefix="GATK results"):
		name = samp.split("_on_")[0]
		try:
			if mode == "UNSUPERVISED":
				res = process_gatk_sample_unsupervised(sample_name=name)
			else:
				res = process_gatk_sample(sample_name=name)
		except Exception as e:
			print "failed on sample", samp
			raise e
		# continue
		yield res


def process_varscan_samples():
	avail_samples = [x for x in os.listdir(XPDIR + "results/varscan/") if x.endswith(".vcf")]
	for _, samp in time_iterator(avail_samples, logger, msg_prefix="VARSCAN results"):
		name = samp.split("_on_")[0]
		try:
			if mode == "UNSUPERVISED":
				res = process_varscan_sample_unsupervised(sample_name=name)
			else:
				res = process_varscan_sample(sample_name=name)
		except Exception as e:
			print "failed on sample", samp
			# continue
			raise e
		yield res


def process_micado_samples():
	# avail_samples = [x for x in os.listdir("../micado_synthetic_results/synthetic/") if x.endswith(".significant_alterations.json")]
	avail_samples = [x for x in os.listdir(XPDIR + "results/micado/") if x.endswith(".significant_alterations.json")]
	for _, samp in time_iterator(avail_samples, logger, msg_prefix="MICADo results"):
		name = samp.split(".")[0]
		try:
			if mode == "UNSUPERVISED":
				res = process_micado_sample_unsupervised(sample_name=name)
			else:
				res = process_micado_sample(sample_name=name)
		except Exception as e:
			print "failed on sample", samp
			raise e
		# continue
		yield res


if __name__ == '__main__':
	if mode == "UNSUPERVISED":
		gatk_aggregated_results = pd.DataFrame.from_records(
			itertools.chain(*[flatten_sample_unsupervised(x, "gatk") for x in process_gatk_samples()]))
		gatk_aggregated_results.to_csv(XPDIR + "summary/agg_unsupervised_gatk_results_on_%s_data.csv" % XPKEY)

		varscan_aggregated_results = pd.DataFrame.from_records(
			itertools.chain(*[flatten_sample_unsupervised(x, "varscan") for x in process_varscan_samples()]))
		varscan_aggregated_results.to_csv(XPDIR + "summary/agg_unsupervised_varscan_results_on_%s_data.csv" % XPKEY)

		micado_aggregated_results = pd.DataFrame.from_records(
			itertools.chain(*[flatten_sample_unsupervised(x, "micado") for x in process_micado_samples()]))
		micado_aggregated_results.to_csv(XPDIR + "summary/agg_unsupervised_micado_results_on_%s_data.csv" % XPKEY)


	else:
		gatk_aggregated_results = pd.DataFrame.from_records([flatten_sample(x, "gatk") for x in process_gatk_samples()])
		gatk_aggregated_results.to_csv(XPDIR + "summary/agg_gatk_results_on_%s_data.csv" % XPKEY)

		varscan_aggregated_results = pd.DataFrame.from_records([flatten_sample(x, "varscan") for x in process_varscan_samples()])
		varscan_aggregated_results.to_csv(XPDIR + "summary/agg_varscan_results_on_%s_data.csv" % XPKEY)

		micado_aggregated_results = pd.DataFrame.from_records([flatten_sample(x, "micado") for x in process_micado_samples()])
		micado_aggregated_results.to_csv(XPDIR + "summary/agg_micado_results_on_%s_data.csv" % XPKEY)
