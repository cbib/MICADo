import difflib
import numpy as np
import os
import simplejson as json
from helpers.helpers import time_iterator
from helpers.logger import init_logger

__author__ = 'hayssam'
import pandas as pd

logger = init_logger("GATKPPROC", {})


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
	logger.info("Processing sample %s", sample_name)
	a_sample = "data/synthetic/results/varscan/%s_on_NM_000546_5.vcf" % sample_name
	alteration_description = json.loads(open("data/synthetic/results/sampler/%s.alterations.json" % sample_name).read())
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
	return classification_report(result_dict, "varscan")


def process_gatk_sample(sample_name):
	logger.info("Processing sample %s", sample_name)
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
	return classification_report(result_dict, "gatk")


def process_micado_sample(sample_name):
	micado_content = json.load(open("../micado_synthetic_results/synthetic/%s.significant_alterations.json" % sample_name))
	alteration_description = json.loads(open("data/synthetic/results/sampler/%s.alterations.json" % sample_name).read())
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
	return classification_report(result_dict, "micado")


def flatten_sample(result_dict, caller_key):
	return {
		"git_revision_hash": result_dict['sampler']['git_revision_hash'],
		"sample": result_dict['sample_name'],
		"n_reads": result_dict['sampler']['parameters']['n_reads'],
		"fraction_altered": result_dict['sampler']['parameters']['fraction_altered'],
		"n_alterations": result_dict['sampler']['parameters']['n_alterations'],
		"n_caller_alterations": len(result_dict[caller_key]),
		"prec": result_dict['precision'],
		'recall': result_dict['recall'],
		"tp": result_dict['tp'],
		"fp": result_dict['fp'],
		'fn': result_dict['fn'],
		"caller": caller_key,
	}


pd.set_option('display.width', 250)


def process_gatk_samples():
	avail_samples = [x for x in os.listdir("data/synthetic/results/gatk/") if x.endswith(".vcf")]
	for _, samp in time_iterator(avail_samples, logger):
		name = samp.split("_on_")[0]
		try:
			res = process_gatk_sample(sample_name=name)
		except Exception as e:
			print "failed on sample", samp
			raise e
		yield res


def process_varscan_samples():
	avail_samples = [x for x in os.listdir("data/synthetic/results/varscan/") if x.endswith(".vcf")]
	for _, samp in time_iterator(avail_samples, logger):
		name = samp.split("_on_")[0]
		try:
			res = process_varscan_sample(sample_name=name)
		except Exception as e:
			print "failed on sample", samp
			raise e
		yield res


def process_micado_samples():
	avail_samples = [x for x in os.listdir("../micado_synthetic_results/synthetic/") if x.endswith(".significant_alterations.json")]
	for _, samp in time_iterator(avail_samples, logger):
		name = samp.split(".")[0]
		try:
			res = process_micado_sample(sample_name=name)
		except Exception as e:
			print "failed on sample", samp
			raise e
		yield res


gatk_aggregated_results = pd.DataFrame.from_records([flatten_sample(x, "gatk") for x in process_gatk_samples()])
gatk_aggregated_results.to_csv("data/synthetic/summary/gatk_results_on_synthetic_data.csv")

varscan_aggregated_results = pd.DataFrame.from_records([flatten_sample(x, "varscan") for x in process_varscan_samples()])
varscan_aggregated_results.to_csv("data/synthetic/summary/varscan_results_on_synthetic_data.csv")

micado_aggregated_results = pd.DataFrame.from_records([flatten_sample(x, "micado") for x in process_micado_samples()])
micado_aggregated_results.to_csv("data/synthetic/summary/micado_results_on_synthetic_data.csv")
