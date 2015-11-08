from collections import defaultdict
import os

import pandas as pd
from helpers.helpers import time_iterator
from helpers.logger import init_logger

pd.set_option('display.width', 250)
import random
import itertools
import simplejson
import binascii

__author__ = 'hayssam'

logger = init_logger("POSPROCESS")


def hash_dict(d):
	return hash(tuple(sorted(d.items())))


# @@@

# find closest alterations

def decorate_alt(alt, origin):
	d = {"origin": origin, "start": alt['start'], "alt_type": alt['alt_type'], "len": alt['alt_length'], "sequence": alt["alt_content"]}
	d['hash'] = hash_dict(d)
	return d


def decorate_micado_alt(alt):
	d = {
		"sequence": alt['alt_sequence'],
		"origin": "micado",
		"start": alt['start'],
		"pvalue_ratio": alt['pvalue_ratio'],
		"alt_type": alt['alt_type'],
		"len": alt['end'] - alt['start'] + 1 if alt['alt_type'] == "I" else 0,
		"cigar": alt['compact_cigar'],
		"ref_seq": alt['alignment'][0],
		"alt_seq": alt['alignment'][1],
		"is_multi": alt['is_multi']
	}
	d['hash'] = hash_dict(d)
	return d


def closest_alteration(alt, alt_list):
	closest = min(alt_list, key=lambda x: abs(x['start'] - alt['start']))

	return closest


def is_match(x):
	return abs(x['src']['start'] - x['closest']['start']) <= 20 and x['src']['alt_type'] == x['closest']['alt_type']


def remove_fp(matched_alterations):
	return list(itertools.ifilter(is_match, matched_alterations))


def label(matched_alterations):
	for m in matched_alterations:
		m['is_match'] = is_match(m)
	return matched_alterations


def build_xp_description(results):
	xp_metadata_keys = ['timestamp', 'git_revision_hash', 'n_reads', 'exec_time']
	sampler_metadata_keys = ['alt_weight', 'fraction_altered', 'max_len', 'n_alterations', 'n_reads', 'seed']
	micado_metadata_keys = ['kmer_length', 'min_support_percentage', 'npermutations']
	xp_metadata = {k: results[k] for k in xp_metadata_keys}
	xp_metadata.update({"tool_sampler_" + k: results['sampler']['parameters'][k] for k in sampler_metadata_keys})
	xp_metadata.update({"tool_micado_" + k: results['parameters'][k] for k in micado_metadata_keys})
	xp_metadata['uuid'] = binascii.hexlify(str(hash_dict(xp_metadata)))
	return xp_metadata


def build_alteration_pair_description(pair):
	if pair['src']['origin'] == "micado":
		micado_res = pair['src']
		if pair['is_match']:
			injected_res = pair['closest']
		else:
			injected_res = defaultdict(lambda: None)
		micado_additional_data = {"micado_" + k: v for k, v in pair['src'].items()}
	else:
		injected_res = pair['src']
		micado_additional_data = {}
		if pair['is_match']:
			micado_res = pair['closest']
		else:
			micado_res = defaultdict(lambda: None)
	record = {
		"micado_hash": micado_res['hash'],
		"injected_hash": injected_res['hash'],
		"injected_pos": injected_res['start'],
		"micado_pos": micado_res['start'],
		"injected_alt_type": injected_res['alt_type'],
		"micado_alt_type": micado_res['alt_type'],
		"micado_len": micado_res['len'],
		"injected_len": injected_res['len'],
		"micado_sequence": micado_res['sequence'],
		"injected_sequence": injected_res['sequence'],
		"is_match": pair['is_match'],
	}
	record.update(micado_additional_data)
	return record


def tabulate_result(results):
	result_table = []
	xp_metadata = build_xp_description(results)

	identified_alterations_dict = map(lambda x: decorate_micado_alt(x), results.get('significant_alterations', []))
	injected_alteration_dict = map(lambda x: decorate_alt(x, "sampler"), results['sampler']['injected_alterations'])

	identified_to_injected = label(map(lambda x: {"src": x, "closest": closest_alteration(x, injected_alteration_dict)}, identified_alterations_dict))

	if len(identified_alterations_dict) < 1:
		injected_to_identified = [{"src": x, "closest": None, "is_match": False} for x in injected_alteration_dict]
	else:
		injected_to_identified = label(map(lambda x: {"src": x, "closest": closest_alteration(x, identified_alterations_dict)}, injected_alteration_dict))

	accounted_alt = set()

	for pair in identified_to_injected + injected_to_identified:
		if pair['src']['hash'] in accounted_alt:
			continue

		record = build_alteration_pair_description(pair)
		record.update(xp_metadata)
		accounted_alt.add(record['micado_hash'])
		accounted_alt.add(record['injected_hash'])
		result_table.append(record)
	return result_table


results_dir = "../micado_synthetic_results/synthetic/"
avail_results = [results_dir + x for x in os.listdir(results_dir) if x.endswith(".json") and "combined" in x]
len(avail_results)

result_table = []
for i, input_json in time_iterator(avail_results, logger=logger):
	# input_json = random.choice(avail_results)
	with open(input_json, "r") as f:
		try:
			result_dict = simplejson.load(f)
		except simplejson.JSONDecodeError:
			logger.critical("Malformed json file %s",input_json)
			continue

	# result_dict.keys()
	# result_dict['sampler']['injected_alterations']
	# result_dict['significant_alterations']
	this_result_table = tabulate_result(result_dict)
	result_table.extend(this_result_table)

all_results = pd.DataFrame.from_records(result_table)
all_results.to_csv("data/summary/results_on_synthetic_data.csv")
