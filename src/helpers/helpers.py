import os
import time
import datetime
# from numpy import mean

__author__ = 'hayssam'

def mean(vec):
	return sum(vec)*1.0/len(vec)

def time_iterator(an_iter, logger, delta_percent=0.01, msg_prefix=None, tot_items=None):
	tot_items = tot_items or len(an_iter)
	percent_increase = max(int(tot_items * delta_percent), 1)
	percent_increase_times = []
	start_time = time.time()
	last_percent = start_time

	for i, doc in enumerate(an_iter):

		yield i, doc

		if (i % percent_increase ) == 0:
			elapsed = time.time() - start_time
			elapsed_since_last = time.time() - last_percent
			last_percent = time.time()
			percent_increase_times.append(elapsed_since_last)
			estimate_completion = mean(percent_increase_times) * (1 / delta_percent) + start_time

			finish_date = datetime.datetime.fromtimestamp(estimate_completion).strftime('%H:%M:%S')
			if msg_prefix:
				logger.info("%s - Iterated %d (%2.f %%) items in %d secs, finish @ %s", msg_prefix, i, i * 100.0 / tot_items, elapsed, finish_date)
			else:
				logger.info("Iterated %d (%2.f %%) items in %d secs, finish @ %s", i, i * 100.0 / tot_items, elapsed, finish_date)
	logger.info("Iterated all %d items", tot_items)


def get_or_create_dir(dirname):
	if not os.path.isdir(dirname):
		os.makedirs(dirname)
	return dirname
