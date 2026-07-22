#!/usr/bin/env python3

import sys
import os
import random
import argparse
import pysam

import multiprocessing as mp

from itertools import groupby
from operator import itemgetter

def replaceIUPAC(seq=None):

	dna_alpha=['A', 'C', 'G', 'T']
	
	return ''.join([c if c in dna_alpha else random.choice(dna_alpha) for c in seq.upper()])

#def getLowBQIntervals(qual, min_bq, max_len_lowbq):
#
#	qual_int = [ord(x)-33 for x in qual]
#	qual_bool = [(x >= min_bq) for x in qual_int]
#	
#	runs_lowBQ = [[i for i, _ in group] for key, group in groupby(enumerate(qual_bool), key=itemgetter(1)) if not key]
#	itv_runs_lowBQ = [(x[0], x[-1]) for x in runs_lowBQ if (len(x)>=max_len_lowbq)]
#	
#	return itv_runs_lowBQ

def getLowBQIntervals(qual, min_bq, max_len_lowbq, max_ratio_low_bq):

	max_count_low_bq = int(max_len_lowbq * max_ratio_low_bq) # Max number of low base quality bases within a sliding window

	qual_int = [ord(x)-33 for x in qual] # Qualities in int format
	qual_bool = [(x >= min_bq) for x in qual_int] # Quality for each position is either >= min_bq (True) or < min_bq (False)

	del qual_int

	qual_count_windows = [0] * (len(qual_bool) - max_len_lowbq + 1)
	qual_count_windows[0] = qual_bool[:max_len_lowbq].count(False)

	for i in range(1, len(qual_bool) - max_len_lowbq):

		qual_count_windows[i] = qual_count_windows[i-1] + int(qual_bool[i+max_len_lowbq-1] == False) - int(qual_bool[i-1] == False)

	del qual_bool

	qual_bool_windows = [(qual_count_windows[i] < max_count_low_bq) for i in range(len(qual_count_windows))]

	del qual_count_windows

	runs_lowBQ = [[i for i, _ in group] for key, group in groupby(enumerate(qual_bool_windows), key=itemgetter(1)) if not key]
	itv_runs_lowBQ = [(x[0], x[-1] + max_len_lowbq) for x in runs_lowBQ]

	if (len(itv_runs_lowBQ) <= 1): return itv_runs_lowBQ
	else:

		# Merge overlapping intervals
		itv_runs_lowBQ.sort(key=lambda x: x[0])
		merged_itv_runs_lowBQ = [itv_runs_lowBQ[0]]

		for curr in itv_runs_lowBQ:

			prev = merged_itv_runs_lowBQ[-1]

			if (curr[0] <= prev[1]): merged_itv_runs_lowBQ[-1] = (prev[0], max(prev[1], curr[1]))
			else: merged_itv_runs_lowBQ.append(curr)

		return merged_itv_runs_lowBQ

def splitRead(title, seq, qual, min_len_seq, min_bq, max_len_lowbq, max_ratio_low_bq):

	out = []
	len_rec = len(seq)

	if (len_rec >= min_len_seq):

		lowBQIntervals = getLowBQIntervals(qual, min_bq, max_len_lowbq, max_ratio_low_bq)
		id_seq = 0

		if (len(lowBQIntervals) == 0): out.append((title + "/" + str(id_seq), replaceIUPAC(seq), qual))
		else:

			prev_pos = 0

			for itv in lowBQIntervals:

				if (itv[0]-prev_pos >= min_len_seq):

					out.append((title + "/" + str(id_seq), replaceIUPAC(seq[prev_pos:itv[0]]), qual[prev_pos:itv[0]]))

					id_seq += 1

				prev_pos = itv[1] + 1

			if (len_rec-prev_pos >= min_len_seq): out.append((title + "/" + str(id_seq), replaceIUPAC(seq[prev_pos:]), qual[prev_pos:]))

	return out

def splitReads(reads, min_len_seq, min_bq, max_len_lowbq, max_ratio_low_bq):

	return [splitRead(x[0], x[1], x[2], min_len_seq, min_bq, max_len_lowbq, max_ratio_low_bq) for x in reads]

def printRecord(title, seq, qual):

	print("@" + title) # Name
	print(seq) # Sequence
	print("+")
	print(qual) # Qualities

def printRecords(l_records):

	for rec in l_records:

		print("@" + rec[0]) # Name
		print(rec[1]) # Sequence
		print("+")
		print(rec[2]) # Qualities

if __name__ == "__main__":

	# Parse arguments
	parser = argparse.ArgumentParser(prog='splitAtLowBQ', description='Split corrected long reads at low base quality positions', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-f', '--in_long_fastq', action='append', help='Filename of a corrected long read FASTQ file', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-t', '--threads', action='store', help='Number of threads to use', type=int, default=1, required=False)
	optional.add_argument('-c', '--chunk_per_thread', action='store', help='Number of FASTQ records per thread', type=int, default=1000, required=False)
	optional.add_argument('-l', '--min_length', action='store', help='Minimum length of output reads', type=int, default=5000, required=False)
	optional.add_argument('-b', '--min_bq', action='store', help='Minimum base quality', type=int, default=5, required=False)
	optional.add_argument('-w', '--min_window_size', action='store', help='Minimum length of sliding windows', type=int, default=100, required=False)
	optional.add_argument('-r', '--max_ratio_bq', action='store', help='Maximum ratio of low base quality within sliding windows for split', type=float, default=0.9, required=False)

	args = parser.parse_args()

	if (args.min_length < 1): raise SystemError("Minimum length of output reads cannot be less than 1.")
	if (args.min_bq < 0): raise SystemError("Cannot use base quality less than 0.")
	if (args.min_bq > 90): raise SystemError("Cannot use base quality greater than 90.")
	if (args.min_window_size < 1): raise SystemError("Minimum length of sliding windows cannot be less than 1.")
	if (args.max_ratio_bq < 0.0): raise SystemError("Maximum ratio of low base quality within sliding windows cannot be less than 0.0.")
	if (args.max_ratio_bq > 1.0): raise SystemError("Maximum ratio of low base quality within sliding windows cannot be greater than 1.0.")

	for fn in args.in_long_fastq:
	
		if (len(fn) == 0): raise SystemError("Error: Specify file name\n")
		if not os.path.exists(fn): raise SystemError("Error: File does not exist\n")

	for fn in args.in_long_fastq:

		with pysam.FastxFile(fn) as fh:

			if (args.threads == 1):

				for rec in fh:

					if (len(rec.sequence) >= args.min_length):

						rec_out = splitRead(rec.name, rec.sequence, rec.quality, args.min_length, args.min_bq, args.min_window_size, args.max_ratio_bq)

						printRecords(rec_out)

			else:

				records_in = [ [] for _ in range(args.threads) ]
				records_out = []

				id_rec = 0

				for rec in fh:

					if (len(rec.sequence) >= args.min_length):

						records_in[id_rec % args.threads].append((rec.name, rec.sequence, rec.quality))

						id_rec += 1

						if (id_rec == args.threads * args.chunk_per_thread):

							thread_pool = mp.Pool(processes = args.threads)

							for i in range(args.threads):

								thread_pool.apply_async(splitReads, args=(records_in[i], args.min_length, args.min_bq, args.min_window_size, args.max_ratio_bq), callback=records_out.append)

							thread_pool.close() 
							thread_pool.join()

							for i in range(args.threads):

								for rec_out in records_out[i]: printRecords(rec_out)

							records_in = [ [] for _ in range(args.threads) ]
							records_out = []

							id_rec = 0

				if (id_rec != 0):

					thread_pool = mp.Pool(processes = args.threads)

					for i in range(args.threads):

						thread_pool.apply_async(splitReads, args=(records_in[i], args.min_length, args.min_bq, args.min_window_size, args.max_ratio_bq), callback=records_out.append)

					thread_pool.close() 
					thread_pool.join()

					for i in range(args.threads):

						for rec_out in records_out[i]: printRecords(rec_out)
