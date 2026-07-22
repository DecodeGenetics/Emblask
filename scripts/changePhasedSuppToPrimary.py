#!/usr/bin/env python3

import sys
import os
import math
import pysam
import argparse

from pathlib import Path
from operator import itemgetter
from intervaltree import Interval, IntervalTree

if __name__ == "__main__":
	
	# Parse arguments
	parser = argparse.ArgumentParser(prog='changePhasedSuppToPrimary', description='Change phased supplementary alignments into primary and set matching unphased primary to supplementary.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-b', '--in_bam', action='store', help='Input BAM filename', required=True)
	required.add_argument('-B', '--out_bam', action='store', help='Output BAM filename', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-t', '--threads', action='store', help='Number of threads for compression/decompressing BAM', type=int, default=1, required=False)

	args = parser.parse_args()

	# Start
	in_bamf = pysam.AlignmentFile(args.in_bam, "rb", threads=args.threads) # Open input BAM

	read_d = {}

	for rec in in_bamf.fetch():

		if (rec.query_length != 0) and (rec.is_unmapped == False) and (rec.is_secondary == False) and rec.has_tag("HP"):

			if (not rec.has_tag("AS")): sys.exit("BAM file contain primary/supplementary alignment(s) without Alignment Score (tag AS)")

			if (rec.query_name not in read_d): read_d[rec.query_name] = (False, -1) # (Whether primary is phased, best alignment score of phased supplementary where -1 means no phased supplementary)

			nb_phase_align = list(read_d[rec.query_name])

			if rec.is_supplementary: nb_phase_align[1] = max(nb_phase_align[1], int(rec.get_tag("AS"))) # Phased supplementary: save which phased supplementary has best alignment score
			else: nb_phase_align[0] = True # Phased primary: save that read has phased primary

			read_d[rec.query_name] = tuple(nb_phase_align)

	in_bamf.close()

	in_bamf = pysam.AlignmentFile(args.in_bam, "rb", threads=args.threads) # Open input BAM
	out_bamf = pysam.AlignmentFile(args.out_bam, "wb", template=in_bamf, threads=args.threads) # Open output BAM

	for rec in in_bamf.fetch(until_eof=True):

		if (rec.query_length != 0) and (rec.is_unmapped == False) and (rec.is_secondary == False) and (rec.query_name in read_d):

			nb_phase_align = read_d[rec.query_name]

			if (nb_phase_align[0] == False) and (nb_phase_align[1] != -1): # Read has no phased primary but at least one phased supplementary

				if rec.is_supplementary and (int(rec.get_tag("AS")) == nb_phase_align[1]) and rec.has_tag("HP"): # Record is phased supplementary with best alignment score seen for that read

					rec.is_supplementary = False # Supplementary becomes primary
					nb_phase_align = read_d[rec.query_name]
					read_d[rec.query_name] = (nb_phase_align[0], 0) # Only one supplementary can become primary so set AS to 0
					
				elif (rec.is_supplementary == False): rec.is_supplementary = True # Record is unphased primary: becomes supplementary

		out_bamf.write(rec)

	in_bamf.close()
	out_bamf.close()
		
	
	
			

