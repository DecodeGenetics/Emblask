#!/usr/bin/env python3

import sys
import os
import math
import pysam
import argparse

from pathlib import Path
from intervaltree import Interval, IntervalTree

def loadReadList(read_list = None):

	readl = {}

	if (read_list != None):

		f_read_list = open(read_list, "r")

		for read_name in f_read_list: readl[read_name.split()[0]] = None

		f_read_list.close()

	return readl

def interval_extract(pos_l):

	out_l = []

	if (len(pos_l) != 0):
	
		prev_pos = pos_l[0]
		range_start = pos_l[0]
	  
		for pos in pos_l[1:]:

			if (pos == prev_pos + 1): prev_pos = pos
			else:
				out_l.append((range_start, prev_pos))
		    		
				prev_pos = pos
				range_start = pos

		out_l.append((range_start, prev_pos))

	return out_l

def loadReadsWithSuppAlign(bamf, read_names, len_supp_align, mapq, bam_threads, sam_coords = None):

	itv_d = {}
	all_read_names = (len(read_names) == 0)

	if (sam_coords is None):

		in_bamf = pysam.AlignmentFile(bamf, "rb", threads=bam_threads) # BAM file containing reads

		for rec in in_bamf.fetch():

			if (rec.query_length != 0) and rec.is_supplementary and (rec.mapping_quality >= mapq) and (rec.query_alignment_length >= len_supp_align) and (all_read_names or (rec.query_name in read_names)):

				if (rec.reference_name not in itv_d): itv_d[rec.reference_name] = IntervalTree()

				itv_d[rec.reference_name][rec.reference_start:rec.reference_end] = None

		in_bamf.close()

	else:

		for sam_coord in sam_coords:

			in_bamf = pysam.AlignmentFile(bamf, "rb", threads=bam_threads) # BAM file containing reads

			for rec in in_bamf.fetch(region=sam_coord):

				if (rec.query_length != 0) and rec.is_supplementary and (rec.mapping_quality >= mapq) and (rec.query_alignment_length >= len_supp_align) and (all_read_names or (rec.query_name in read_names)):

					if (rec.reference_name not in itv_d): itv_d[rec.reference_name] = IntervalTree()

					itv_d[rec.reference_name][rec.reference_start:rec.reference_end] = None

			in_bamf.close()

	for contig_name, contig_itv in itv_d.items():

		contig_itv_m = contig_itv.copy()

		contig_itv_m.merge_overlaps()

		for itv_m in contig_itv_m:

			counts = [0] * (itv_m.end - itv_m.begin)

			for itv in contig_itv[itv_m.begin:itv_m.end]:

				for i in range(itv.begin-itv_m.begin, itv.end-itv_m.begin): counts[i] += 1

			pos = [i for i,x in enumerate(counts) if x >= 2]

			del counts

			pos_itv = interval_extract(pos) # Returns inclusive intervals

			del pos

			for x in pos_itv: print(contig_name + "\t" + str(itv_m.begin + x[0]) + "\t" + str(itv_m.begin + x[1] + 1))
					
if __name__ == "__main__":

	parser = argparse.ArgumentParser(prog='getSuppAlignCoord', description='Extract the coordinates to which supplementary alignments ', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-b', '--in_bam', action='store', help='Input BAM filename', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-r', '--in_read_list', action='store', help='Input list of reads without @ or > prefix', required=False)
	optional.add_argument('-o', '--len_supp_align', action='store', help='Minimum length of supplementary alignments to consider', type=int, default=2000, required=False)
	optional.add_argument('-m', '--mapq', action='store', help='Minimum mapping quality of supplementary alignments to consider', type=int, default=30, required=False)
	optional.add_argument('-s', '--sam_coord', action='append', help='SAM coordinates in which supplementary alignments are considered', required=False)
	optional.add_argument('-t', '--threads', action='store', help='Number of threads for compressing/decompressing BAM', type=int, default=1, required=False)

	args = parser.parse_args()

	read_names = loadReadList(args.in_read_list)

	loadReadsWithSuppAlign(args.in_bam, read_names, args.len_supp_align, args.mapq, args.threads, args.sam_coord)
