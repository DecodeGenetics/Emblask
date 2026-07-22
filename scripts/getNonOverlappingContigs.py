#!/usr/bin/env python3

import sys
import os
import math
import pysam
import argparse

from pathlib import Path
from operator import itemgetter
from intervaltree import Interval, IntervalTree

def loadVCF(in_vcf, min_gq):

	var_d = {}
	vcf_fn = pysam.VariantFile(in_vcf)
	samples = list(vcf_fn.header.samples)

	if (len(samples) != 1): sys.exit("Input VCF must contain exactly one sample name but contains instead " + str(len(samples)) + "sample name(s).")

	sample = samples[0]

	for rec in vcf_fn.fetch():

		filt = list(rec.filter)

		if (len(filt) == 1) and ((filt[0] == "PASS") or (filt[0] == ".")): # Only keep pass variants

			if rec.samples[sample].phased and (len(rec.ref) == 1) and (len(rec.alts) == 1) and (len(rec.alts[0]) == 1): # Only keep phased SNPs with exactly one alternative allele

				gq = int(rec.samples[sample]["GQ"])

				if (gq >= min_gq):

					if (rec.chrom not in var_d): var_d[rec.chrom] = IntervalTree()

					var_d[rec.chrom][(rec.pos-1):(rec.pos+len(rec.ref)-1)] = gq

	return var_d

if __name__ == "__main__":
	
	# Parse arguments
	parser = argparse.ArgumentParser(prog='getNonOverlappingContigs', description='Given an haplotype asm HA aligned to another hap asm HR, list HA contigs which do not overlap HR or diverge too much from HR.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-b', '--in_bam', action='store', help='Input BAM filename', required=True)
	required.add_argument('-v', '--in_vcf', action='store', help='Input *phased* VCF filename', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-m', '--min_mapq', action='store', help='Minimum mapping quality of alignments to consider', type=int, default=5, required=False)
	optional.add_argument('-g', '--min_gq', action='store', help='Minimum genotype quality of variants to consider', type=int, default=30, required=False)
	optional.add_argument('-d', '--min_div_rate', action='store', help='Minimum divergence rate of an HA contig to report it as non-overlapping', type=float, default=0.1, required=False)
	optional.add_argument('-i', '--min_sv_len', action='store', help='Minimum SV length in an HA contig to report it as non-overlapping', type=int, default=8000, required=False)
	optional.add_argument('-u', '--min_unaligned_len', action='store', help='Minimum unaligned length of an an HA contig to report it as non-overlapping', type=int, default=10000, required=False)
	optional.add_argument('-e', '--min_extend_len', action='store', help='Minimum increased length of an an HA contig to report it as non-overlapping', type=int, default=3000, required=False)
	optional.add_argument('-t', '--threads', action='store', help='Number of threads for compression/decompressing BAM', type=int, default=1, required=False)

	args = parser.parse_args()

	d_hr_len = {} # HR haplotig <name, lengths>
	d_map = {} # HA haplotigs which map to HR and we want to keep
	d_unmap = {} # HA haplotigs which do not map to HR and we want to keep
	d_var = loadVCF(args.in_vcf, args.min_gq) # Load phased SNPs for haplotigs

	in_bamf = pysam.AlignmentFile(args.in_bam, "rb", threads=args.threads) # Open input BAM

	# 0 - Extract HR contig names and lengths
	in_header_d = in_bamf.header.to_dict()

	for contigs in in_header_d["SQ"]: d_hr_len[contigs["SN"]] = int(contigs["LN"])

	# 1 - Gather Alternate Haplotype information
	for rec in in_bamf.fetch(until_eof=True):

		if (rec.query_length != 0) and (rec.is_unmapped == False) and (rec.is_secondary == False):

			if (not rec.has_tag("AS")): sys.exit("BAM file contain alignment(s) without alignment score (tag AS)")
			if (not rec.has_tag("NM")): sys.exit("BAM file contain alignment(s) without the number of mismatches and gaps (tag NM)")

			has_sv = False
			query_length = rec.query_length

			for operation, length in rec.cigartuples:

				if ((operation == 1) or (operation == 2)) and (length >= args.min_sv_len): # Large insertion or deletion

					has_sv = True
					break

			# The query length never includes hard clips so we had it
			if (rec.cigartuples[0][0] == 5): query_length += rec.cigartuples[0][1]
			if (rec.cigartuples[-1][0] == 5): query_length += rec.cigartuples[-1][1]

			rec_info = (rec.reference_name, rec.reference_start, rec.reference_end, rec.query_alignment_start, rec.query_alignment_end, query_length, rec.mapping_quality, has_sv, rec.get_tag("NM"), rec.get_tag("AS"), not rec.is_supplementary, not rec.is_reverse)

			if (rec.query_name not in d_map): d_map[rec.query_name] = []

			d_map[rec.query_name].append(rec_info)

		elif (rec.is_unmapped == True): d_unmap[rec.query_name] = None

	in_bamf.close()

	# 2 - Filter contigs
	hr2ha_d = {}

	for q_n, q_info in d_map.items():

		q_info_sort = sorted(q_info, key=itemgetter(10, 9, 8, 7, 0, 1), reverse=True) # Sort alignments by 1) primary first, supplementary second, 2) alignment score 3) number of errors 4) ref contig name 4) alignment position
		q_info_sort_keep = [] # non-overlapping segments we keep

		r_itv_d = {} # Interval trees (one per chr) of reference intervals we selected
		q_itv = IntervalTree() # Interval tree of query intervals we selected

		for segment in q_info_sort:

			r_n, r_s, r_e, q_s, q_e, q_l, mapq, has_sv, nm_tag, as_tag, is_prim, strand = segment # Segment information
			valid_segment = True # Condition to keep this segment is that it does not overlap another segment on the reference or within the query and it has a good enough mapq

			if not strand:
		
				rev_q_s = q_l - q_e
				rev_q_e = q_l - q_s

				q_s = rev_q_s
				q_e = rev_q_e
			
			if not is_prim: # If this is the primary alignment segment, we keep

				if mapq < args.min_mapq: valid_segment = False
				elif (r_n in r_itv_d) and (len(r_itv_d[r_n][r_s:r_e]) >= 1): valid_segment = False # Overlaps the reference
				elif (len(q_itv[q_s:q_e]) >= 1): valid_segment = False # Overlaps the query

			if valid_segment:

				if (r_n not in r_itv_d): r_itv_d[r_n] = IntervalTree()

				q_info_sort_keep.append(segment)
				q_itv[q_s:q_e] = None
				r_itv_d[r_n][r_s:r_e] = None

		tot_q_align_len = 0
		tot_q_nm = 0
		tot_has_sv = False

		q_length = 0

		r_nb_gq40 = 0
		r_nb_gq30 = 0

		q_nb_gq40 = 0
		q_nb_gq30 = 0

		for segment in q_info_sort_keep:

			r_n, r_s, r_e, q_s, q_e, q_l, mapq, has_sv, nm_tag, as_tag, is_prim, strand = segment # Segment information

			if (q_n == "contig_10128_1-39079"):

				print(str(r_n) + "\t" + str(r_s) + "\t" + str(r_e) + "\t" + str(q_n) + "\t" + str(q_s) + "\t" + str(q_e) + "\t" + str(strand))

			if not strand:
		
				rev_q_s = q_l - q_e
				rev_q_e = q_l - q_s

				q_s = rev_q_s
				q_e = rev_q_e

			tot_q_align_len += q_e - q_s
			tot_q_nm += nm_tag

			q_length = max(q_length, q_l)

			if has_sv: tot_has_sv = True

			if (r_n in d_var):

				for var in d_var[r_n][r_s:r_e]:

					r_nb_gq30 += 1

					if (var.data >= 40): r_nb_gq40 += 1

			if (q_n in d_var):

				for var in d_var[q_n][q_s:q_e]:

					q_nb_gq30 += 1

					if (var.data >= 40): q_nb_gq40 += 1

		#if (r_nb_gq40 >= 1) or (r_nb_gq30 >= 2) or (q_nb_gq40 >= 1) or (q_nb_gq30 >= 2): print(q_n + "\t+") # Haplotig maps to another haplotig which is maybe a contig as it has alignments supporting a heterozygous variant
		#elif (tot_has_sv == True): print(q_n + "\t+") # Contig has an SV with length above threshold
		#elif (tot_q_nm / q_length >= args.min_div_rate): print(q_n + "\t+") # Diverge ratio between the 2 haplotypes is too large
		#elif (q_length - tot_q_align_len >= args.min_unaligned_len): print(q_n + "\t+") # Unaligned part of the contig is too large
		#elif (len(r_itv_d) == 1): # HA haplotig maps to a single HR haplotig
		#
		#	hr_n = list(r_itv_d.keys())[0]
		#
		#	if hr_n not in hr2ha_d: hr2ha_d[hr_n] = []
		#
		#	hr2ha_d[hr_n].append((q_n, q_length))

	#for hr_n, hr_segments in hr2ha_d.items():

	#	if (len(hr_segments) == 1) and (hr_segments[0][1] >= d_hr_len[hr_n] + args.min_extend_len):

	#		print(hr_segments[0][0] + "\t+\n")
	#		print(hr_n + "\t-\n")

	#for q_n in d_unmap.keys(): print(q_n + "\t+") # Contig does not align
