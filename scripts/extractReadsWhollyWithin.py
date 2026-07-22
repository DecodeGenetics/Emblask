#!/usr/bin/env python3

import sys
import os
import math
import pysam
import argparse

from pathlib import Path
from intervaltree import Interval, IntervalTree

def loadBEDregions(bed_fn):

	regions_d = {}
	f_bed_fn = open(bed_fn, "r")

	for record in f_bed_fn:

		rec_tokens = record.split()

		contig = rec_tokens[0]
		pos_s = int(rec_tokens[1])
		pos_e = int(rec_tokens[2])

		if (contig not in regions_d): regions_d[contig] = IntervalTree()

		regions_d[contig][pos_s:pos_e] = None

	f_bed_fn.close()

	return regions_d

def extractReads(bam_fn, regions_d, mapq, min_len_align, include_supplementary, include_secondary, bam_output, bam_threads):

	seq = {} # Reads to output
	bam_out = (bam_output != "")

	in_bamf = pysam.AlignmentFile(bam_fn, "rb", threads=bam_threads)

	if bam_out: out_bamf = pysam.AlignmentFile(bam_output, "wb", template=in_bamf, threads=bam_threads)

	for rec in in_bamf.fetch():

		# Record is not null, unmapped, below minimum MAPQ, below minimum aligned length
		if (rec.query_length != 0) and (not rec.is_unmapped) and (rec.mapping_quality >= mapq) and (rec.reference_length >= min_len_align) and (rec.query_name not in seq):

			# Record is either primary or supplementary/secondary alignment if required
			if (not rec.is_supplementary and not rec.is_secondary) or (include_supplementary and rec.is_supplementary) or (include_secondary and rec.is_secondary):

				# Contig is in the provided regions and there exist a BED region fully covering that read alignment
				#if (rec.reference_name in regions_d) and (len(list(regions_d[rec.reference_name].envelop(rec.reference_start, rec.reference_end))) != 0):
				if (rec.reference_name in regions_d):

						itvs = regions_d[rec.reference_name]
						is_wholly_within = False

						for itv in itvs:

							if (rec.reference_start >= itv.begin) and (rec.reference_end <= itv.end):

								is_wholly_within = True
								break

						if is_wholly_within:

							if bam_out: out_bamf.write(rec) # Write to output BAM
							else:

								# Add record name or mate to dictionnary 
								seq[rec.query_name] = None

								if rec.query_qualities: # Record has base qualities -> FASTQ output

									print("@" + rec.query_name + "\n" + rec.query_sequence + "\n+\n" + ''.join(map(lambda x: chr(x+33), rec.query_qualities)))

								else: # Record has no base qualities -> FASTQ output

									print(">" + rec.query_name + "\n" + rec.query_sequence)



	in_bamf.close()

	if bam_out: out_bamf.close()

if __name__ == "__main__":

	# Defaults values
	parser = argparse.ArgumentParser(prog='extractReadsWhollyWithin', description='Extract reads mapping wholly within one or multiple regions', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-b', '--in_bam', action='store', help='Input BAM filename', required=True)
	required.add_argument('-r', '--in_bed', action='store', help='Input BED filename', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-m', '--mapq', action='store', help='Minimum mapping quality of supplementary alignments to consider', type=int, default=0, required=False)
	optional.add_argument('-l', '--min_len_align', action='store', help='Minimum length of an aligned reads to be extracted', type=int, default=3000, required=False)
	optional.add_argument('-s', '--include_supplementary', action='store_true', help='Extract reads from supplementary alignments', default=False, required=False)
	optional.add_argument('-S', '--include_secondary', action='store_true', help='Extract reads from secondary alignments', default=False, required=False)
	optional.add_argument('-B', '--bam_output', action='store', help='Write to BAM file instead of FASTA/FASTQ on stdout', default="", required=False)
	optional.add_argument('-t', '--threads', action='store', help='Number of threads for compressing/decompressing BAM', type=int, default=1, required=False)

	args = parser.parse_args()

	regions_d = loadBEDregions(args.in_bed)

	extractReads(args.in_bam, regions_d, args.mapq, args.min_len_align, args.include_supplementary, args.include_secondary, args.bam_output, args.threads)
