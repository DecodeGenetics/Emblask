#!/usr/bin/env python3

import sys
import os
import math
import pysam
import argparse

from intervaltree import Interval, IntervalTree

def getOverlap(a, b):

	return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def loadTSVbinning(in_tsv_fn):

	bin_d = {}
	f_tsv_fn = open(in_tsv_fn, "r")

	for record in f_tsv_fn:

		rec_tokens = record.split()

		contig = rec_tokens[0]
		hp = int(rec_tokens[1])
		phase_conf = int(rec_tokens[2])

		bin_d[contig] = (hp, phase_conf)

	f_tsv_fn.close()

	return bin_d

def detectHapCollision(in_bam_fn, bin_d, threads_bam):

	# 1 
	in_bamf = pysam.AlignmentFile(in_bam_fn, "rb", threads=threads_bam)

	genome_d = {}

	for rec in in_bamf.fetch():

		if (rec.query_length != 0) and (not rec.is_unmapped) and (not rec.is_supplementary) and (not rec.is_secondary): # Record is primary

			if (rec.query_name in bin_d) and (rec.query_alignment_length >= 0.99 * rec.query_length): # Record was trio-binned and aligns on >= 99% of its length

				if (rec.reference_name not in genome_d): genome_d[rec.reference_name] = IntervalTree()
				
				genome_d[rec.reference_name][rec.reference_start:rec.reference_end] = rec.query_name # Save record coordinates

	in_bamf.close()

	# 2
	in_bamf = pysam.AlignmentFile(in_bam_fn, "rb", threads=threads_bam)

	for rec in in_bamf.fetch():

		if (rec.query_length != 0) and (not rec.is_unmapped) and (not rec.is_supplementary) and (not rec.is_secondary): # Record is primary

			if (rec.query_name in bin_d) and (rec.query_alignment_length >= 0.99 * rec.query_length): # Record was trio-binned and aligns on >= 99% of its length

				query_hap, query_phase_conf = bin_d[rec.query_name]
				overlap_contigs_l = list(genome_d[rec.reference_name].overlap(rec.reference_start, rec.reference_end)) # Find contigs overlapping same reference coordinates

				valid_ovlp = []

				for contig_hap in overlap_contigs_l:

					if (rec.query_name != contig_hap.data):

						ovlp = getOverlap((rec.reference_start, rec.reference_end), (contig_hap.begin, contig_hap.end))

						if (ovlp >= (rec.reference_end - rec.reference_start) * 0.99): valid_ovlp.append(contig_hap) # If contig overlaps query on 99% of it, keep

				if (len(valid_ovlp) == 1): # Of query is overlapped by one single contig (multiple would mean something else is going on at that locus)

					for contig_hap in valid_ovlp:

						ovlp_hap, ovlp_phase_conf = bin_d[contig_hap.data]
						
						if (query_hap != 0) and (query_hap == ovlp_hap) and (query_phase_conf > ovlp_phase_conf): # Query and overlap contig share same haplotype -> haplotype collision

							correct_hap = 1 if (query_hap == 2) else 2

							print(rec.query_name + "\t" + str(query_hap) + "\t" + str(correct_hap))
							break

						if (query_hap == 0) and (query_hap != ovlp_hap): # Query is not binned, overlap is binned -> query can be binned

							correct_hap = 1 if (ovlp_hap == 2) else 2

							print(rec.query_name + "\t" + str(query_hap) + "\t" + str(correct_hap))
							break

	in_bamf.close()

if __name__ == "__main__":

	# Defaults values
	parser = argparse.ArgumentParser(prog='detectHaplotypeCollision', description='Detect haplotigs which bin to the same haplotype but should not', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-b', '--in_bam', action='store', help='Input BAM filename', required=True)
	required.add_argument('-p', '--in_phase_tsv', action='store', help='Input phasing TSV filename', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-t', '--threads', action='store', help='Number of threads for compressing/decompressing BAM', type=int, default=1, required=False)

	args = parser.parse_args()

	bin_tsv_d = loadTSVbinning(args.in_phase_tsv)

	detectHapCollision(args.in_bam, bin_tsv_d, args.threads)
