#!/usr/bin/env python3

import sys
import pysam
import argparse

from intervaltree import Interval, IntervalTree

def getOverlap(a, b):

	return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def loadBED(bed_fn):

	f = open(bed_fn, "r")
	regions_l = []

	for line in f:

		tokens = line.split()

		chrom = tokens[0]
		p_start = int(tokens[1]) # 0-based
		p_end = int(tokens[2])

		regions_l.append((chrom, p_start, p_end))

	f.close()

	return regions_l

#def loadBAM(bam_fn, threads_bam, regions_l):
#
#	in_bamf = pysam.AlignmentFile(bam_fn, "rb", threads=threads_bam)
#
#	for chrom, p_start, p_end in regions_l:
#
#		nb_softclip_7kb = 0
#
#		for rec in in_bamf.fetch(chrom, p_start, p_end):
#
#			if ((rec.cigartuples[0][0] == 4) and (int(rec.cigartuples[0][1]) >= 7000)) or ((rec.cigartuples[-1][0] == 4) and (int(rec.cigartuples[-1][1]) >= 7000)): nb_softclip_7kb += 1
#
#			if (nb_softclip_7kb >= 2):
#
#				print(chrom + "\t" + str(p_start) + "\t" + str(p_end))
#
#				break

def loadBAM(bam_fn, regions_l, min_len_align, threads_bam):

	in_bamf = pysam.AlignmentFile(bam_fn, "rb", threads=threads_bam)

	for chrom, p_start, p_end in regions_l:

		sa_d = {}
		valid_ovlp = False

		for rec in in_bamf.fetch(chrom, p_start, p_end): # iterate over alignments of the region

			if not rec.is_secondary and not rec.is_unmapped: # Only iterate over primary and supplementary alignments

				# Alignment has supplementary alignments and at least one soft clipping longer than required alignment length
				if (rec.has_tag("SA")) and (((rec.cigartuples[0][0] == 4) and (int(rec.cigartuples[0][1]) >= min_len_align)) or ((rec.cigartuples[-1][0] == 4) and (int(rec.cigartuples[-1][1]) >= min_len_align))):

					sa_l = rec.get_tag("SA").split(";") # Get list of supplementary alignment from SA tag

					for sa in sa_l: # Each SA tag element is (rname, pos, strand, CIGAR, mapQ, NM)

						if (len(sa) != 0): # Make sure SA elem is not empty

							sa_rn, sa_pos_str, sa_strand, sa_cigar, sa_mapq, sa_nm = sa.split(",") # Split SA tag into (rname, pos, strand, CIGAR, mapQ, NM) tokens

							sa_pos = int(sa_pos_str)

							if (sa_rn != rec.reference_name) or (sa_pos != rec.reference_start): # Make sure SA element is not the currently analyzed alignment

								# Compute from SA element CIGAR the reference length 
								len_align_ref = 0
								curr_len_op = ""

								for c in sa_cigar:

									if c.isdigit(): curr_len_op += c
									else:

										if (c == "M") or (c == "D"): len_align_ref += int(curr_len_op)
										
										curr_len_op = ""

								if (len_align_ref >= min_len_align): # If reference length longer than required alignment length

									if sa_rn not in sa_d: # New contig/haplotig 

										sa_d[sa_rn] = IntervalTree()
										sa_d[sa_rn][sa_pos:sa_pos+len_align_ref] = None

									else: # Go over list of positions for that contig/haplotig that were already covered by other SA elements

										sa_itv = list(sa_d[sa_rn][sa_pos:sa_pos+len_align_ref])

										for itv in sa_itv:

											if (getOverlap((sa_pos, sa_pos+len_align_ref), (itv.begin, itv.end)) >= min_len_align): # If an overlap of the required length is found

												valid_ovlp = True
												break

										if not valid_ovlp: sa_d[sa_rn][sa_pos:sa_pos+len_align_ref] = None

				if valid_ovlp:

					print(chrom + "\t" + str(p_start) + "\t" + str(p_end))
					break

if __name__ == "__main__":

	parser = argparse.ArgumentParser(prog='test', description='test', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-b', '--in_bam', action='store', help='Input BAM filename', required=True)
	required.add_argument('-r', '--in_regions', action='store', help='Input BED filename', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-t', '--threads', action='store', help='Number of threads for compressing/decompressing BAM', type=int, default=1, required=False)
	optional.add_argument('-l', '--min_len_align', action='store', help='Minimum alignment length of supplementary alignments to other contigs/haplotigs', type=int, default=3000, required=False)

	args = parser.parse_args()

	regions_l = loadBED(args.in_regions)
	loadBAM(args.in_bam, regions_l, args.min_len_align, args.threads)
