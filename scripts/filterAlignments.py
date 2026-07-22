#!/usr/bin/env python3

import sys
import pysam
import argparse

if __name__ == "__main__":

	# Parse arguments
	parser = argparse.ArgumentParser(prog='filterAlignments', description='Filter alignments', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-i', '--in_bam', action='store', help='Input BAM filename', required=True)
	required.add_argument('-o', '--out_bam', action='store', help='Output BAM filename', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-m', '--min_mapq', action='store', help='Filter out alignments with less than this mapping quality.', type=int, default=0, required=False)
	optional.add_argument('-l', '--min_len', action='store', help='Filter out alignments shorter than this length if MAPQ!=60.', type=int, default=0, required=False)
	optional.add_argument('-e', '--min_error_rate', action='store', help='Filter out alignments with error rate larger than this ratio ("de" tag required in BAM records). 1.0 is no filtering.', type=float, default=1.0, required=False)
	optional.add_argument('-u', '--keep_unmapped', action='store_true', help='Keep unmapped alignments in output BAM.', default=False, required=False)
	optional.add_argument('-t', '--threads', action='store', help='Number of threads for compressing/decompressing BAM', type=int, default=1, required=False)

	args = parser.parse_args()

	# Main program
	in_bamf = pysam.AlignmentFile(args.in_bam, "rb", threads=args.threads)
	out_bamf = pysam.AlignmentFile(args.out_bam, "wb", template=in_bamf, threads=args.threads)

	count_rm_map = 0
	count_rm_umap = 0

	for rec in in_bamf.fetch(until_eof=args.keep_unmapped):

		# Only keep unmapped records if required 
		if rec.is_unmapped:

			if args.keep_unmapped: out_bamf.write(rec)
			else: count_rm_umap += 1

		# Keep alignment if:
		# - Aligned length and MAPQ are greater or equal to minimum required
		# - Aligned length is less than minimum required but MAPQ is 60 (the maximum)
		elif ((rec.query_alignment_length >= args.min_len) and (rec.mapping_quality >= args.min_mapq)) or ((rec.query_alignment_length < args.min_len) and (rec.mapping_quality == 60)):

			if (args.min_error_rate >= 1.0): out_bamf.write(rec) # No filtering on error rate
			else:

				if not rec.has_tag("de"): sys.exit("Record found without the de tag which is required for error rate filtering. Abort.")

				if (float(rec.get_tag("de")) < args.min_error_rate): out_bamf.write(rec)
				else: count_rm_map += 1

		else: count_rm_map += 1

	in_bamf.close()
	out_bamf.close()

	print("Removed " + str(count_rm_map) + " aligned records.")
	print("Removed " + str(count_rm_umap) + " unaligned records.")
