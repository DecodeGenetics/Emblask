#!/usr/bin/env python3

import sys
import pysam
import argparse

def computePhaseCoordinates(bam_fns, mapq, bam_threads, bed_fn=None):

	contigs = {}
	contig_lens = {}

	for bam_fn in bam_fns:

		in_bamf = pysam.AlignmentFile(bam_fn, "rb", threads=bam_threads)

		if (bed_fn == None) or (bed_fn == ""):

			for contig_name, contig_len in zip(in_bamf.references, in_bamf.lengths):

				if (contig_name not in contigs): contigs[contig_name] = [0, 0, 0, contig_len]
				elif (contigs[contig_name][3] != contig_len): sys.exit("Input BAM files use different references")

			for rec in in_bamf.fetch():

				if (rec.query_length != 0) and (rec.is_unmapped == False) and (rec.is_secondary == False) and (rec.mapping_quality >= mapq):

					if (rec.has_tag("HP")): contigs[rec.reference_name][rec.get_tag("HP")] += rec.reference_length
					else: contigs[rec.reference_name][0] += rec.reference_length

		else:

			for contig_name, contig_len in zip(in_bamf.references, in_bamf.lengths):

				if (contig_name not in contig_lens): contig_lens[contig_name] = contig_len
				elif (contig_lens[contig_name] != contig_len): sys.exit("Input BAM files use different references")

			bed = open(bed_fn, "r")

			for line in bed:

				tokens = line.split()

				if (tokens[0] in contig_lens) and (int(tokens[1]) < int(tokens[2])):

					cov = [0, 0, 0]
					coord = (tokens[0], max(int(tokens[1]), 0), min(int(tokens[2]), contig_lens[tokens[0]]))
					coord_len = coord[2] - coord[1]

					if (coord not in contigs): contigs[coord] = [0, 0, 0]

					for rec in in_bamf.fetch(coord[0], coord[1], coord[2]):

						if (rec.query_length != 0) and (rec.is_unmapped == False) and (rec.is_secondary == False) and (rec.mapping_quality >= mapq):

							ref_len = min(rec.reference_end, coord[2]) - max(rec.reference_start, coord[1])

							if (rec.has_tag("HP")): contigs[coord][rec.get_tag("HP")] += ref_len
							else: contigs[coord][0] += ref_len

			bed.close()

		in_bamf.close()

	if (bed_fn == None):

		for contig_name, contig_cov in contigs.items():

			no_cov, hp1_cov, hp2_cov, contig_len = contig_cov

			print(contig_name + "\t0\t" + str(contig_len) + "\t" + str(no_cov/contig_len) + "\t" + str(hp1_cov/contig_len) + "\t" + str(hp2_cov/contig_len))

	else:

		for contig_coord, contig_cov in contigs.items():

			coord_len = contig_coord[2] - contig_coord[1]

			print(contig_coord[0] + "\t" + str(contig_coord[1]) + "\t" + str(contig_coord[2]) + "\t" + str(contig_cov[0]/coord_len) + "\t" + str(contig_cov[1]/coord_len) + "\t" + str(contig_cov[2]/coord_len))

if __name__ == "__main__":

	# Parse arguments
	parser = argparse.ArgumentParser(prog='getPhaseCoveragePerContig', description='Compute the phase coverage per contig', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-b', '--in_bam', action='append', help='BAM filename', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-m', '--mapq', action='store', help='Only consider alignments with this minimum quality.', type=int, default=0, required=False)
	optional.add_argument('-r', '--in_regions', action='store', help='Compute phase coverage from all reads overlapping the given BED coordinates', default=None, required=False)
	optional.add_argument('-t', '--threads', action='store', help='Number of threads for compression/decompressing BAM', type=int, default=1, required=False)

	args = parser.parse_args()

	computePhaseCoordinates(args.in_bam, args.mapq, args.threads, args.in_regions)

