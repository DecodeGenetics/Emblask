#!/usr/bin/env python3

import sys
import pysam
import argparse

from intervaltree import Interval, IntervalTree

def computeHapIDperContig(vcf_fn):

	ps_d = {}
	pys_vcf = pysam.VariantFile(vcf_fn)

	for rec in pys_vcf:

		rec_filter = rec.filter.keys()

		if (len(rec_filter) == 1) and ((rec_filter[0] == "PASS") or (rec_filter[0] == ".")): # Make sure variant is PASS

			rec_samples = rec.samples.keys()

			if (len(rec_samples) == 1): # Does not process multi-sample VCF

				sample = rec_samples[0]

				if ("PS" in rec.samples[sample]): # Variant is phased

					ps = rec.samples[sample]["PS"]
					gt = rec.samples[sample]["GT"]

					chrom = rec.chrom
					pos = int(rec.pos) - 1 # 0-based

					ps_id = (chrom, ps)

					if (ps_id not in ps_d): ps_d[ps_id] = [[pos, pos+1], [0, 0]]

					coord, counts = ps_d[ps_id]

					coord[0] = min(coord[0], pos)
					coord[1] = max(coord[1], pos+1)

					if (gt == (0, 1)): counts = (counts[0]+1, counts[1])
					elif (gt == (1, 0)): counts = (counts[0], counts[1]+1)

					ps_d[ps_id] = [coord, counts]


	pys_vcf.close()

	for ps_id, ps_coord_count in ps_d.items():

		contig = ps_id[0]
		coord, counts = ps_coord_count

		print(contig + "\t" + str(coord[0]) + "\t" + str(coord[1]) + "\t" + str(ps_id[1]) + "\t" + str(counts[0]) + "\t" + str(counts[1]))

	return ps_d

if __name__ == "__main__":

	# Parse arguments
	parser = argparse.ArgumentParser(prog='getPhasesetInfo', description='Get informartion about each phaseset of a VCF', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-v', '--in_vcf', action='store', help='Phased VCF filename (must contain Phase Set IDs in a PS tag)', required=True)

	args = parser.parse_args()

	ps = computeHapIDperContig(args.in_vcf)
