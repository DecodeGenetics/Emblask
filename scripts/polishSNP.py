#!/usr/bin/env python3

import sys
import pysam
import argparse

if __name__ == "__main__":

	# Parse arguments
	parser = argparse.ArgumentParser(prog='polishSNP', description='Polish the SNPs of an assembly', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-a', '--in_asm', action='store', help='Assembly filename in FASTQ format', required=True)
	required.add_argument('-v', '--in_vcf', action='store', help='Variant filename in VCF format', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-g', '--min_gq_snp', action='store', help='Minimum Genotype Quality (GQ) of SNPs to use for polishing.', type=int, default=20, required=False)
	optional.add_argument('-G', '--min_gq_indel', action='store', help='Minimum Genotype Quality (GQ) of indels to use for polishing.', type=int, default=20, required=False)
	optional.add_argument('-d', '--min_dp', action='store', help='Minimum depth of variant position to polish.', type=int, default=5, required=False)
	optional.add_argument('-D', '--max_dp', action='store', help='Maximum depth of variant position to polish. -1 is no limit.', type=int, default=-1, required=False)
	optional.add_argument('-f', '--min_vaf', action='store', help='Minimum VAF of variant to consider.', type=float, default=0.35, required=False)

	args = parser.parse_args()

	# Get started
	asm_fn = pysam.FastxFile(args.in_asm)
	vcf_fn = pysam.VariantFile(args.in_vcf)

	for rec_asm in asm_fn:
	
		seq = list(rec_asm.sequence)

		delta_pos = 0
	
		for rec_vcf in vcf_fn.fetch(rec_asm.name):
	
			if (len(rec_vcf.alts) == 1): # Only use bi-allelic variants
	
				rec_filter = rec_vcf.filter.keys()
	
				if (len(rec_filter) == 1) and ((rec_filter[0] == "PASS") or (rec_filter[0] == ".")): # Only use "PASS" or "." variants
	
					rec_samples = rec_vcf.samples.keys()
	
					if (len(rec_samples) == 1): # Discard record if multi-sample
	
						sample = rec_samples[0]
	
						gq = rec_vcf.samples[sample]["GQ"]
						dp = rec_vcf.samples[sample]["DP"]
						vaf = rec_vcf.samples[sample]["VAF"][0]
	
						if (dp >= args.min_dp) and ((args.max_dp == -1) or (dp <= args.max_dp)) and (vaf >= args.min_vaf): # Discard record if too much depth (most likely collapsed region)

							if (len(rec_vcf.ref) == 1) and (len(rec_vcf.alts[0]) == 1): # Record is a SNP

								if (gq >= args.min_gq_snp): seq[rec_vcf.pos + delta_pos - 1] = rec_vcf.alts[0]

							elif (gq >= args.min_gq_indel): # Record is an indel so check GQ now

								if (len(rec_vcf.ref) > 1) and (len(rec_vcf.alts[0]) == 1): # Record is a deletion

									del seq[rec_vcf.pos + delta_pos : rec_vcf.pos + delta_pos + len(rec_vcf.ref) - 1]

									delta_pos -= len(rec_vcf.ref) - 1

								elif (len(rec_vcf.ref) == 1) and (len(rec_vcf.alts[0]) > 1): # Record is an insertion

									seq[rec_vcf.pos + delta_pos : rec_vcf.pos + delta_pos] = list(rec_vcf.alts[0])[1:]
									
									delta_pos += len(rec_vcf.alts[0]) - 1

								

		print(">" + rec_asm.name) # Name
		print(''.join(seq)) # Sequence

		
