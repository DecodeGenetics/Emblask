#!/usr/bin/env python3

import sys
import pysam
import argparse

from intervaltree import Interval, IntervalTree

# Extract reads matching the reference or alternate allele from the primary alignments of a phased BAM file
# IMPORTANT1: the reads must be aligned with soft-clipping to output the original reads (-Y option in minimap2)
# IMPORTANT2: we assume the phasing tool did not synchronize phasing for supplementary alignments, i.e, primary/supplementary alignments of a read can have different HP in different PS
# IMPORTANT3: because of IMPORTANT2, when using include_supplementary, the output BAM file can contain supplementary alignments WITHOUT any primary alignments

def extractReads(bam_fn, ps, mapq, include_supplementary, include_unmapped_mate, only_include_phased, bam_output, threads_bam):

	seq = {} # Reads to output
	bam_out = (bam_output != "")

	in_bamf = pysam.AlignmentFile(bam_fn, "rb", threads=threads_bam)

	if bam_out: out_bamf = pysam.AlignmentFile(bam_output, "wb", template=in_bamf, threads=threads_bam)

	# 1 - A contig can have multiple PS
	contig2pos_d = {}

	for ps_id, ps_coord_count in ps.items():

		contig = ps_id[0]
		coord = ps_coord_count[0]

		if (contig not in contig2pos_d): contig2pos_d[contig] = IntervalTree()

		contig2pos_d[contig][coord[0]:coord[1]] = ps_id

	for ps_id, ps_coord_count in ps.items():

		contig = ps_id[0]
		coord, counts = ps_coord_count

		# Determine which phaseset we want to extract (HP1 or HP2)
		hp = 0

		if (counts[0] >= counts[1]): hp = 1
		else: hp = 2

		# Output reads matching phaseset
		for rec in in_bamf.fetch(contig, coord[0], coord[1]):

			if (rec.query_length != 0):

				is_primary = (not rec.is_supplementary and not rec.is_secondary and not rec.is_unmapped)

				if is_primary or (include_supplementary and rec.is_supplementary):

					if (rec.has_tag("HP") and (rec.get_tag("HP") == hp)):

						# Determine if this read/alignment overlaps this PS or another PS.
						# A read/alignment is attributed to a PS if it overlaps the most that PS compared to the other PS it overlaps
						is_valid_ps = False
						itvs = contig2pos_d[contig][rec.reference_start:rec.reference_end]

						if (len(itvs) > 0):

							if (len(itvs) == 1): is_valid_ps = True # Read overlaps only one PS so it is automatically valid
							else: # Read overlaps multiple PS, only validate alignment if it overlaps this PS more than the others

								ps_overlap = 0
								len_overlap = 0

								for itv in itvs:

									l_len_overlap = max(0, min(itv.end, rec.reference_end) - max(itv.begin, rec.reference_start))

									if (l_len_overlap >= len_overlap):

										len_overlap = l_len_overlap
										ps_overlap = itv.data

								if (ps_overlap == ps_id[1]): is_valid_ps = True

						if is_valid_ps:

							paired = rec.is_read1 or rec.is_read2
							curr = (False, False)

							if (rec.query_name in seq): curr = seq[rec.query_name]

							if bam_out:

								if is_primary and paired: seq[rec.query_name] = (rec.is_read1 or curr[0], rec.is_read2 or curr[1]) # Record that this read is in the haplotype we wanted for that PS

								out_bamf.write(rec)

							elif is_primary and (not paired and not curr[0]) or (rec.is_read1 and not curr[0]) or (rec.is_read2 and not curr[1]): # Record name was unseen before or mate was unseen before

								# Add record name or mate to dictionnary 
								seq[rec.query_name] = (not paired or rec.is_read1 or curr[0], rec.is_read2 or curr[1])

								if rec.query_qualities: print("@" + rec.query_name + "\n" + rec.query_sequence + "\n+\n" + ''.join(map(lambda x: chr(x+33), rec.query_qualities))) # FASTQ output
								else: print(">" + rec.query_name + "\n" + rec.query_sequence) # FASTA output
	in_bamf.close()

	if not only_include_phased or include_unmapped_mate: # Unmapped reads are assumed to be unphased

		in_bamf = pysam.AlignmentFile(bam_fn, "rb", threads=threads_bam)

		for rec in in_bamf.fetch(until_eof=include_unmapped_mate):

			if (rec.query_length != 0):

				is_primary = (not rec.is_supplementary and not rec.is_secondary and not rec.is_unmapped)

				# 1. if extract unphased reads:
				# - Always extract primary alignment
				# - Extract supplementary if required by input CL argument
				# 2. secondary alignments are always unphased so extract them if required by input CL argument
				# 3. unmapped reads are always unphased so extract them if required by input CL argument
				if (not only_include_phased and (is_primary or (include_supplementary and rec.is_supplementary))) or (include_unmapped_mate and rec.is_unmapped and (rec.query_name in seq)):

					is_unphased = (not rec.has_tag("HP")) or (rec.get_tag("HP") == 0)

					#if (not only_include_phased and is_unphased) or (not is_unphased and not is_primary and (rec.query_name in seq)):
					if is_unphased:

						if bam_out: out_bamf.write(rec)
						else:
							paired = rec.is_read1 or rec.is_read2
							curr = (False, False)

							if (rec.query_name in seq): curr = seq[rec.query_name]

							# Record name was unseen before or mate was unseen before
							if (not paired and not curr[0]) or (rec.is_read1 and not curr[0]) or (rec.is_read2 and not curr[1]):

								# Add record name or mate to dictionnary 
								seq[rec.query_name] = (not paired or rec.is_read1 or curr[0], rec.is_read2 or curr[1])

								if rec.query_qualities: print("@" + rec.query_name + "\n" + rec.query_sequence + "\n+\n" + ''.join(map(lambda x: chr(x+33), rec.query_qualities))) # FASTQ output
								else: print(">" + rec.query_name + "\n" + rec.query_sequence) # FASTA output

		in_bamf.close()

	if bam_out: out_bamf.close()

def computeHapIDperContig(vcf_fn, alt, use_indels):

	ps_d = {}

	pys_vcf = pysam.VariantFile(vcf_fn)
	it_pys_vcf = pys_vcf.fetch() # Iterator over alignments

	for rec in it_pys_vcf:

		if (len(rec.alts) == 1) and (use_indels or ((len(rec.ref) == 1) and (len(rec.alts[0]) == 1))): # Only use bi-allelic SNPs by default, maybe bi-allelic indels if required

			rec_filter = rec.filter.keys()

			if (len(rec_filter) == 1) and ((rec_filter[0] == "PASS") or (rec_filter[0] == ".")): # Make sure variant is PASS

				rec_samples = rec.samples.keys()

				if (len(rec_samples) != 1): sys.exit("Cannot use multi-sample VCF") # Does not process multi-sample VCF

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

					if not alt:

						if (gt == (0, 1)): counts = (counts[0]+1, counts[1])
						elif (gt == (1, 0)): counts = (counts[0], counts[1]+1)

					else:
						if (gt == (0, 1)): counts = (counts[0], counts[1]+1)
						elif (gt == (1, 0)): counts = (counts[0]+1, counts[1])

					ps_d[ps_id] = [coord, counts]
					
	return ps_d

def computeHapIDperContigBED(vcf_fn, bed, use_indels):

	# 1 - Extract regions
	coord_phase_d = {}

	bed_fn = open(bed, "r")

	for rec in bed_fn:

		contig, pos_s_str, pos_e_str, phase = rec.split()

		pos_s = int(pos_s_str)
		pos_e = int(pos_e_str)

		if (contig not in coord_phase_d): coord_phase_d[contig] = IntervalTree()
		if (len(list(coord_phase_d[contig][pos_s:pos_e])) != 0): sys.exit("BED file must contain non-overlapping regions. Abort.")

		if (phase=="alt") or (phase=="ALT"): coord_phase_d[contig][pos_s:pos_e] = 1
		elif (phase=="ref") or (phase=="REF"): coord_phase_d[contig][pos_s:pos_e] = 0
		else: sys.exit("Unrecognized phase annotation in BED file: " + phase + ". Must be ALT or REF. Abort.")

	# 2 - Extract variants
	ps_d = {}

	pys_vcf = pysam.VariantFile(vcf_fn)
	it_pys_vcf = pys_vcf.fetch() # Iterator over alignments

	for rec in it_pys_vcf:

		if (len(rec.alts) == 1) and (use_indels or ((len(rec.ref) == 1) and (len(rec.alts[0]) == 1))): # Only use bi-allelic SNPs by default, maybe bi-allelic indels if required

			rec_filter = rec.filter.keys()

			if (len(rec_filter) == 1) and ((rec_filter[0] == "PASS") or (rec_filter[0] == ".")): # Make sure variant is PASS

				rec_samples = rec.samples.keys()

				if (len(rec_samples) != 1): sys.exit("Cannot use multi-sample VCF") # Does not process multi-sample VCF

				sample = rec_samples[0]

				if ("PS" in rec.samples[sample]): # Variant is phased

					rec_ps = rec.samples[sample]["PS"]
					rec_gt = rec.samples[sample]["GT"]

					chrom = rec.chrom
					pos = int(rec.pos) - 1 # 0-based
					extract_phase = list(coord_phase_d[chrom][pos:pos+1])

					if (len(extract_phase) == 1):

						phase = extract_phase[0].data # 0 (REF) or 1 (ALT)

						ps_id = (chrom, rec_ps)

						if (ps_id not in ps_d): ps_d[ps_id] = [[pos, pos+1], [0, 0]]

						coord, counts = ps_d[ps_id]

						coord[0] = min(coord[0], pos)
						coord[1] = max(coord[1], pos+1)

						if phase == 0:

							if (rec_gt == (0, 1)): counts = (counts[0]+1, counts[1])
							elif (rec_gt == (1, 0)): counts = (counts[0], counts[1]+1)

						else:
							if (rec_gt == (0, 1)): counts = (counts[0], counts[1]+1)
							elif (rec_gt == (1, 0)): counts = (counts[0]+1, counts[1])

						ps_d[ps_id] = [coord, counts]

	return ps_d

if __name__ == "__main__":

	# Parse arguments
	parser = argparse.ArgumentParser(prog='extractHapReadsFromBAM', description='Extract reads matching the reference (default) or alternate allele (--alt) from the alignments of a phased BAM file', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-b', '--in_bam', action='store', help='Phased BAM filename', required=True)
	required.add_argument('-v', '--in_vcf', action='store', help='Phased VCF filename (must contain Phase Set IDs in a PS tag)', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-m', '--mapq', action='store', help='Only consider alignments with this minimum quality.', type=int, default=0, required=False)
	optional.add_argument('-a', '--alt', action='store_true', help='Extract reads for the Alternate allele.', default=False, required=False)
	optional.add_argument('-s', '--include_supplementary', action='store_true', help='Output reads from supplementary alignments in addition to primary.', default=False, required=False)
	optional.add_argument('-u', '--include_unmapped_mate', action='store_true', help='Output unmapped mate if its other mapped is 1/ mapped and 2/ ref or alt phasing (using -a or not).', default=False, required=False)
	optional.add_argument('-p', '--only_include_phased', action='store_true', help='Only output phased (haplotagged) reads.', default=False, required=False)
	optional.add_argument('-B', '--bam_output', action='store', help='Write to BAM file instead of FASTA/FASTQ on stdout', default="", required=False)
	optional.add_argument('-r', '--regions_annot_bed', action='store', help='BED file of annotated regions (REF or ALT) to extract.', default="", required=False)
	optional.add_argument('-t', '--threads', action='store', help='Number of threads for compressing/decompressing BAM', type=int, default=1, required=False)
	optional.add_argument('-i', '--use_indels', action='store_true', help='Use phased indels too (SNPs only by default).', default=False, required=False)

	args = parser.parse_args()

	if ((args.regions_annot_bed != "") and args.alt): sys.exit("Cannot use -a/--alt with -r/--regions_annot_bed. Abort.")

	if (args.regions_annot_bed == ""): ps = computeHapIDperContig(args.in_vcf, args.alt, args.use_indels)
	else: ps = computeHapIDperContigBED(args.in_vcf, args.regions_annot_bed, args.use_indels)

	extractReads(args.in_bam, ps, args.mapq, args.include_supplementary, args.include_unmapped_mate, args.only_include_phased, args.bam_output, args.threads)
