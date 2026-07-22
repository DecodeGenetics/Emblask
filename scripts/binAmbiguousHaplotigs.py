#!/usr/bin/env python3

import sys
import pysam
import argparse

from intervaltree import Interval, IntervalTree

# Extract reads matching the reference or alternate allele from the primary alignments of a phased BAM file
# IMPORTANT1: the reads must be aligned with soft-clipping to output the original reads (-Y option in minimap2)
# IMPORTANT2: we assume the phasing tool did not synchronize phasing for supplementary alignments, i.e, primary/supplementary alignments of a read can have different HP in different PS
# IMPORTANT3: because of IMPORTANT2, when using include_supplementary, the output BAM file can contain supplementary alignments WITHOUT any primary alignments

def loadBinning(bin_list_fn):

	bin_d = {}
	f = open(bin_list_fn, "r")

	for line in f:

		hap, binning = line.split()
	
		if (hap in bin_d): sys.exit("Haplotig/Contig " + hap + " is duplicated in input binning list")

		bin_d[hap] = int(binning)

	f.close()

	return bin_d

def binHaplotig(bin_d, bam_fn, ps, min_mapq, threads_bam):

	level_certainty = 0
	stable_bin = False

	new_bin_d = {k:(v,level_certainty) for k,v in bin_d.items()}

	in_bamf = pysam.AlignmentFile(bam_fn, "rb", threads=threads_bam)

	# 1 - Extract name of reads mapping to homref phase set(s) and unphased reads
	hap2ps_coord_d = {}

	for ps_id, ps_coord_count in ps.items():

		hapn = ps_id[0]
		coord = ps_coord_count[0]

		if (hapn not in hap2ps_coord_d): hap2ps_coord_d[hapn] = IntervalTree()

		hap2ps_coord_d[hapn][coord[0]:coord[1]] = ps_id

	while not stable_bin:

		# List of ambiguous haplotigs to dictionnary
		hapn_d = {}

		for hap, binning in new_bin_d.items():

			if (binning[0] == 0): hapn_d[hap] = None

		# 2 - Iterate over ambiguous or undetermined binning haplotigs, see where the supplementary alignments of their reads map to
		hap2sa = {}

		level_certainty += 1

		for hapn in hapn_d.keys():

			hap_len = in_bamf.get_reference_length(hapn)
			coords = ((0,1), (hap_len-1, hap_len)) # First coord of haplotig and last coords of haplotig

			hap2sa[hapn] = [[], []]

			for i in range(0, len(coords)):

				for rec in in_bamf.fetch(hapn, coords[i][0], coords[i][1]):

					is_primary = (not rec.is_supplementary and not rec.is_secondary)

					# Record maps and is primary or supplementary alignment
					if (rec.query_length != 0) and (is_primary or rec.is_supplementary) and not rec.is_unmapped:

						if rec.has_tag("SA"): # Record is not interesting if it doesn't have supplementary alignments

							isValidRec = not rec.has_tag("HP") # Valid record if unphased

							if (not isValidRec) and (hapn in hap2ps_coord_d): # Record is phased, must check if read his homref haplotype

								hp = rec.get_tag("HP")
								ps_itvs = list(hap2ps_coord_d[hapn][rec.reference_start:rec.reference_end])

								if (len(ps_itvs) > 0):

									ps_id = 0

									if (len(ps_itvs) == 1): ps_id = ps_itvs[0].data # Read overlaps only one PS
									else: # Read overlaps multiple PS, only validate alignment if it overlaps this PS more than the others

										len_overlap = 0

										for itv in ps_itvs:

											l_len_overlap = max(0, min(itv.end, rec.reference_end) - max(itv.begin, rec.reference_start))

											if (l_len_overlap >= len_overlap):

												len_overlap = l_len_overlap
												ps_id = itv.data

									ps_counts = ps[ps_id][1]
									
									if (ps_counts[0] >= ps_counts[1]) and (hp == 1): isValidRec = True
									if (ps_counts[0] < ps_counts[1]) and (hp == 2): isValidRec = True

							if isValidRec:

								sa_tags = rec.get_tag("SA")
								sa_best_align = ("", "", 0, 0)

								for sa_tag in sa_tags.split(";")[:-1]:

									ref_n, ref_pos, strand, cigar, mapq, nm = sa_tag.split(",")

									if (ref_n not in hapn_d): # Supp. alignment is to a haplotig which as a confident binning

										mapq = int(mapq)
										ref_pos_s = int(ref_pos)
										ref_pos_e = ref_pos_s

										if (mapq >= min_mapq): # MAPQ is good enough

											op_ll = []

											for c in cigar: # Parse CIGAR to determine end position of alignment on reference

												if c.isnumeric(): op_ll.append(c)
												else:
													if (c == 'M') or (c == 'D'): ref_pos_e += int(''.join(op_ll))

													op_ll = []

											if ((ref_pos_e - ref_pos_s) > (sa_best_align[3] - sa_best_align[2])): sa_best_align = (rec.query_name, ref_n, ref_pos_s, ref_pos_e)

								if (sa_best_align[3] - sa_best_align[2] > 0): hap2sa[hapn][i].append(sa_best_align)

		# 3 - Iterate over SA, see to which haplotype they bin to
		hapn2align_len = {}
		
		for hapn, info_align_l in hap2sa.items():

			hapn2align_len[hapn] = [[], []]

			for i in range(len(hapn2align_len[hapn])):

				tot_align_len = [0, 0] # (H1, H2)

				for sa in info_align_l[i]:

					query_n, ref_n, ref_pos_s, ref_pos_e = sa

					for rec in in_bamf.fetch(ref_n, ref_pos_s, ref_pos_e):

						if (rec.query_name == query_n) and (rec.reference_start == ref_pos_s-1) and (rec.reference_end == ref_pos_e-1):

							is_primary = (not rec.is_supplementary and not rec.is_secondary)

							if (rec.query_length != 0) and (is_primary or rec.is_supplementary) and not rec.is_unmapped:

								ref_align_len = rec.reference_end - rec.reference_start
								ps_itvs = []
								ps_id = 0

								if (ref_n in hap2ps_coord_d): ps_itvs = list(hap2ps_coord_d[ref_n][rec.reference_start:rec.reference_end])

								if (len(ps_itvs) == 0): tot_align_len[new_bin_d[ref_n][0]-1] += ref_align_len # No phasing set here
								elif rec.has_tag("HP"):

									if (len(ps_itvs) == 1): ps_id = ps_itvs[0].data # Read overlaps only one PS
									else: # Read overlaps multiple PS, only validate alignment if it overlaps this PS more than the others

										len_overlap = 0

										for itv in ps_itvs:

											l_len_overlap = max(0, min(itv.end, rec.reference_end) - max(itv.begin, rec.reference_start))

											if (l_len_overlap >= len_overlap):

												len_overlap = l_len_overlap
												ps_id = itv.data

									ps_counts = ps[ps_id][1]
									hp_homref = 0

									if (ps_counts[0] >= ps_counts[1]): hp_homref = 1
									else: hp_homref = 2

									if (rec.get_tag("HP") == hp_homref): tot_align_len[new_bin_d[ref_n][0]-1] += ref_align_len
									elif (new_bin_d[ref_n][0] == 1): tot_align_len[1] += ref_align_len
									else: tot_align_len[0] += ref_align_len

				hapn2align_len[hapn][i] = tot_align_len

		#print(str(hap2sa))
		#print(" ")
		#print(str(hapn2align_len))

		# 4 - TODO
		stable_bin = True

		for hapn, align_hap in hapn2align_len.items():

			left_h = 0
			right_h = 0
			selected_h = 0

			sum_left = align_hap[0][0] + align_hap[0][1]
			sum_right = align_hap[1][0] + align_hap[1][1]

			if (sum_left != 0):

				if ((align_hap[0][0] / sum_left) >= 0.8): left_h = 1
				elif ((align_hap[0][1] / sum_left) >= 0.8): left_h = 2

			if (sum_right != 0):

				if ((align_hap[1][0] / sum_right) >= 0.8): right_h = 1
				elif ((align_hap[1][1] / sum_right) >= 0.8): right_h = 2

			if (left_h == right_h): selected_h = left_h
			elif ((left_h != 0) and (right_h == 0)): selected_h = left_h
			elif ((left_h == 0) and (right_h != 0)): selected_h = right_h

			if (selected_h != new_bin_d[hapn][0]):

				new_bin_d[hapn] = (selected_h, level_certainty)
				stable_bin = False

			elif (selected_h==0): new_bin_d[hapn] = (0, level_certainty)

	for hapn, bin_h in new_bin_d.items(): print(hapn + "\t" + str(bin_h[0]) + "\t" + str(bin_h[1]))

def computeHapIDperContig(vcf_fn, use_indels):

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

					if (gt == (0, 1)): counts = (counts[0]+1, counts[1])
					elif (gt == (1, 0)): counts = (counts[0], counts[1]+1)

					ps_d[ps_id] = [coord, counts]
					
	return ps_d

if __name__ == "__main__":

	# Parse arguments
	parser = argparse.ArgumentParser(prog='binAmbiguousHaplotigs', description='Bin ambiguous haplotigs using the reads overlapping their boundaries and mapping to non-ambiguously binned haplotigs', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-b', '--in_bam', action='store', help='Phased BAM filename', required=True)
	required.add_argument('-v', '--in_vcf', action='store', help='Phased VCF filename (must contain Phase Set IDs in a PS tag)', required=True)
	required.add_argument('-B', '--in_binning', action='store', help='Haplotig binning in TSV. Each row is haplotig name and haplotype ID (0 for undefined or ambiguous, otherwise {1|2})', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-m', '--mapq', action='store', help='Only consider supplementary alignments with this minimum quality.', type=int, default=60, required=False)
	optional.add_argument('-t', '--threads', action='store', help='Number of threads for compressing/decompressing BAM', type=int, default=1, required=False)
	optional.add_argument('-i', '--use_indels', action='store_true', help='Use phased indels too (SNPs only by default).', default=False, required=False)

	args = parser.parse_args()

	ps_d = computeHapIDperContig(args.in_vcf, args.use_indels)
	bin_d = loadBinning(args.in_binning)

	binHaplotig(bin_d, args.in_bam, ps_d, args.mapq, args.threads)
