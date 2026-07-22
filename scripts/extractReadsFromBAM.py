#!/usr/bin/env python3

import sys
import pysam
import argparse

# Extract all reads from primary, supplementary, secondary and/or unmapped alignment records from a BAM file
# IMPORTANT: the reads must be aligned with soft-clipping to output the original reads (-Y option in minimap2)

def extractReads(bam_fn, tag_values, no_tag, mapq, include_supplementary, include_secondary, include_unmapped_mate, include_unmapped_all, bam_output, include_read_list, exclude_read_list, threads_bam):

	seq = {} # Reads to output
	bam_out = (bam_output != "") # Whether output file is BAM

	rl_d = {} # Read list to include/exclude
	include_rl = False
	exclude_rl = False

	in_bamf = pysam.AlignmentFile(bam_fn, "rb", threads=threads_bam) # Open input BAM file

	if bam_out: out_bamf = pysam.AlignmentFile(bam_output, "wb", template=in_bamf, threads=threads_bam) # Open output BAM file

	if (include_read_list != ""):

		rl_fn = open(include_read_list, "r")
		include_rl = True

		for line in rl_fn: rl_d[line.strip()] = None

		rl_fn.close()

	elif (exclude_read_list != ""):

		rl_fn = open(exclude_read_list, "r")
		exclude_rl = True

		for line in rl_fn: rl_d[line.strip()] = None

		rl_fn.close()

	if (len(tag_values) == 0) and (len(no_tag) == 0): # No tag filter to consider

		for rec in in_bamf.fetch(until_eof=include_unmapped_all):

			# Record is not unmapped or below minimum MAPQ
			if (rec.query_length != 0) and (((not rec.is_unmapped) and (rec.mapping_quality >= mapq)) or (include_unmapped_all and rec.is_unmapped)):

				# Record is not on an include/exclude list or record is on a list and matches a name on that list
				if ((not include_rl) and (not exclude_rl)) or (include_rl and rec.query_name in rl_d) or (exclude_rl and rec.query_name not in rl_d):

					is_primary = (not rec.is_supplementary and not rec.is_secondary and not rec.is_unmapped)

					# Record is either primary or supplementary/secondary/unmapped alignment record if required
					if is_primary or (include_supplementary and rec.is_supplementary) or (include_secondary and rec.is_secondary) or (include_unmapped_all and rec.is_unmapped):

						paired = (rec.is_read1 or rec.is_read2)
						curr = (False, False)

						if (rec.query_name in seq): curr = seq[rec.query_name]

						# Record name was unseen before or mate was unseen before
						if bam_out or (not paired and not curr[0]) or (rec.is_read1 and not curr[0]) or (rec.is_read2 and not curr[1]):

							# Add record name or mate to dictionnary 
							seq[rec.query_name] = (not paired or rec.is_read1 or curr[0], rec.is_read2 or curr[1])

							if bam_out: out_bamf.write(rec) # Write to output BAM
							elif rec.query_qualities: print("@" + rec.query_name + "\n" + rec.query_sequence + "\n+\n" + ''.join(map(lambda x: chr(x+33), rec.query_qualities))) # FASTQ output
							else: print(">" + rec.query_name + "\n" + rec.query_sequence) # FASTA output

	else:

		tag_values_parsed = [] # (key_values)

		for tv in tag_values:

			if (len(tv.split(":")) == 3): #Tag must be of the form TAG_NAME:TAG_VALUE_TYPE:TAG_VALUE

				tag_name, tag_type, tag_value = tv.split(":")

				if (tag_type == 'i'): tag_values_parsed.append((tag_name, int(tag_value)))
				elif (tag_type == 'f'): tag_values_parsed.append((tag_name, float(tag_value)))
				elif (tag_type == 'Z'): tag_values_parsed.append((tag_name, tag_value))
				else: sys.exit("Unspported tag value type: " + tag_type)

			else: sys.exit("Unspported tag: " + tag)

		for rec in in_bamf.fetch(until_eof=include_unmapped_all):

			if (rec.query_length != 0) and (((not rec.is_unmapped) and (rec.mapping_quality >= mapq)) or (include_unmapped_all and rec.is_unmapped)):

				# Record is not on an include/exclude list or record is on a list and matches a name on that list
				if ((not include_rl) and (not exclude_rl)) or (include_rl and rec.query_name in rl_d) or (exclude_rl and rec.query_name not in rl_d):

					is_primary = (not rec.is_supplementary and not rec.is_secondary and not rec.is_unmapped)

					# Record is either primary or supplementary/secondary alignment if required
					if is_primary or (include_supplementary and rec.is_supplementary) or (include_secondary and rec.is_secondary):

						paired = (rec.is_read1 or rec.is_read2)
						curr = (False, False)

						if (rec.query_name in seq): curr = seq[rec.query_name]

						# Record name was unseen before or mate was unseen before
						if bam_out or ((not paired and not curr[0]) or (rec.is_read1 and not curr[0]) or (rec.is_read2 and not curr[1])):

							has_tag_value = False
							has_no_tag = False

							for tag_name, tag_value in tag_values_parsed: # Check if read has any of the tagss

								if (rec.has_tag(tag_name) and (rec.get_tag(tag_name) == tag_value)):

									has_tag_value = True
									break

							if (not has_tag_value) and (not rec.is_supplementary) and (not rec.is_secondary):

								for tag_name in no_tag:

									if not rec.has_tag(tag_name):

										has_no_tag = True
										break		

							if has_tag_value or has_no_tag:

								# Add record name or mate to dictionnary 
								seq[rec.query_name] = (not paired or rec.is_read1 or curr[0], rec.is_read2 or curr[1])

								if bam_out: out_bamf.write(rec) # Write to output BAM
								elif rec.query_qualities: print("@" + rec.query_name + "\n" + rec.query_sequence + "\n+\n" + ''.join(map(lambda x: chr(x+33), rec.query_qualities))) # FASTQ output
								else: print(">" + rec.query_name + "\n" + rec.query_sequence) # FASTA output

					elif (include_unmapped_all and rec.is_unmapped): # Record is unmapped

						if bam_out: out_bamf.write(rec) # Write to output BAM
						elif rec.query_qualities: print("@" + rec.query_name + "\n" + rec.query_sequence + "\n+\n" + ''.join(map(lambda x: chr(x+33), rec.query_qualities))) # FASTQ output
						else: print(">" + rec.query_name + "\n" + rec.query_sequence) # FASTA output

	in_bamf.close()

	if (include_unmapped_mate): # Get unmapped reads for which the other mate of the pair has a mapped primary alignment

		in_bamf = pysam.AlignmentFile(bam_fn, "rb", threads=threads_bam)

		for rec in in_bamf.fetch(until_eof=True):

			# Record is unmapped, not null or below minimum MAPQ
			if rec.is_unmapped and (rec.query_length != 0) and (rec.query_name in seq):

				# Record is not on an include/exclude list or record is on a list and matches a name on that list
				if ((not include_rl) and (not exclude_rl)) or (include_rl and rec.query_name in rl_d) or (exclude_rl and rec.query_name not in rl_d):

					if bam_out: out_bamf.write(rec) # Write to output BAM
					elif rec.query_qualities: print("@" + rec.query_name + "\n" + rec.query_sequence + "\n+\n" + ''.join(map(lambda x: chr(x+33), rec.query_qualities))) # FASTQ output
					else: print(">" + rec.query_name + "\n" + rec.query_sequence) # FASTA output

		in_bamf.close()

	if bam_out: out_bamf.close()

if __name__ == "__main__":

	# Parse arguments
	parser = argparse.ArgumentParser(prog='extractReadsFromBAM', description='Extract all reads from primary, supplementary and secondary alignments from a BAM file', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	required = parser.add_argument_group('Required arguments')

	required.add_argument('-b', '--in_bam', action='store', help='BAM filename', required=True)

	optional = parser.add_argument_group('Optional arguments')

	optional.add_argument('-l', '--include_read_list', action='store', help='Name of reads to include', default="", required=False)
	optional.add_argument('-L', '--exclude_read_list', action='store', help='Name of reads to exclude', default="", required=False)
	optional.add_argument('-i', '--include_tag_value', action='append', help='Only output reads for which at least one alignment has this tag name and value. Must be "name:type:value".', default=[], required=False)
	optional.add_argument('-I', '--include_no_tag', action='append', help='Only output reads for which at least one alignment does NOT have this tag name', default=[], required=False)
	optional.add_argument('-m', '--mapq', action='store', help='Only consider alignments with this minimum quality.', type=int, default=0, required=False)
	optional.add_argument('-s', '--include_supplementary', action='store_true', help='Output reads from supplementary alignments in addition to primary.', default=False, required=False)
	optional.add_argument('-S', '--include_secondary', action='store_true', help='Output reads from secondary alignments in addition to primary.', default=False, required=False)
	optional.add_argument('-u', '--include_unmapped_mate', action='store_true', help='Output unmapped mate if its other mapped is mapped (and included one tag value of -i if was used).', default=False, required=False)
	optional.add_argument('-U', '--include_unmapped_all', action='store_true', help='Output all unmapped reads (regardless of tag values if -i was used).', default=False, required=False)
	optional.add_argument('-B', '--bam_output', action='store', help='Write to BAM file instead of FASTA/FASTQ on stdout', default="", required=False)
	optional.add_argument('-t', '--threads', action='store', help='Number of threads for compressing/decompressing BAM', type=int, default=1, required=False)

	args = parser.parse_args()

	if ((args.include_read_list != "") and (args.exclude_read_list != "")): sys.exit("Cannot use -l and -L (--include_read_list and --exclude_read_list) at the same time.")
	if (args.include_unmapped_mate and args.include_unmapped_all): sys.exit("Cannot use -u and -U (--include_unmapped_mate and --include_unmapped_all) at the same time.")

	# Do what you got to do
	extractReads(args.in_bam, args.include_tag_value, args.include_no_tag, args.mapq, args.include_supplementary, args.include_secondary, args.include_unmapped_mate, args.include_unmapped_all, args.bam_output, args.include_read_list, args.exclude_read_list, args.threads)
