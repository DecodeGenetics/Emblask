#!/bin/bash

THREADS=0
LR_BAM=""
SR_BAM=""
OUT_PREFIX=""
COV_PREFIX=""

PRINT_HELP=0
ONLY_OUTPUT_COVERAGE=0
OUT_UNMAPPED=1

SZ_WINDOW=100000 # Size of a window in bp
L_THREADS=4 # Number of cores to use per job
MAX_BQ=40 # Maximum base quality of the input long reads

# Parse options from the command line
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -t|--global-threads) THREADS=$2; shift ;;
        -T|--local-threads) L_THREADS=$2; shift ;;
        -l|--lr-bam) LR_BAM="$2"; shift ;;
	-s|--sr-bam) SR_BAM="$2"; shift ;;
	-o|--out-pref) OUT_PREFIX="$2"; shift ;;
	-w|--window-sz) SZ_WINDOW=$2; shift ;;
	-C|--cov-pref) COV_PREFIX="$2"; shift ;;
	-Q|--max-bq) MAX_BQ=$2; shift ;;
	-c|--out-cov) ONLY_OUTPUT_COVERAGE=1; ;;
	-u|--no-umap) OUT_UNMAPPED=0; ;;
	-h|--help) PRINT_HELP=1; ;;
        *) echo "Unknown argument used: $1"; exit 1 ;;
    esac
    shift
done

if [ ${PRINT_HELP} -eq 1 ]
then

	echo "Usage:"
	echo ""
	echo "bash $0 -t <NB_THREADS> -l <LONG_READ_BAM> -s <SHORT_READ_BAM> -o <OUTPUT_PREFIX>"
	echo ""
	echo "<Mandatory>: "
	echo "-t <NB_THREADS>: Number of threads to use globally."
	echo "-l <LONG_READ_BAM>: BAM file of long reads mapped to the assembly."
	echo "-s <SHORT_READ_BAM>: BAM file of short reads mapped to the assembly."
	echo "- <OUTPUT_PREFIX>: Output prefix."
	echo ""
	echo "[Optional]: "
	echo "-w [WINDOW_SIZE]: Split regions in windows of that size for correction. Default is 100000."
	echo "-T [LOCAL_NB_THREADS]: Number of threads to use per region correction. Default is 4."
	echo "-C [COVERAGE_FILES_PREFIX]: Load coverage files PREFIX.dict, PREFIX.bed, PREFIX.hc.bed, PREFIX.lmc.bed"
	echo "-Q [MAX_BASE_QUALITY]: Maximum base quality of the input long reads. Default is 40"
	echo "-u: Do not output unmapped (uncorrected) reads."
	echo "-c: Only output coverage files, do not perform correction."
	echo ""

	exit 0
fi

# Arguments check
if [ ${THREADS} -le 0 ]
then
	1>&2 echo "Minimum number of threads cannot be 0 or less. Abort."
	exit 1

elif [ ${L_THREADS} -le 0 ]
then
	1>&2 echo "Minimum number of local threads cannot be 0 or less. Abort."
	exit 1

elif [ ${L_THREADS} -gt ${THREADS} ]
then
	1>&2 echo "Number of local threads cannot be greater than number of global threads. Abort."
	exit 1

elif [ "${LR_BAM}" == "" ] || [ ! -s "${LR_BAM}" ]
then
	1>&2 echo "Empty LR BAM filename given or file doesn't exist or file is empty. Abort."
	exit 1

elif [ "${SR_BAM}" == "" ] || [ ! -s "${SR_BAM}" ]
then
	1>&2 echo "Empty SR BAM filename given or file doesn't exist or file is empty. Abort."
	exit 1

elif [ ${SZ_WINDOW} -le 0 ]
then
	1>&2 echo "Window size cannot be 0 or less. Abort."
	exit 1

elif [ ${MAX_BQ} -le 0 ]
then
	1>&2 echo "Maximum base quality of the input long reads cannot be 0 or less. Abort."
	exit 1

elif [ ${MAX_BQ} -gt 90 ]
then
	1>&2 echo "Maximum base quality of the input long reads be larger than 90. Abort."
	exit 1

elif [ ! "${COV_PREFIX}" == "" ]
then
	if [ ! -e ${COV_PREFIX}.dict ] || [ ! -e ${COV_PREFIX}.bed ] || [ ! -e ${COV_PREFIX}.hc.bed ] || [ ! -e ${COV_PREFIX}.lmc.bed ]
	then
		1>&2 echo "One or more of the following input coverage files does not exist: ${COV_PREFIX}.dict, ${COV_PREFIX}.bed, ${COV_PREFIX}.hc.bed, ${COV_PREFIX}.lmc.bed. Abort."
		exit 1
	fi

	if [ ! -s "${COV_PREFIX}.dict" ] || [ ! -s "${COV_PREFIX}.bed" ]
	then
		1>&2 echo "One or more of the following input coverage files is empty: ${COV_PREFIX}.dict, ${COV_PREFIX}.bed. Abort."
		exit 1
	fi
fi

NB_JOBS=$((THREADS/L_THREADS)) # Number of jobs

# Make temp directory
TMP_DIR=$(mktemp -d)

# Global variables need to be exported to be accessible within functions called by GNU Parallel
export LR_BAM
export SR_BAM
export TMP_DIR

export SZ_WINDOW
export MAX_BQ

export RATATOSK_K64
export RATATOSK_K96
export RATATOSK_K128

export SAMTOOLS
export BEDTOOLS
export SEQTK
export PIGZ_BIN

export EXTRACT_READS_PY
export GET_SUPP_ALIGN_PY

getCoverageRegions() {

	# Create BED file ref coordinates
	${SAMTOOLS} view -H ${LR_BAM} | awk '{if ($1=="@SQ") {CONTIG=substr($2,4,length($2)-3); LEN=substr($3,4,length($3)-3); print CONTIG "\t0\t" LEN}}' > ${TMP_DIR}/contigs.bed

	# Create dict file ref coordinates
	awk '{print $1 "\t" $3}' ${TMP_DIR}/contigs.bed > ${TMP_DIR}/contigs.dict

	# Compute coverage for all primary and supplementary alignments with good mapq
	${SAMTOOLS} depth -@ ${THREADS} -aa -Q 20 -J -G 3844 ${LR_BAM} > ${TMP_DIR}/cov.tsv

	# Compute mean, low and high coverage
	MEAN_COV=$(awk '{SUM+=$3; COUNT+=1} END {print SUM/COUNT}' ${TMP_DIR}/cov.tsv)
	HIGH_COV=$(echo -e "${MEAN_COV}" | awk '{print 1.4*$1}')

	# Create BED file of low and high coverage regions
	awk -v hicov=${HIGH_COV} '{if ($3>=hicov){print $1 "\t" ($2-1) "\t" $2}}' ${TMP_DIR}/cov.tsv > ${TMP_DIR}/hc.bed
	rm -rf ${TMP_DIR}/cov.tsv

	if [ -s ${TMP_DIR}/hc.bed ]
	then
		# Merge high coverage BED file
		${BEDTOOLS} merge -i ${TMP_DIR}/hc.bed > ${TMP_DIR}/hc.merge.bed

		# Keep only regions which are 1kb long minimum
		awk '{if ($3-$2>=1000){print $0}}' ${TMP_DIR}/hc.merge.bed > ${TMP_DIR}/hc.merge.2.bed
		mv -f ${TMP_DIR}/hc.merge.2.bed ${TMP_DIR}/hc.merge.bed

		if [ -s ${TMP_DIR}/hc.merge.bed ]
		then
			# Merge regions which are within 1kb of each other
			${BEDTOOLS} merge -d 1000 -i ${TMP_DIR}/hc.merge.bed > ${TMP_DIR}/hc.merge.2.bed
			mv -f ${TMP_DIR}/hc.merge.2.bed ${TMP_DIR}/hc.merge.bed

			# Remove HC region if it is less than 3 times the reads N50 in the region
			while read BED_REGION
			do

				CONTIG=$(echo -e "${BED_REGION}" | awk '{print $1}')
				START=$(echo -e "${BED_REGION}" | awk '{print $2}')
				END=$(echo -e "${BED_REGION}" | awk '{print $3}')

				LEN_REGION=$((END-START))

				N50=$(${SAMTOOLS} view -@ $((THREADS-2)) -F 3844 ${LR_BAM} "${CONTIG}:${START}-${END}" | awk '{print length($10)}' | sort -rn | \
					awk '{ sum += $0; print $0, sum }' | tac | awk 'NR==1 { halftot=$2/2 } lastsize>halftot && $2<halftot { print $1 } { lastsize=$2 }')
				N50_DIV_3=$((N50/3))

				if [ ${LEN_REGION} -gt ${N50_DIV_3} ]
				then
					echo -e "${BED_REGION}"
				fi 	

			done < ${TMP_DIR}/hc.merge.bed > ${TMP_DIR}/hc.merge.2.bed

			mv -f ${TMP_DIR}/hc.merge.2.bed ${TMP_DIR}/hc.merge.bed

			# Extend high coverage regions of 3kb
			if [ -s ${TMP_DIR}/hc.merge.bed ]
			then
				${BEDTOOLS} slop -i ${TMP_DIR}/hc.merge.bed -g ${TMP_DIR}/contigs.dict -b 3000 | sort -k1,1 -k2,2n > ${TMP_DIR}/hc.merge.extend.bed
				${BEDTOOLS} merge -i ${TMP_DIR}/hc.merge.extend.bed | awk '{if ($2<=5000) {print $1 "\t0\t" $3} else {print $0}}' > ${TMP_DIR}/hc.merge.bed # Regions starting near the contig start are extended to the contig start
				rm -rf ${TMP_DIR}/hc.merge.extend.bed
			
				# Get regions which are not high coverage
				${BEDTOOLS} subtract -a ${TMP_DIR}/contigs.bed -b ${TMP_DIR}/hc.merge.bed > ${TMP_DIR}/lmc.merge.bed

				# Split into sub-regions for multi-threading
				awk -v window_sz=${SZ_WINDOW} '{for (i=$2; i<$3; i+=window_sz) {MAX=$3; if ((i+window_sz)<$3){MAX=i+window_sz}; print $1 "\t" i "\t" MAX}}' ${TMP_DIR}/hc.merge.bed > ${TMP_DIR}/hc.merge.split.bed
				mv -f ${TMP_DIR}/hc.merge.split.bed ${TMP_DIR}/hc.merge.bed
				awk -v window_sz=${SZ_WINDOW} '{for (i=$2; i<$3; i+=window_sz) {MAX=$3; if ((i+window_sz)<$3){MAX=i+window_sz}; print $1 "\t" i "\t" MAX}}' ${TMP_DIR}/lmc.merge.bed > ${TMP_DIR}/lmc.merge.split.bed
				mv -f ${TMP_DIR}/lmc.merge.split.bed ${TMP_DIR}/lmc.merge.bed
			else
				awk -v window_sz=${SZ_WINDOW} '{for (i=$2; i<$3; i+=window_sz) {MAX=$3; if ((i+window_sz)<$3){MAX=i+window_sz}; print $1 "\t" i "\t" MAX}}' ${TMP_DIR}/contigs.bed > ${TMP_DIR}/lmc.merge.bed
			fi
		else
			awk -v window_sz=${SZ_WINDOW} '{for (i=$2; i<$3; i+=window_sz) {MAX=$3; if ((i+window_sz)<$3){MAX=i+window_sz}; print $1 "\t" i "\t" MAX}}' ${TMP_DIR}/contigs.bed > ${TMP_DIR}/lmc.merge.bed
		fi
	else
		awk -v window_sz=${SZ_WINDOW} '{for (i=$2; i<$3; i+=window_sz) {MAX=$3; if ((i+window_sz)<$3){MAX=i+window_sz}; print $1 "\t" i "\t" MAX}}' ${TMP_DIR}/contigs.bed > ${TMP_DIR}/lmc.merge.bed
	fi
}

correctContigSLR() {

	CONTIG=$1
	START=$2
	END=$3

	L_THREADS=$4

	for LR_START_0 in $(seq ${START} ${SZ_WINDOW} $((END-1)))
	do

		LR_END_0=$((LR_START_0+SZ_WINDOW)) # 0-based (BED)
		LR_END_0=$((LR_END_0>END ? END : LR_END_0)) # 0-based (BED)

		LR_START_1=$((LR_START_0+1)) # 1-based (SAM)
		LR_END_1=$((LR_END_0+1)) # 1-based (SAM)

		ID_FILE="${CONTIG}_${LR_START_0}_${LR_END_0}"
		LR_FILE_PREFIX="${TMP_DIR}/lr.${ID_FILE}.corrected"

		echo "Processing ${CONTIG}:${LR_START_1}-${LR_END_1}"

		LR_COV=$(${SAMTOOLS} coverage -r "${CONTIG}:${LR_START_1}-${LR_END_1}" ${LR_BAM} | tail -n 1 | awk '{printf("%.0f\n", $7)}')

		# Extract LR
		${SAMTOOLS} view -b -@ ${L_THREADS} ${LR_BAM} "${CONTIG}:${LR_START_1}-${LR_END_1}" | ${SAMTOOLS} bam2fq -n -@ ${L_THREADS} - 2> /dev/null | ${PIGZ_BIN} -p ${L_THREADS} -c > ${TMP_DIR}/lr.${ID_FILE}.fastq.gz

		if [ ${LR_COV} -gt 500 ] # Coverage is greater than 500
		then
			${PIGZ_BIN} -p ${L_THREADS} -c -d ${TMP_DIR}/lr.${ID_FILE}.fastq.gz > ${LR_FILE_PREFIX}.fastq

		elif [ ! $(zcat ${TMP_DIR}/lr.${ID_FILE}.fastq.gz | head -n 1 | wc -c) -eq 0 ] # LR file is not empty
		then
		
			${PIGZ_BIN} -p ${L_THREADS} -c -d ${TMP_DIR}/lr.${ID_FILE}.fastq.gz | awk '{if (NR%4==1){print substr($0,2,length($0)-1)}}' | sort  > ${TMP_DIR}/lr.${ID_FILE}.rnames

			# Extract BAM of LR overlapping left and right boundaries of region
			${SAMTOOLS} view -b -@ ${L_THREADS} -F 3844 ${LR_BAM} "${CONTIG}:${LR_START_0}-${LR_START_1}" > ${TMP_DIR}/lr.${ID_FILE}.lo.bam # Left-overlapping reads
			${SAMTOOLS} view -b -@ ${L_THREADS} -F 3844 ${LR_BAM} "${CONTIG}:${LR_END_0}-${LR_END_1}" > ${TMP_DIR}/lr.${ID_FILE}.ro.bam # Right-overlapping reads

			# Compute length of longest segment overlapping the left and right boundaries of region
			MAX_LEN_READ_LBORDER=$(${SAMTOOLS} view -@ ${L_THREADS} ${TMP_DIR}/lr.${ID_FILE}.lo.bam | awk -v P_START=${LR_START_1} 'BEGIN {MAXL=0} {L=P_START-$4; if (L>MAXL) {MAXL=L}} END {printf "%3.0f\n", MAXL*1.1}')
			MAX_LEN_READ_RBORDER=$(${SAMTOOLS} view -@ ${L_THREADS} ${TMP_DIR}/lr.${ID_FILE}.ro.bam | awk -v P_END=${LR_END_1} 'BEGIN {MAXL=0} {L=$4+length($10)-P_END; if (L>MAXL) {MAXL=L}} END {printf "%3.0f\n", MAXL*1.1}')

			# Compute list of read names for alignments overlapping left boundary
			#${SAMTOOLS} bam2fq -n -@ ${L_THREADS} ${TMP_DIR}/lr.${ID_FILE}.lo.bam 2> /dev/null | awk '{if (NR%4==1){print substr($0,2,length($0)-1)}}' | sort > ${TMP_DIR}/lr.${ID_FILE}.lo.rnames
			#comm -23 ${TMP_DIR}/lr.${ID_FILE}.rnames ${TMP_DIR}/lr.${ID_FILE}.lo.rnames > ${TMP_DIR}/lr.${ID_FILE}.diff.rnames # Create list of reads not already corrected

			# Compute list of read names for alignments starting within the region
			if [ ${LR_START_0} -eq 0 ]
			then
				cp ${TMP_DIR}/lr.${ID_FILE}.rnames ${TMP_DIR}/lr.${ID_FILE}.diff.rnames # Left boundary position is 0 -> start of the contig, no reads to discard
			else
				${SAMTOOLS} view -@ ${L_THREADS} ${TMP_DIR}/lr.${ID_FILE}.lo.bam | awk -v P_START=${LR_START_1} '{if ($4<P_START) {print $1}}' | sort > ${TMP_DIR}/lr.${ID_FILE}.lo.rnames
				comm -23 ${TMP_DIR}/lr.${ID_FILE}.rnames ${TMP_DIR}/lr.${ID_FILE}.lo.rnames > ${TMP_DIR}/lr.${ID_FILE}.diff.rnames
			fi

			# Clean-up
			rm -rf ${TMP_DIR}/lr.${ID_FILE}.*.bam ${TMP_DIR}/lr.${ID_FILE}.lo.rnames

			# SR base coordinates
			SR_END_0=$((LR_END_0+MAX_LEN_READ_RBORDER))

			SR_START_0=$((LR_START_0-MAX_LEN_READ_LBORDER))
			SR_START_0=$((SR_START_0>=0 ? SR_START_0 : 0))

			# Get SR extended coordinates from secondary alignment coordinates in LR
			echo -e "${CONTIG}\t${SR_START_0}\t${SR_END_0}" > ${TMP_DIR}/sr.${ID_FILE}.bed
			python ${GET_SUPP_ALIGN_PY} -t ${L_THREADS} -b ${LR_BAM} -r ${TMP_DIR}/lr.${ID_FILE}.rnames -s ${CONTIG} >> ${TMP_DIR}/sr.${ID_FILE}.bed
			${BEDTOOLS} slop -b 1000 -i ${TMP_DIR}/sr.${ID_FILE}.bed -g ${TMP_DIR}/contigs.dict | sort -k1,1 -k2,2n > ${TMP_DIR}/sr.${ID_FILE}.2.bed
			${BEDTOOLS} merge -d 1000 -i ${TMP_DIR}/sr.${ID_FILE}.2.bed > ${TMP_DIR}/sr.${ID_FILE}.bed
			rm -rf ${TMP_DIR}/sr.${ID_FILE}.2.bed
			
			# Extract main BAMs
			${SAMTOOLS} view -@ ${L_THREADS} -b -M -L ${TMP_DIR}/sr.${ID_FILE}.bed -q 1 -o ${TMP_DIR}/sr.${ID_FILE}.1.bam -U ${TMP_DIR}/sr.${ID_FILE}.2.bam ${SR_BAM}
			${SAMTOOLS} index -@ ${L_THREADS} ${TMP_DIR}/sr.${ID_FILE}.1.bam
			${SAMTOOLS} index -@ ${L_THREADS} ${TMP_DIR}/sr.${ID_FILE}.2.bam
			# Extract FASTQs
			python ${EXTRACT_READS_PY} -t ${L_THREADS} -b ${TMP_DIR}/sr.${ID_FILE}.1.bam -m 1 --include_supplementary --include_secondary --include_unmapped_mate | ${PIGZ_BIN} -p ${L_THREADS} -c > ${TMP_DIR}/sr.${ID_FILE}.1.fastq.gz
			python ${EXTRACT_READS_PY} -t ${L_THREADS} -b ${TMP_DIR}/sr.${ID_FILE}.2.bam -m 0 --include_supplementary --include_secondary --include_unmapped_mate | ${PIGZ_BIN} -p ${L_THREADS} -c > ${TMP_DIR}/sr.${ID_FILE}.2.fastq.gz

			#${SAMTOOLS} view -@ ${L_THREADS} -b -M -L ${TMP_DIR}/sr.${ID_FILE}.bed ${SR_BAM} -o ${TMP_DIR}/sr.${ID_FILE}.1.bam
			#${SAMTOOLS} index -@ ${L_THREADS} ${TMP_DIR}/sr.${ID_FILE}.1.bam
			#python ${EXTRACT_READS_PY} -t ${L_THREADS} -b ${TMP_DIR}/sr.${ID_FILE}.1.bam -m 0 --include_supplementary --include_secondary --include_unmapped_mate | ${PIGZ_BIN} -p ${L_THREADS} -c > ${TMP_DIR}/sr.${ID_FILE}.1.fastq.gz

			# Clean-up
			rm -rf ${TMP_DIR}/sr.${ID_FILE}.*.bam* ${TMP_DIR}/sr.${ID_FILE}.*bed ${TMP_DIR}/lr.${ID_FILE}.rnames

			HELPER_SR_RATATOSK="-u ${TMP_DIR}/sr.${ID_FILE}.2.fastq.gz"

			if [ ! $(zcat ${TMP_DIR}/sr.${ID_FILE}.1.fastq.gz | head -n 1 | wc -c) -eq 0 ] # SR file 1 is not empty
			then

				if [ $(zcat ${TMP_DIR}/sr.${ID_FILE}.2.fastq.gz | head -n 1 | wc -c) -eq 0 ] # SR file 2 is empty
				then
					HELPER_SR_RATATOSK=""
				fi

				${PIGZ_BIN} -p ${L_THREADS} -c -d ${TMP_DIR}/lr.${ID_FILE}.fastq.gz > ${LR_FILE_PREFIX}.3.fastq

				# Correction SR 1
				${RATATOSK_K64} correct -c ${L_THREADS} -k 31 -K 63 -w 5000 -Q ${MAX_BQ} --force-io-order \
				-s ${TMP_DIR}/sr.${ID_FILE}.1.fastq.gz ${HELPER_SR_RATATOSK} -l ${LR_FILE_PREFIX}.3.fastq -o ${LR_FILE_PREFIX}
				mv -f ${LR_FILE_PREFIX}.fastq ${LR_FILE_PREFIX}.3.fastq

				# Correction SR 2
				${RATATOSK_K96} correct -c ${L_THREADS} -k 63 -K 95 -w 5000 -Q ${MAX_BQ} --force-io-order \
				-s ${TMP_DIR}/sr.${ID_FILE}.1.fastq.gz ${HELPER_SR_RATATOSK} -l ${LR_FILE_PREFIX}.3.fastq -o ${LR_FILE_PREFIX}
				mv -f ${LR_FILE_PREFIX}.fastq ${LR_FILE_PREFIX}.3.fastq

				# Correction LR 1
				${RATATOSK_K128} correct -2 -c ${L_THREADS} -k 125 -K 127 -w 25000 -Q ${MAX_BQ} -i 20000 --no-snp-correction --force-io-order \
				-s ${LR_FILE_PREFIX}.3.fastq -l ${LR_FILE_PREFIX}.3.fastq -L ${TMP_DIR}/lr.${ID_FILE}.fastq.gz -o ${LR_FILE_PREFIX}
				mv -f ${LR_FILE_PREFIX}.fastq ${LR_FILE_PREFIX}.3.fastq

				${SEQTK} subseq ${LR_FILE_PREFIX}.3.fastq ${TMP_DIR}/lr.${ID_FILE}.diff.rnames > ${LR_FILE_PREFIX}.fastq # Only keep reads starting within the region

				# Clean up
				rm -rf ${TMP_DIR}/sr.${ID_FILE}.*fastq.gz ${LR_FILE_PREFIX}.trim.fastq.gz ${LR_FILE_PREFIX}.3.fastq
				rm -rf ${LR_FILE_PREFIX}.index.k*.fasta
			else
				${SEQTK} subseq ${TMP_DIR}/lr.${ID_FILE}.fastq.gz ${TMP_DIR}/lr.${ID_FILE}.diff.rnames > ${LR_FILE_PREFIX}.fastq # Only keep reads starting within the region
			fi

			rm -rf ${TMP_DIR}/lr.${ID_FILE}.diff.rnames ${TMP_DIR}/lr.${ID_FILE}.rnames
		fi

		# Clean-up
		rm -rf ${TMP_DIR}/lr.${ID_FILE}.fastq.gz
	done
}

concatReads() {

	BED_FILE_ALL=$1
	BED_FILE_SUB=$2
	OUT_FILE_PREFIX=$3

	# Make sure concat file does not already exist because of a previously failed run
	while read BED_REGION_ALL
	do
		CONTIG_ALL=$(echo -e "${BED_REGION_ALL}" | awk '{print $1}')
		START_0_ALL=$(echo -e "${BED_REGION_ALL}" | awk '{print $2}')
		END_0_ALL=$(echo -e "${BED_REGION_ALL}" | awk '{print $3}')

		ID_FILE_ALL="${CONTIG_ALL}_${START_0_ALL}_${END_0_ALL}"
		CONCAT_LR_PREFIX="${TMP_DIR}/${OUT_FILE_PREFIX}.${ID_FILE_ALL}"

		rm -rf "${CONCAT_LR_PREFIX}.fastq"
		rm -rf "${CONCAT_LR_PREFIX}.bed"

		${BEDTOOLS} intersect -wb -a <(echo -e "${CONTIG_ALL}\t${START_0_ALL}\t${END_0_ALL}") -b ${BED_FILE_SUB} > "${CONCAT_LR_PREFIX}.bed"

		while read BED_REGION_SUB
		do
			CONTIG_SUB=$(echo -e "${BED_REGION_SUB}" | awk '{print $1}')
			START_0_SUB=$(echo -e "${BED_REGION_SUB}" | awk '{print $2}')
			END_0_SUB=$(echo -e "${BED_REGION_SUB}" | awk '{print $3}')

			ID_FILE_SUB="${CONTIG_SUB}_${START_0_SUB}_${END_0_SUB}"
			LR_FILE_PREFIX="${TMP_DIR}/lr.${ID_FILE_SUB}"

			if [ -s "${LR_FILE_PREFIX}.corrected.fastq" ] # File exists and is not empty
			then
				awk '{if (NR%4==1){print substr($0,2,length($0)-1)}}' "${LR_FILE_PREFIX}.corrected.fastq" | sort >> "${CONCAT_LR_PREFIX}.rnames"

				if [ -f "${CONCAT_LR_PREFIX}.fastq" ]
				then
					cat "${LR_FILE_PREFIX}.corrected.fastq" >> "${CONCAT_LR_PREFIX}.fastq"
				else
					mv -f "${LR_FILE_PREFIX}.corrected.fastq" "${CONCAT_LR_PREFIX}.fastq"
				fi
			fi

		done < "${CONCAT_LR_PREFIX}.bed"

		rm -rf "${CONCAT_LR_PREFIX}.bed"

		if [ -s "${CONCAT_LR_PREFIX}.rnames" ] # Sort list of read names
		then
			sort "${CONCAT_LR_PREFIX}.rnames" > "${CONCAT_LR_PREFIX}.sorted.rnames"
			mv -f "${CONCAT_LR_PREFIX}.sorted.rnames" "${CONCAT_LR_PREFIX}.rnames"
		fi
		
	done < ${BED_FILE_ALL}
}

rmDuplicatesAndConcatAll() {

	BED_FILE_ALL=$1
	PREFIX_FILE_HC=$2
	PREFIX_FILE_LMC=$3
	OUT_FILE_PREFIX=$4

	# Make sure concat file does not already exist because of a previously failed run
	while read BED_REGION_ALL
	do
		CONTIG_ALL=$(echo -e "${BED_REGION_ALL}" | awk '{print $1}')
		START_0_ALL=$(echo -e "${BED_REGION_ALL}" | awk '{print $2}')
		END_0_ALL=$(echo -e "${BED_REGION_ALL}" | awk '{print $3}')

		ID_FILE_ALL="${CONTIG_ALL}_${START_0_ALL}_${END_0_ALL}"
		CONCAT_LR_PREFIX="${TMP_DIR}/${OUT_FILE_PREFIX}"

		PREFIX_FILE_HC_ALL="${TMP_DIR}/${PREFIX_FILE_HC}.${ID_FILE_ALL}"
		PREFIX_FILE_LMC_ALL="${TMP_DIR}/${PREFIX_FILE_LMC}.${ID_FILE_ALL}"

		if [ ! -s "${PREFIX_FILE_HC_ALL}.fastq" ] # HC file does not exist or is empty
		then
			if [ -s "${PREFIX_FILE_LMC_ALL}.fastq" ] # LMC file exists and is not empty
			then
				cat "${PREFIX_FILE_LMC_ALL}.fastq" >> "${CONCAT_LR_PREFIX}.fastq"
			fi

		else # HC file exists and is not empty

			if [ ! -s "${PREFIX_FILE_LMC_ALL}.fastq" ] # LMC file does not exist or is empty
			then
				cat "${PREFIX_FILE_HC_ALL}.fastq" >> "${CONCAT_LR_PREFIX}.fastq"
			else

				# Remove from LMC all the LR that were corrected in the HC list
				comm -13 "${PREFIX_FILE_HC_ALL}.rnames" "${PREFIX_FILE_LMC_ALL}.rnames" > "${CONCAT_LR_PREFIX}.diff.rnames"
				${SEQTK} subseq "${PREFIX_FILE_LMC_ALL}.fastq" "${CONCAT_LR_PREFIX}.diff.rnames" >> "${CONCAT_LR_PREFIX}.fastq"

				# Add HC reads to final output
				cat "${PREFIX_FILE_HC_ALL}.fastq" >> "${CONCAT_LR_PREFIX}.fastq"

				# Clean-up
				rm -rf "${CONCAT_LR_PREFIX}.diff.rnames"
			fi
		fi

		# Final clean-up
		rm -rf ${PREFIX_FILE_HC_ALL}.* ${PREFIX_FILE_LMC_ALL}.*

	done < ${BED_FILE_ALL}
}

getUnmappedReads() {

	${SAMTOOLS} view -@ $((THREADS/2)) -b -f 4 ${LR_BAM} | ${SAMTOOLS} bam2fq -n -@ $((THREADS/2)) - 2> /dev/null >> ${OUT_PREFIX}.fastq;
}

# Function names need to be exported to be accessible by GNU Parallel
export -f getCoverageRegions
export -f correctContigSLR
export -f concatReads
export -f rmDuplicatesAndConcatAll
export -f getUnmappedReads

# ============= Coverage ============= #

if [ ! "${COV_PREFIX}" == "" ]
then

	echo -e "- Loading coverage"

	ln -s $(readlink -f ${COV_PREFIX}.dict) ${TMP_DIR}/contigs.dict # Dictionnary of all contigs involved (pair<contig_name,contig_length>)
	ln -s $(readlink -f ${COV_PREFIX}.bed) ${TMP_DIR}/contigs.bed # BED file of merged HC and LMC regions
	ln -s $(readlink -f ${COV_PREFIX}.hc.bed) ${TMP_DIR}/hc.merge.bed # BED file of HC regions
	ln -s $(readlink -f ${COV_PREFIX}.lmc.bed) ${TMP_DIR}/lmc.merge.bed # BED file of LMC regions
else

	echo -e "- Computing coverage"
	getCoverageRegions

	if [ ${ONLY_OUTPUT_COVERAGE} -eq 1 ]
	then 
		mv ${TMP_DIR}/contigs.dict ${OUT_PREFIX}.dict # Dictionnary of all contigs involved (pair<contig_name,contig_length>)
		mv ${TMP_DIR}/contigs.bed ${OUT_PREFIX}.bed # BED file of merged HC and LMC regions
		mv ${TMP_DIR}/hc.merge.bed ${OUT_PREFIX}.hc.bed # BED file of HC regions
		mv ${TMP_DIR}/lmc.merge.bed ${OUT_PREFIX}.lmc.bed # BED file of LMC regions

		# Clean-up
		rm -rf ${TMP_DIR}

		exit 0
	fi
fi
# =========== Correction =========== #

if [ -s ${TMP_DIR}/hc.merge.bed ]
then

	echo -e "- Correcting High Coverage regions in parallel using SR"

	if [ $(wc -l ${TMP_DIR}/hc.merge.bed | awk '{print $1}') -eq 1 ]
	then
		CONTIG=$(cat ${TMP_DIR}/hc.merge.bed | awk '{print $1}')
		START=$(cat ${TMP_DIR}/hc.merge.bed | awk '{print $2}')
		END=$(cat ${TMP_DIR}/hc.merge.bed | awk '{print $3}')

		correctContigSLR ${CONTIG} ${START} ${END} ${THREADS}
	else
		xargs -P ${NB_JOBS} -I {} -a ${TMP_DIR}/hc.merge.bed bash -c \
		"CONTIG=\$(echo \"{}\" | awk '{print \$1}'); START=\$(echo \"{}\" | awk '{print \$2}'); END=\$(echo \"{}\" | awk '{print \$3}'); correctContigSLR \$CONTIG \$START \$END ${L_THREADS}" _
	fi

	echo -e "- Merging High Coverage read partitions"
	concatReads "${TMP_DIR}/contigs.bed" "${TMP_DIR}/hc.merge.bed" "lr.hc"
fi

if [ -s ${TMP_DIR}/lmc.merge.bed ]
then

	echo -e "- Correcting Low-Medium Coverage regions in parallel using SR"

	if [ $(wc -l ${TMP_DIR}/lmc.merge.bed | awk '{print $1}') -eq 1 ]
	then
		CONTIG=$(cat ${TMP_DIR}/lmc.merge.bed | awk '{print $1}')
		START=$(cat ${TMP_DIR}/lmc.merge.bed | awk '{print $2}')
		END=$(cat ${TMP_DIR}/lmc.merge.bed | awk '{print $3}')

		correctContigSLR ${CONTIG} ${START} ${END} ${THREADS}
	else
		xargs -P ${NB_JOBS} -I {} -a ${TMP_DIR}/lmc.merge.bed bash -c \
		"CONTIG=\$(echo \"{}\" | awk '{print \$1}'); START=\$(echo \"{}\" | awk '{print \$2}'); END=\$(echo \"{}\" | awk '{print \$3}'); correctContigSLR \$CONTIG \$START \$END ${L_THREADS}" _
	fi

	echo -e "- Merging Low-Medium read partitions"
	concatReads "${TMP_DIR}/contigs.bed" "${TMP_DIR}/lmc.merge.bed" "lr.lmc"
fi

echo -e "- Final merging of read partitions"
rmDuplicatesAndConcatAll "${TMP_DIR}/contigs.bed" "lr.hc" "lr.lmc" "lr.merged"

# =========== Merging =========== #
mv ${TMP_DIR}/lr.merged.fastq ${OUT_PREFIX}.fastq

if [ ${OUT_UNMAPPED} -eq 1 ]
then
	echo -e "- Adding unmapped (uncorrected) reads"
	getUnmappedReads
fi

# Clean-up
rm -rf ${TMP_DIR}
