#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { mapLRtoAsm as mapLRtoAsm_1; mapLRtoAsm as mapLRtoAsm_2 } from 'module/main.nf'
include { var_CallFilterPhase as var_CallFilterPhase_1; var_CallFilterPhase as var_CallFilterPhase_2 } from 'module/main.nf'
include { varCall_hifi as hapResAsmPolish_phase_varCall } from 'module/main.nf'
include { varCall_ontR9_trainedModels as mixHapAsm_filter_call; varCall_ontR9_trainedModels as hapResAsmPolish_polishDualAsm_call } from 'module/main.nf'
include { varPhase as mixHapAsm_filter_phase; varPhase as hapResAsmPolish_phase_varPhase1; varPhase as hapResAsmPolish_phase_varPhase2 } from 'module/main.nf'
include { mapPairedIllumina as mapPairedIllumina_correct; mapPairedIllumina as mapPairedIllumina_polish } from 'module/main.nf'
include { mergePairedIlluminaBAM as mergePairedIlluminaBAM_correct; mergePairedIlluminaBAM as mergePairedIlluminaBAM_polish } from 'module/main.nf'

/*
Basic quality filtering on corrected ONT reads stored in FASTQ or BAM.
By default, the script detects and removes windows >= 100 bp for which >90% of the bases have QUAL<=5.
Hence N50 of the filtered reads is decreased but quality is increased.
Output reads are compressed
*/
process filterONT {

	label 'small_node'

	input:
		path bam_in
		path fastq_in

	output:
		path 'reads.filtered.fq.gz'

	shell '/bin/bash', '-euo', 'pipefail'

	script:

		if (bam_in) {
			"""
			samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
			pigz=\${PIGZ_BIN:-${params.tools.pigz.bin}}

			filter_fastq_py=\${FILTER_FASTQ_PY:-${params.python.script.filter_fastq}}

			MAX_BQ=\$(bc -l <<< \"${params.pipeline.fastq_filter.min_baseq_ratio} * ${params.max_lr_bq}\" | awk '{printf(\"%.0f\", \$0)}') 

			\${samtools} bam2fq -n -@ ${task.cpus} ${bam_in} | \${pigz} -p ${task.cpus} -c > reads.fq.gz;
			python \${filter_fastq_py} -f reads.fq.gz -t ${task.cpus} -c 1000 -l ${params.pipeline.min_len_read} -b \${MAX_BQ} \
			-w ${params.pipeline.fastq_filter.len_window_baseq_filter} -r ${params.pipeline.fastq_filter.ratio_low_baseq_window_filter} | \
			\${pigz} -p ${task.cpus} -c > reads.filtered.fq.gz;
			rm -rf reads.fq.gz;
			"""
		}
		else if (fastq_in) {
			"""
			pigz=\${PIGZ_BIN:-${params.tools.pigz.bin}}

			filter_fastq_py=\${FILTER_FASTQ_PY:-${params.python.script.filter_fastq}}

			MAX_BQ=\$(bc -l <<< \"${params.pipeline.fastq_filter.min_baseq_ratio} * ${params.max_lr_bq}\" | awk '{printf(\"%.0f\", \$0)}') 

			python \${filter_fastq_py} -f ${fastq_in} -t ${task.cpus} -c 1000 -l ${params.pipeline.min_len_read} -b \${MAX_BQ} \
			-w ${params.pipeline.fastq_filter.len_window_baseq_filter} -r ${params.pipeline.fastq_filter.ratio_low_baseq_window_filter} | \
			\${pigz} -p ${task.cpus} -c > reads.filtered.fq.gz;
			"""
		}
		else error("No corrected reads in BAM or FASTQ format provided in input")
}

/*
Extract Illumina reads from a BAM file (not in pairs).
Build k=31 de Bruijn graph and index from Illumina reads.
Compress output.
*/
process extractIllumina_buildDBG {

	label 'medium_node'

	input:
		tuple val(trio_member), path(bam_in), path(fastq_in) // trio_member can be 'proband', 'father' or 'mother'

	output:
		tuple path("${trio_member}.fq.gz"), path("${trio_member}.gfa.gz"), path("${trio_member}.bfi")

	shell '/bin/bash', '-euo', 'pipefail'

	script:

		if (bam_in) {
			"""
			samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
			bifrost=\${BIFROST:-${params.tools.bifrost.bin}}

			\${samtools} fastq -@ ${task.cpus} -o ${trio_member}.fq.gz ${bam_in};
			\${bifrost} build -v -k 31 -t ${task.cpus} -s ${trio_member}.fq.gz -o ${trio_member};
			"""
		}
		else if (fastq_in) {
			"""
			bifrost=\${BIFROST:-${params.tools.bifrost.bin}}

			ln -s ${fastq_in} ${trio_member}.fq.gz
			\${bifrost} build -v -k 31 -t ${task.cpus} -s ${trio_member}.fq.gz -o ${trio_member};
			"""
		}
		else error("No corrected reads in BAM or FASTQ format provided in input")
}

/*
Extract interleaved Illumina reads from a BAM file.
Compress output.
*/
process extractPairedIllumina {

	label 'medium_node'

	input:
		path(sr_bam_in)

	output:
		path "sr.fq.gz"

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}

	\${samtools} collate -@ ${task.cpus} -u -O -n 128 ${sr_bam_in} /tmp/collate | \${samtools} fastq -o sr.fq.gz -;
	"""
}

/*
Split a compressed FASTQ file of Illumina reads into uncompressed chunks of equal number of lines.
*/
process splitPairedIllumina {

	label 'medium_node'

	input:
		path sr_fq

	output:
		path "sr.fq.part_*"

	shell '/bin/bash', '-euo', 'pipefail'

	// Split paired short read file into chunks of 80M lines (or 10M pairs)
	"""
	pigz=\${PIGZ_BIN:-${params.tools.pigz.bin}}

	\${pigz} -p ${task.cpus} -d -c ${sr_fq} | split -l 80000000 -d - sr.fq.part_;
	"""
}

process mixHapAsm {

	label 'large_node'

	input:
		path filtered_lr_fq

	output:
		tuple path('assembly.fasta'), path('assembly.fasta.fai'), path('assembly.bed')

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	flye=\${FLYE:-${params.tools.flye.bin}}
	seqkit=\${SEQKIT:-${params.tools.seqkit.bin}}

	FLYE_SUBSAMPLE_PARAM=\"\";
	if [ ${params.genome_size} -gt 0 ]; then FLYE_SUBSAMPLE_PARAM=\"--asm-coverage ${params.pipeline.assembly.subsampling_cov} --genome-size ${params.genome_size}\"; fi;

	\${flye} --nano-corr ${filtered_lr_fq} -o . -t ${task.cpus} \${FLYE_SUBSAMPLE_PARAM} --read-error ${params.pipeline.assembly.max_error_rate.lenient} --scaffold --extra-params max_bubble_length=300000;
	if [ ! \$? -eq 0 ]; then echo \"Flye has exited prematurely\" 1>&2; exit 1; fi;
	ASM_COMPLETED=\$({ grep \"INFO: Final assembly\" <(tail -n 50 flye.log) || true; } | wc -l);
	if [ \${ASM_COMPLETED} -eq 0 ] || [ ! -s assembly.fasta ]; then echo \"Flye log incomplete or empty assembly\" 1>&2; exit 1; fi;

	\${seqkit} seq -w0 assembly.fasta > assembly.tmp.fasta
	mv -f assembly.tmp.fasta assembly.fasta
	\${samtools} faidx assembly.fasta
	awk '{print \$1 \"\\t0\\t\" \$2}' assembly.fasta.fai | sort -k1,1 -k2,2n > assembly.bed

	rm -rf 00-assembly 10-consensus 20-repeat 30-contigger 40-polishing;
	"""
}

process mixHapAsm_filter {

	label 'medium_node'

	input:
		tuple path('lr.asm.hap.bam'), path('lr.asm.hap.bam.bai')
		tuple path('asm.fasta'), path('asm.fasta.fai'), path('asm.bed')
		tuple path('lr.asm.vcf.gz'), path('lr.asm.vcf.gz.tbi')
		path asm_cov

	output:
		tuple path('assembly.no_haplotig.fasta'), path('assembly.no_haplotig.fasta.fai'), path('assembly.no_haplotig.bed')

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	bcftools=\${BCFTOOLS:-${params.tools.bcftools.bin}}
	seqtk=\${SEQTK:-${params.tools.seqtk.bin}}

	get_phase_cov_py=\${GET_PHASE_COV_PY:-${params.python.script.get_phase_cov}}

	MEAN_COV=\$(tail -n 1 ${asm_cov} | cut -f3)
	head -n -1 ${asm_cov} | awk -v meancov=\$MEAN_COV '{if (\$3<0.6*meancov) {print \$1}}' > lr.asm.low_cov.list
	\${bcftools} query -i \"(FILTER=='PASS') && ((GT=='AR') & (GQ>=${params.pipeline.variant_filter.min_gq}) & (AD>=${params.pipeline.variant_filter.min_ad}) & \
	(VAF>=${params.pipeline.variant_filter.min_vaf}) & (VAF<=${params.pipeline.variant_filter.max_vaf}))\" -f '%CHROM\\n' lr.asm.vcf.gz | sort | uniq -c > pepper.var.count

	while read CONTIG_NAME
	do
		VARIANT_COUNT=\$(awk -v contig=\"\${CONTIG_NAME}\" '{if (\$2==contig) {print \$1}}' pepper.var.count)
		if [ -z \${VARIANT_COUNT} ] || [ \${VARIANT_COUNT} -lt 2 ]; then echo -e \"\${CONTIG_NAME}\"; fi

	done < lr.asm.low_cov.list > haplotig.1.list

	while read CONTIG_NAME
	do
		awk -v contig=\"\${CONTIG_NAME}\" '{if (\$1==contig) {print \$0}}' asm.bed;

	done < lr.asm.low_cov.list > lr.asm.low_cov.bed

	python \${get_phase_cov_py} -t ${task.cpus} -b lr.asm.hap.bam -m ${params.pipeline.min_mapq.strict} -r lr.asm.low_cov.bed > MARGIN_PHASED.haplotagged.cov.list
	awk '{COV=\$4+\$5+\$6; COV_THRESHOLD=0.8*COV; if ((\$5>=COV_THRESHOLD) || (\$6>=COV_THRESHOLD)) {print \$1}}' MARGIN_PHASED.haplotagged.cov.list > haplotig.2.list
	cat haplotig.1.list haplotig.2.list | sort | uniq > haplotigs.list; comm -23 <(awk '{print \$1}' asm.bed | sort -k1,1) haplotigs.list > contig.collapsed.list
	\${seqtk} subseq asm.fasta contig.collapsed.list > assembly.no_haplotig.fasta
	\${samtools} faidx assembly.no_haplotig.fasta
	awk '{print \$1 \"\\t0\\t\" \$2}' assembly.no_haplotig.fasta.fai | sort -k1,1 -k2,2n > assembly.no_haplotig.bed
	"""
}

process prepLocalCorrection  {

	// This process makes the assumption that less than 10 billion jobs will be used :)
	// Hence, all files IDs are encoded using 9 digits

	label 'medium_node'

	input:
		tuple path('lr.bam'), path('lr.bam.bai')
		tuple path('sr.bam'), path('sr.bam.bai')

	output:
		path("split/lr.corrected.?????????.*") // Should include all the lr.corrected.?????????.hc.bed, split/lr.corrected.?????????.lmc.bed, split/lr.corrected.?????????.bed, split/lr.corrected.?????????.hc.dict


	shell '/bin/bash', '-euo', 'pipefail'

	"""
	bedtools=\${BEDTOOLS:-${params.tools.bedtools.bin}}

	local_ratatosk_sh=\${RATATOSK_LOCAL_CORRECTION:-${params.scripts.local_ratatosk}}

	mkdir -p split
	bash \${local_ratatosk_sh} -t ${task.cpus} -c -l lr.bam -s sr.bam -o lr.corrected
	cat <(sed 's/\$/\\t1/' lr.corrected.hc.bed) <(sed 's/\$/\\t2/' lr.corrected.lmc.bed) | sort -k1,1 -k2,2n | \
	awk 'BEGIN {SUM=0; ID=0; FS=\"\\t\"; OFS=\"\\t\"} {print \$1,\$2,\$3,\$4,ID,\"lr.corrected.\" ID \".bed\"; SUM+=\$3-\$2; \
	if (SUM>=${params.pipeline.correction.job_len_region}) {ID+=1; SUM=0;}}' > lr.corrected.split.bed
	cd split
	awk '{print>>\$6}' ../lr.corrected.split.bed
	ID_MAX=\$(tail -n 1 ../lr.corrected.split.bed | awk '{print \$5}')
	rm -rf lr.corrected.split.bed

	for i in \$(seq 0 \${ID_MAX})
	do
		ID_JOB=\$(printf \"%09d\" \$i)
		awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {if (\$4==1){print \$1,\$2,\$3}}' lr.corrected.\${i}.bed > lr.corrected.\${ID_JOB}.hc.bed
		awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {if (\$4==2){print \$1,\$2,\$3}}' lr.corrected.\${i}.bed > lr.corrected.\${ID_JOB}.lmc.bed
		rm -rf lr.corrected.\${i}.bed
		\${bedtools} merge -i <(cat lr.corrected.\${ID_JOB}.hc.bed lr.corrected.\${ID_JOB}.lmc.bed | sort -k1,1 -k2,2n) > lr.corrected.\${ID_JOB}.bed
		join ../lr.corrected.dict <(cut -f1 lr.corrected.\${ID_JOB}.bed | sort | uniq) | awk '{print \$1 \"\\t\" \$2}' > lr.corrected.\${ID_JOB}.dict
	done
	"""
}

process localCorrection  {

	label 'medium_node'

	input:
		//tuple val(id), path('in.lr.corr.bed'), path('in.lr.corr.dict'), path('in.lr.corr.hc.bed'), path('in.lr.corr.lmc.bed')
		tuple val(id), path("*")
		tuple path('lr.bam'), path('lr.bam.bai')
		tuple path('sr.bam'), path('sr.bam.bai')

	output:
		path "${id}.lr.corr.fastq.gz"

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	pigz=\${PIGZ_BIN:-${params.tools.pigz.bin}}

	local_ratatosk_sh=\${RATATOSK_LOCAL_CORRECTION:-${params.scripts.local_ratatosk}}

	if [ -s lr.corrected.${id}.bed ]
	then

		bash \${local_ratatosk_sh} -t ${task.cpus} -u -Q ${params.max_lr_bq} -l lr.bam -s sr.bam -C lr.corrected.${id} -o ${id}.lr.corr

		if [ ! \$? -eq 0 ]; then echo \"Something went wrong during the local correction of a long read chunk\" 1>&2; exit 1; fi

		\${pigz} -p ${task.cpus} ${id}.lr.corr.fastq
		\${pigz} -p ${task.cpus} -t ${id}.lr.corr.fastq.gz

		if [ ! \$? -eq 0 ]; then echo \"File ${id}.lr.corr.fastq.gz is malformed\" 1>&2; exit 1; fi
	fi
	"""
}

process mergeLocalCorrection {

	label 'medium_node'

	input:
		path lr_chunks_fq_gz
		tuple path('lr.asm.bam'), path('lr.asm.bam.bai')

	output:
		path 'lr.corrected.fastq.gz'

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	pigz=\${PIGZ_BIN:-${params.tools.pigz.bin}}
	seqkit=\${SEQKIT:-${params.tools.seqkit.bin}}

	cat *.fastq.gz > lr.corrected.fastq.gz
	\${seqkit} rmdup -j ${task.cpus} -n -o lr.corrected.nodup.fastq lr.corrected.fastq.gz
	rm -rf lr.corrected.fastq.gz
	\${pigz} -p ${task.cpus} lr.corrected.nodup.fastq
	mv -f lr.corrected.nodup.fastq.gz lr.corrected.fastq.gz
	\${samtools} view -@ ${task.cpus} -b -f 4 lr.asm.bam | \${samtools} bam2fq -n -@ ${task.cpus} - 2> /dev/null | \${pigz} -p ${task.cpus} -c >> lr.corrected.fastq.gz
	\${pigz} -p ${task.cpus} -t lr.corrected.fastq.gz
	if [ ! \$? -eq 0 ]; then echo \"File lr.corrected.fastq.gz is malformed\" 1>&2; exit 1; fi
	"""
}

process hapResAsm {

	label 'large_node'

	input:
		path filtered_lr_fq_gz

	output:
		tuple path('assembly.fasta'), path('assembly.fasta.fai'), path('assembly.bed')

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	flye=\${FLYE:-${params.tools.flye.bin}}
	seqkit=\${SEQKIT:-${params.tools.seqkit.bin}}

	FLYE_SUBSAMPLE_PARAM=\"\"
	if [ ${params.genome_size} -gt 0 ]; then FLYE_SUBSAMPLE_PARAM=\"--asm-coverage ${params.pipeline.assembly.subsampling_cov} --genome-size ${params.genome_size}\"; fi

	\${flye} --nano-hq ${filtered_lr_fq_gz} -o . -t ${task.cpus} --keep-haplotypes --read-error ${params.pipeline.assembly.max_error_rate.strict} --scaffold \${FLYE_SUBSAMPLE_PARAM}
	if [ ! \$? -eq 0 ]; then echo \"Flye has exited prematurely\" 1>&2; exit 1; fi
	ASM_COMPLETED=\$( { grep \"INFO: Final assembly\" <(tail -n 50 flye.log) || true; } | wc -l)
	if [ \${ASM_COMPLETED} -eq 0 ] || [ ! -s assembly.fasta ]; then echo \"Flye log incomplete or empty assembly\" 1>&2; exit 1; fi

	\${seqkit} seq -w0 assembly.fasta > assembly.tmp.fasta;
	mv -f assembly.tmp.fasta assembly.fasta;

	\${samtools} faidx assembly.fasta
	awk '{print \$1 \"\\t0\\t\" \$2}' assembly.fasta.fai | sort -k1,1 -k2,2n > assembly.bed

	rm -rf 00-assembly 10-consensus 20-repeat 30-contigger 40-polishing
	"""
}

process hapResAsm_filter {

	label 'medium_node'

	input:
		path lr_fq_gz
		tuple path(asm_fa), path(asm_fai), path(asm_bed)

	output:
		tuple path('assembly.filtered.fasta'), path('assembly.filtered.fasta.fai'), path('assembly.filtered.bed')

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	minimap2=\${MINIMAP2:-${params.tools.minimap2.bin}}
	seqtk=\${SEQTK:-${params.tools.seqtk.bin}}

	TASK_MEM=\$(echo -e \"${task.memory}\" | cut -d \" \" -f1)
	MEM_PER_THREADS_SORT=\$(bc -l <<< \"(\${TASK_MEM} / ${task.cpus}) * ${params.tools.samtools.sort.mem_safety_ratio} * 1000\" | awk '{printf(\"%.0f\", \$0)}')

	\${minimap2} -t ${task.cpus} ${params.tools.minimap2.param.ont_reads} -Y ${asm_fa} ${lr_fq_gz} > lr.asm.sam
	\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M lr.asm.sam > lr.asm.bam
	if [ ! -s lr.asm.bam ]; then echo \"File lr.asm.bam does not exist or is empty\" 1>&2; exit 1; fi
	\${samtools} quickcheck lr.asm.bam
	if [ ! \$? -eq 0 ]; then echo \"File lr.asm.bam is malformed\" 1>&2; exit 1; fi
	\${samtools} index -@ ${task.cpus} lr.asm.bam; rm -rf lr.asm.sam

	\${samtools} depth -@ ${task.cpus}-J -Q ${params.pipeline.min_mapq_strict} -G ${params.pipeline.samtools.depth.filter_sec} lr.asm.bam | \
	awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"; CONTIG=\"\"; POS_S=-1; POS_E=-1} {if (\$3>=${params.pipeline.variant_filter.min_depth}) {if ((\$1==CONTIG) && (\$2==POS_E)) {POS_E+=1} else \
	{ if ((POS_S!=-1) && (POS_E-POS_S>=${params.pipeline.min_len_contig})) {print CONTIG, (POS_S-1), (POS_E-1)}; CONTIG=\$1; POS_S=\$2; POS_E=\$2+1;}}} \
	END {if ((POS_S!=-1) && (POS_E-POS_S>=${params.pipeline.min_len_contig})) {print CONTIG, (POS_S-1), (POS_E-1)}}' > lr.asm.filtered.bed
	\${seqtk} subseq assembly.fasta lr.asm.filtered.bed | tr ':' '_' > assembly.filtered.fasta
	\${samtools} faidx assembly.filtered.fasta
	awk '{print \$1 \"\\t0\\t\" \$2}' assembly.filtered.fasta.fai | sort -k1,1 -k2,2n > assembly.filtered.bed
	"""
}

process hapResAsm_polish1_1 {

	label 'medium_node'

	input:
		path lr_fq_gz
		tuple path(asm_fa), path(asm_fai), path(asm_bed)
		tuple path('lr.bam'), path('lr.bam.bai')
		path lr_cov

	output:
		tuple path('lr.asm.minCovLen.filtered.2.bam'), path('lr.asm.minCovLen.filtered.2.bam.bai'), emit: bam
		path 'lr.asm.minCovLen.filtered.2.collapsed.bed', emit: bed

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	minimap2=\${MINIMAP2:-${params.tools.minimap2.bin}}

	filter_align_py=\${FILTER_ALIGN_PY:-${params.python.script.filter_align}}

	TASK_MEM=\$(echo -e \"${task.memory}\" | cut -d \" \" -f1)
	MEM_PER_THREADS_SORT=\$(bc -l <<< \"(\${TASK_MEM} / ${task.cpus}) * ${params.tools.samtools.sort.mem_safety_ratio} * 1000\" | awk '{printf(\"%.0f\", \$0)}')

	\${minimap2} -t ${task.cpus} -ax map-hifi -Y ${asm_fa} ${lr_fq_gz} > lr.asm.minCovLen.sam
	\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M lr.asm.minCovLen.sam > lr.asm.minCovLen.bam
	if [ ! -s lr.asm.minCovLen.bam ]; then echo \"File lr.asm.minCovLen.bam does not exist or is empty\" 1>&2; exit 1; fi
	\${samtools} quickcheck lr.asm.minCovLen.bam
	if [ ! \$? -eq 0 ]; then echo \"File lr.asm.minCovLen.bam is malformed\" 1>&2; exit 1; fi
	\${samtools} index -@ ${task.cpus} lr.asm.minCovLen.bam
	rm -rf lr.asm.minCovLen.sam

	python \${filter_align_py} -t ${task.cpus} -i lr.asm.minCovLen.bam -o lr.asm.minCovLen.tmp.bam -l ${params.pipeline.min_len_alignment} \
	-e ${params.pipeline.assembly.max_error_rate.lenient} -m ${params.pipeline.min_mapq.strict}
	mv -f lr.asm.minCovLen.tmp.bam lr.asm.minCovLen.ba
	\${samtools} index -@ ${task.cpus} lr.asm.minCovLen.bam

	\${samtools} depth -@ ${task.cpus} -aa -J -G ${params.pipeline.samtools.depth.filter_supp_sec} lr.asm.minCovLen.bam | \
	awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"; CONTIG=\"\"; SUM=0; COUNT=0;} \
	{if (\$1!=CONTIG) {if (CONTIG!=\"\") {MEAN=0; if (COUNT>=1) {MEAN=SUM/COUNT}; print CONTIG, MEAN}; CONTIG=\$1; SUM=0; COUNT=0}; SUM+=\$3; COUNT+=1} \
	END {if (CONTIG!=\"\") {MEAN=0; if (COUNT>=1) {MEAN=SUM/COUNT}; print CONTIG, MEAN}}' | sort -k1,1 > lr.asm.minCovLen.prim.cov;
	\${samtools} depth -@ ${task.cpus} -aa -J -G ${params.pipeline.samtools.depth.filter_sec} lr.asm.minCovLen.bam | \
	awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"; CONTIG=\"\"; SUM=0; COUNT=0;} \
	{if (\$1!=CONTIG) {if (CONTIG!=\"\") {MEAN=0; if (COUNT>=1) {MEAN=SUM/COUNT}; print CONTIG, MEAN}; CONTIG=\$1; SUM=0; COUNT=0}; SUM+=\$3; COUNT+=1} \
	END {if (CONTIG!=\"\") {MEAN=0; if (COUNT>=1) {MEAN=SUM/COUNT}; print CONTIG, MEAN}}' | sort -k1,1 > lr.asm.minCovLen.primsup.cov;
	join lr.asm.minCovLen.prim.cov lr.asm.minCovLen.primsup.cov | awk '{if ((\$3-\$2)<(2*\$2)){print \$1}}' > lr.asm.keep.1.list;
	join lr.asm.keep.1.list <(sort -k1,1 ${asm_bed}) | sort -k1,1 -k2,2n > lr.asm.keep.1.bed;
	\${samtools} view -@ ${task.cpus} -b -M -L lr.asm.keep.1.bed lr.asm.minCovLen.bam > lr.asm.minCovLen.filtered.1.bam;
	\${samtools} index -@ ${task.cpus} lr.asm.minCovLen.filtered.1.bam;

	COV=\$(tail -n 1 ${lr_cov} | awk '{printf(\"%.0f\\n\", \$3)}');
	COLLAPSED_COV_MIN=\$(echo -e \"\${COV}\" | awk '{printf(\"%.0f\\n\", ${params.pipeline.ratio_cov.lower_bound.strict}*\$0)}');
	COLLAPSED_COV_MAX=\$(echo -e \"\${COV}\" | awk '{printf(\"%.0f\\n\", ${params.pipeline.ratio_cov.upper_bound.strict}*\$0)}');
	SUPER_COV_MIN=\$(echo -e \"\${COV}\" | awk '{printf(\"%.0f\\n\", ${params.pipeline.ratio_cov.lower_bound.lenient}*(\$0/2))}');
	SUPER_COV_MAX=\$(echo -e \"\${COV}\" | awk '{printf(\"%.0f\\n\", ${params.pipeline.ratio_cov.upper_bound.lenient}*\$0)}');

	\${samtools} depth -@ ${task.cpus} -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} lr.asm.minCovLen.filtered.1.bam | \
	awk -v sup_min_cov=\${SUPER_COV_MIN} -v sup_max_cov=\${SUPER_COV_MAX} 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"; CONTIG=\"\"; POS_S=-1; POS_E=-1} \
	{if ((\$3>=sup_min_cov) && (\$3<=sup_max_cov)) {if ((\$1==CONTIG) && (\$2==POS_E)) {POS_E+=1} else { if (POS_S!=-1) {print CONTIG, (POS_S-1), (POS_E-1)}; CONTIG=\$1; POS_S=\$2; POS_E=\$2+1;}}} \
	END {if (POS_S!=-1) {print CONTIG, (POS_S-1), (POS_E-1)}}' > lr.asm.keep.2.bed;
	\${samtools} view -@ ${task.cpus} -b -M -L lr.asm.keep.2.bed lr.asm.minCovLen.filtered.1.bam > lr.asm.minCovLen.filtered.2.bam;
	\${samtools} index -@ ${task.cpus} lr.asm.minCovLen.filtered.2.bam;
	\${samtools} depth -@ ${task.cpus} -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} lr.asm.minCovLen.filtered.2.bam | \
	awk -v min_cov=\${COLLAPSED_COV_MIN} -v max_cov=\${COLLAPSED_COV_MAX} 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"; CONTIG=\"\"; POS_S=-1; POS_E=-1} \
	{if ((\$3>=min_cov) && (\$3<=max_cov)) {if ((\$1==CONTIG) && (\$2==POS_E)) {POS_E+=1} else { if (POS_S!=-1) {print CONTIG, (POS_S-1) , (POS_E-1)}; CONTIG=\$1; POS_S=\$2; POS_E=\$2+1;}}} \
	END {if (POS_S!=-1) {print CONTIG, (POS_S-1), (POS_E-1)}}' > lr.asm.minCovLen.filtered.2.collapsed.bed;
	"""
}

process hapResAsm_polish1_2 {

	label 'medium_node'

	input:
		tuple path(asm_fa), path(asm_fai), path(asm_bed)
		tuple path('lr.asm.minCovLen.filtered.2.bam'), path('lr.asm.minCovLen.filtered.2.bam.bai')
		path 'MARGIN_PHASED.phaseset.bed'
		path lr_cov

	output:
		tuple path('assembly.filtered.splitPS.fasta'), path('assembly.filtered.splitPS.fasta.fai'), path('assembly.filtered.splitPS.bed'), emit: asm
		tuple path('lr.asm.minCovLen.splitPS.bam'), path('lr.asm.minCovLen.splitPS.bam.bai'), emit: bam
		path 'lr.asm.minCovLen.splitPS.collapsed.bed', emit: bed

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	minimap2=\${MINIMAP2:-${params.tools.minimap2.bin}}
	seqtk=\${SEQTK:-${params.tools.seqtk.bin}}
	pigz=\${PIGZ_BIN:-${params.tools.pigz.bin}}

	filter_align_py=\${FILTER_ALIGN_PY:-${params.python.script.filter_align}}

	TASK_MEM=\$(echo -e \"${task.memory}\" | cut -d \" \" -f1)
	MEM_PER_THREADS_SORT=\$(bc -l <<< \"(\${TASK_MEM} / ${task.cpus}) * ${params.tools.samtools.sort.mem_safety_ratio} * 1000\" | awk '{printf(\"%.0f\", \$0)}')

	COV=\$(tail -n 1 ${lr_cov} | awk '{printf(\"%.0f\\n\", \$3)}')
	COLLAPSED_COV_MIN=\$(echo -e \"\${COV}\" | awk '{printf(\"%.0f\\n\", ${params.pipeline.ratio_cov.lower_bound.strict}*\$0)}')
	COLLAPSED_COV_MAX=\$(echo -e \"\${COV}\" | awk '{printf(\"%.0f\\n\", ${params.pipeline.ratio_cov.upper_bound.strict}*\$0)}')

	sort -k1,1 -k2,2n MARGIN_PHASED.phaseset.bed > MARGIN_PHASED.phaseset.tmp.bed
	mv -f MARGIN_PHASED.phaseset.tmp.bed MARGIN_PHASED.phaseset.bed
	awk '{if (\$4!=\"ContigEnd\"){print \$1}}' MARGIN_PHASED.phaseset.bed | cut -f1 | sort | uniq > MARGIN_PHASED.phaseset.inconsistent.list
	join MARGIN_PHASED.phaseset.bed MARGIN_PHASED.phaseset.inconsistent.list | cut -d ' ' -f1,2,3 | tr ' ' '\\t' > MARGIN_PHASED.phaseset.inconsistent.bed
	comm -13 MARGIN_PHASED.phaseset.inconsistent.list <(awk '{print \$1}' ${asm_bed} | sort -k1,1) > MARGIN_PHASED.phaseset.consistent.list
	join MARGIN_PHASED.phaseset.inconsistent.bed <(sort -k1,1 ${asm_bed}) | cut -d ' ' -f1,2,3,4,5 | tr ' ' '\\t' | \
	awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"; CHROM=\"\"; SPOS=0; EPOS=0; EPOS_CONTIG=0} \
	{if (\$1!=CHROM){if (CHROM!=\"\"){print CHROM, SPOS, EPOS_CONTIG}; CHROM=\$1; SPOS=0; EPOS=\$3; EPOS_CONTIG=\$5} else \
	{EPOS+=int((\$2-EPOS)/2); print CHROM, SPOS, EPOS; CHROM=\$1; SPOS=EPOS; EPOS=\$3; EPOS_CONTIG=\$5}} \
	END {if (CHROM!=\"\"){print CHROM, SPOS, EPOS_CONTIG}}' > MARGIN_PHASED.phaseset.inconsistent.split.bed
	cat <(\${seqtk} subseq ${asm_fa} MARGIN_PHASED.phaseset.consistent.list) <(\${seqtk} subseq ${asm_fa} MARGIN_PHASED.phaseset.inconsistent.split.bed) | \
	awk '{if (NR%2==1){print \">contig_\" ID; ID+=1} else {print \$0}}' > assembly.filtered.splitPS.fasta
	\${samtools} faidx assembly.filtered.splitPS.fasta
	awk '{print \$1 \"\\t0\\t\" \$2}' assembly.filtered.splitPS.fasta.fai | sort -k1,1 -k2,2n > assembly.filtered.splitPS.bed

	\${samtools} bam2fq -@ ${task.cpus} -n lr.asm.minCovLen.filtered.2.bam | \${pigz} -p ${task.cpus} -c > lr.asm.minCovLen.filtered.2.fq.gz
	\${minimap2} -t ${task.cpus} -ax map-hifi -Y assembly.filtered.splitPS.fasta lr.asm.minCovLen.filtered.2.fq.gz > lr.asm.minCovLen.splitPS.sam
	\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M lr.asm.minCovLen.splitPS.sam > lr.asm.minCovLen.splitPS.bam
	if [ ! -s lr.asm.minCovLen.splitPS.bam ]; then echo \"File lr.asm.minCovLen.splitPS.bam does not exist or is empty\" 1>&2; exit 1; fi
	\${samtools} quickcheck lr.asm.minCovLen.splitPS.bam
	if [ ! \$? -eq 0 ]; then echo \"File lr.asm.minCovLen.splitPS.bam is malformed\" 1>&2; exit 1; fi
	\${samtools} index -@ ${task.cpus} lr.asm.minCovLen.splitPS.bam
	rm -rf lr.asm.minCovLen.splitPS.sam lr.asm.minCovLen.filtered.2.fq.gz

	python \${filter_align_py} -t ${task.cpus} -i lr.asm.minCovLen.splitPS.bam -o lr.asm.minCovLen.splitPS.tmp.bam -l ${params.pipeline.min_len_alignment} \
	-e ${params.pipeline.assembly.max_error_rate.lenient} -m ${params.pipeline.min_mapq.strict}
	mv -f lr.asm.minCovLen.splitPS.tmp.bam lr.asm.minCovLen.splitPS.bam
	\${samtools} index -@ ${task.cpus} lr.asm.minCovLen.splitPS.bam
	\${samtools} depth -@ ${task.cpus} -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} lr.asm.minCovLen.splitPS.bam | \
	awk -v min_cov=\${COLLAPSED_COV_MIN} -v max_cov=\${COLLAPSED_COV_MAX} 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"; CONTIG=\"\"; POS_S=-1; POS_E=-1} \
	{if ((\$3>=min_cov) && (\$3<=max_cov)) {if ((\$1==CONTIG) && (\$2==POS_E)) {POS_E+=1} else { if (POS_S!=-1) {print CONTIG, (POS_S-1), (POS_E-1)}; CONTIG=\$1; POS_S=\$2; POS_E=\$2+1;}}} \
	END {if (POS_S!=-1) {print CONTIG, (POS_S-1), (POS_E-1)}}' > lr.asm.minCovLen.splitPS.collapsed.bed
	"""
}

process hapResAsm_polish2 {

	label 'medium_node'

	input:
		tuple path(asm_fa), path(asm_fai), path(asm_bed)
		tuple path('hap_lr_splitPS.bam'), path('hap_lr_splitPS.bam.bai')
		tuple path('hap_lr_splitPS.vcf.gz'), path('hap_lr_splitPS.vcf.gz.tbi')

	output:
		tuple path('detectCorruptPS/ref-to-ref-map/polished_1.filtered.fasta'), path('detectCorruptPS/ref-to-ref-map/polished_1.filtered.fasta.fai'), path('detectCorruptPS/ref-to-ref-map/polished_1.filtered.bed')
		

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	minimap2=\${MINIMAP2:-${params.tools.minimap2.bin}}
	flye=\${FLYE:-${params.tools.flye.bin}}
	seqkit=\${SEQKIT:-${params.tools.seqkit.bin}}
	seqtk=\${SEQTK:-${params.tools.seqtk.bin}}
	bedtools=\${BEDTOOLS:-${params.tools.bedtools.bin}}
	pigz=\${PIGZ_BIN:-${params.tools.pigz.bin}}

	extract_ref_reads_py=\${EXTRACT_REF_READS_PY:-${params.python.script.extract_ref_reads}}
	get_ps_info_py=\${GET_PS_INFO_PY:-${params.python.script.get_ps_info}}
	extract_reads_within_py=\${EXTRACT_READS_WITHIN_PY:-${params.python.script.extract_reads_within}}

	TASK_MEM=\$(echo -e \"${task.memory}\" | cut -d \" \" -f1)
	MEM_PER_THREADS_SORT=\$(bc -l <<< \"(\${TASK_MEM} / ${task.cpus}) * ${params.tools.samtools.sort.mem_safety_ratio} * 1000\" | awk '{printf(\"%.0f\", \$0)}')

	python \${extract_ref_reads_py} -t ${task.cpus} -m ${params.pipeline.min_mapq.strict} -b hap_lr_splitPS.bam \
	-v hap_lr_splitPS.vcf.gz -B MARGIN_PHASED.haplotagged.polish.bam --include_supplementary
	\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M MARGIN_PHASED.haplotagged.polish.bam > MARGIN_PHASED.haplotagged.polish.sorted.bam
	mv -f MARGIN_PHASED.haplotagged.polish.sorted.bam MARGIN_PHASED.haplotagged.polish.bam
	if [ ! -s MARGIN_PHASED.haplotagged.polish.bam ]; then echo \"File MARGIN_PHASED.haplotagged.polish.bam does not exist or is empty\" 1>&2; exit 1; fi
	\${samtools} quickcheck MARGIN_PHASED.haplotagged.polish.bam
	if [ ! \$? -eq 0 ]; then echo \"File MARGIN_PHASED.haplotagged.polish.bam is malformed\" 1>&2; exit 1; fi
	\${samtools} index -@ ${task.cpus} MARGIN_PHASED.haplotagged.polish.bam
	
	mkdir -p detectCorruptPS
	cd detectCorruptPS
	python \${get_ps_info_py} -v ../hap_lr_splitPS.vcf.gz | \
	awk '{SUM=\$5+\$6; if (((\$5>\$6) && (\$5/SUM<0.75)) || ((\$6>=\$5) && (\$6/SUM<0.75))) {print \$1 \"\\t\" \$2 \"\\t\" \$3}}' > MARGIN_PHASED.corrupt_ps.bed
	\${bedtools} subtract -a ../${asm_bed} -b MARGIN_PHASED.corrupt_ps.bed | sort -k1,1 -k2,2n > MARGIN_PHASED.no_corrupt_ps.bed

	if [ -s MARGIN_PHASED.corrupt_ps.bed ]
	then
		\${samtools} view -@ ${task.cpus} -b -M -L MARGIN_PHASED.corrupt_ps.bed ../MARGIN_PHASED.haplotagged.polish.bam > MARGIN_PHASED.haplotagged.polish.corruptPS.bam
		\${samtools} index -@ ${task.cpus} MARGIN_PHASED.haplotagged.polish.corruptPS.bam
		python \${extract_reads_within_py} -t ${task.cpus} -b ../MARGIN_PHASED.haplotagged.polish.bam -r MARGIN_PHASED.no_corrupt_ps.bed -B MARGIN_PHASED.haplotagged.polish.noCorruptPS.bam -s -S
		\${samtools} index -@ ${task.cpus} MARGIN_PHASED.haplotagged.polish.noCorruptPS.bam
		\${flye} --polish-target ../${asm_fa} --nano-hq MARGIN_PHASED.haplotagged.polish.corruptPS.bam -t ${task.cpus} -o flye.polish.corruptPS
		\${seqkit} seq -w0 -m${params.pipeline.min_len_collapsed_segment} flye.polish.corruptPS/polished_1.fasta > flye.polish.corruptPS/polished_1.tmp.fasta
		mv -f flye.polish.corruptPS/polished_1.tmp.fasta flye.polish.corruptPS/polished_1.fasta
		\${flye} --polish-target ../${asm_fa} --nano-hq MARGIN_PHASED.haplotagged.polish.noCorruptPS.bam -t ${task.cpus} -o flye.polish.noCorruptPS
		\${seqkit} seq -w0 -m${params.pipeline.min_len_collapsed_segment} flye.polish.noCorruptPS/polished_1.fasta > flye.polish.noCorruptPS/polished_1.tmp.fasta
		mv -f flye.polish.noCorruptPS/polished_1.tmp.fasta flye.polish.noCorruptPS/polished_1.fasta
		cat flye.polish.corruptPS/polished_1.fasta flye.polish.noCorruptPS/polished_1.fasta | awk 'BEGIN {ID=0} {if (substr(\$0,1,1)==\">\") {print \">contig_\" ID; ID+=1} else {print \$0}}' > polished_1.fasta
	else
		\${flye} --polish-target ../${asm_fa} --nano-hq ../MARGIN_PHASED.haplotagged.polish.bam -t ${task.cpus} -o flye.polish.noCorruptPS
		\${seqkit} seq -w0 -m${params.pipeline.min_len_collapsed_segment} flye.polish.noCorruptPS/polished_1.fasta > flye.polish.noCorruptPS/polished_1.tmp.fasta
		mv -f flye.polish.noCorruptPS/polished_1.tmp.fasta flye.polish.noCorruptPS/polished_1.fasta
		awk 'BEGIN {ID=0} {if (substr(\$0,1,1)==\">\") {print \">contig_\" ID; ID+=1} else {print \$0}}' flye.polish.noCorruptPS/polished_1.fasta > polished_1.fasta
	fi

	\${samtools} faidx polished_1.fasta
	awk '{print \$1 \"\\t0\\t\" \$2}' polished_1.fasta.fai | sort -k1,1 -k2,2n > polished_1.bed
	\${samtools} bam2fq -@ ${task.cpus} -n ../MARGIN_PHASED.haplotagged.polish.bam | \${pigz} -p ${task.cpus} -c > MARGIN_PHASED.haplotagged.polish.fastq.gz
	\${minimap2} -t ${task.cpus} -ax map-hifi -Y polished_1.fasta MARGIN_PHASED.haplotagged.polish.fastq.gz > lr.asm.sam
	\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M lr.asm.sam > lr.asm.bam
	if [ ! -s lr.asm.bam ]; then echo \"File lr.asm.bam does not exist or is empty\" 1>&2; exit 1; fi
	\${samtools} quickcheck lr.asm.bam
	if [ ! \$? -eq 0 ]; then echo \"File lr.asm.bam is malformed\" 1>&2; exit 1; fi
	\${samtools} index -@ ${task.cpus} lr.asm.bam; rm -rf lr.asm.sam MARGIN_PHASED.haplotagged.polish.fastq.gz

	\${samtools} depth -@ ${task.cpus} -aa -J -Q 0 -G ${params.pipeline.samtools.depth.filter_sec} lr.asm.bam | \
	awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"; CONTIG=\"\"; POS_S=-1; POS_E=-1} \
	{if (\$3<${params.pipeline.variant_filter.min_depth}) {if ((\$1==CONTIG) && (\$2==POS_E)) {POS_E+=1} else \
	{ if ((POS_S!=-1) && (POS_E-POS_S>=${params.pipeline.min_len_collapsed_segment})) {print CONTIG, (POS_S-1), (POS_E-1)}; CONTIG=\$1; POS_S=\$2; POS_E=\$2+1;}}} \
	END {if ((POS_S!=-1) && (POS_E-POS_S>=${params.pipeline.min_len_collapsed_segment})) {print CONTIG, (POS_S-1), (POS_E-1)}}' > lr.asm.discard.bed
	\${bedtools} subtract -a polished_1.bed -b lr.asm.discard.bed | sort -k1,1 -k2,2n > lr.asm.keep.bed
	\${seqtk} subseq polished_1.fasta lr.asm.keep.bed | tr ':' '_' > polished_1.filtered.fasta
	
	mkdir ref-to-ref-map; cd ref-to-ref-map

	\${minimap2} -t ${task.cpus} ${params.tools.minimap2.param.asm} ../polished_1.filtered.fasta ../polished_1.filtered.fasta > ref2ref.sam
	\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M ref2ref.sam > ref2ref.bam
	if [ ! -s ref2ref.bam ]; then echo \"File ref2ref.bam does not exist or is empty\" 1>&2; exit 1; fi
	\${samtools} quickcheck ref2ref.bam
	if [ ! \$? -eq 0 ]; then echo \"File ref2ref.bam is malformed\" 1>&2; exit 1; fi
	\${samtools} index -@ ${task.cpus} ref2ref.bam; rm -rf ref2ref.sam

	\${samtools} view -@ ${task.cpus} -F 2048 --tag de ref2ref.bam | awk '{if ((\$5==0) && (\$1!=\$3)) {for (i=6; i<=NF; i+=1){ if (substr(\$i,1,5)==\"de:f:\") \
	{ERR=substr(\$i,6,length(\$i)-5); if (ERR<=0.001) {print \$1}; break;}}}}' | sort | uniq > duplicates.list
	awk '{if (substr(\$0,1,1)==\">\"){print substr(\$0,2,length(\$0)-1)}}' ../polished_1.filtered.fasta | sort > contigs.all.list
	comm -23 contigs.all.list duplicates.list > contigs.keep.list
	\${seqtk} subseq ../polished_1.filtered.fasta contigs.keep.list > polished_1.fasta
	\${seqtk} seq -L 3000 polished_1.fasta > polished_1.filtered.fasta
	\${samtools} faidx polished_1.filtered.fasta
	awk '{print \$1 \"\\t0\\t\" \$2}' polished_1.filtered.fasta.fai | sort -k1,1 -k2,2n > polished_1.filtered.bed
	"""
}

process hapResAsmPolish_phase_mapFilter {

	label 'medium_node'

	input:

		path lr_fq_gz
		tuple path(asm_fa), path(asm_fai), path(asm_bed)

	output:

		tuple path('lr.asm.bam'), path('lr.asm.bam.bai')
		
	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	minimap2=\${MINIMAP2:-${params.tools.minimap2.bin}}

	filter_align_py=\${FILTER_ALIGN_PY:-${params.python.script.filter_align}}

	TASK_MEM=\$(echo -e \"${task.memory}\" | cut -d \" \" -f1)
	MEM_PER_THREADS_SORT=\$(bc -l <<< \"(\${TASK_MEM} / ${task.cpus}) * ${params.tools.samtools.sort.mem_safety_ratio} * 1000\" | awk '{printf(\"%.0f\", \$0)}')

	\${minimap2} -t ${task.cpus} -ax map-hifi -Y ${asm_fa} ${lr_fq_gz} > lr.asm.sam
	\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M lr.asm.sam > lr.asm.bam
	if [ ! -s lr.asm.bam ]; then echo \"File lr.asm.bam does not exist or is empty\" 1>&2; exit 1; fi
	\${samtools} quickcheck lr.asm.bam
	if [ ! \$? -eq 0 ]; then echo \"File lr.asm.bam is malformed\" 1>&2; exit 1; fi
	\${samtools} index -@ ${task.cpus} lr.asm.bam
	rm -rf lr.asm.sam

	python \${filter_align_py} -t ${task.cpus} -i lr.asm.bam -o lr.asm.tmp.bam \
	-l ${params.pipeline.min_len_alignment} -e ${params.pipeline.assembly.max_error_rate.lenient} -m ${params.pipeline.min_mapq.strict}
	mv -f lr.asm.tmp.bam lr.asm.bam
	\${samtools} index -@ ${task.cpus} lr.asm.bam
	"""
}

process hapResAsmPolish_phase_varFilter {

	label 'medium_node'

	input:
		tuple path('lr.asm.bam'), path('lr.asm.bam.bai')
		tuple path('lr.asm.vcf.gz'), path('lr.asm.vcf.gz.tbi')
		path lr_cov

	output:
		tuple path('lr.asm.filtered.vcf.gz'), path('lr.asm.filtered.vcf.gz.tbi')
		

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	bcftools=\${BCFTOOLS:-${params.tools.bcftools.bin}}

	\${bcftools} query -i \"(type=='snp') && (FILTER=='PASS') && ((GT=='AR') & (GQ>=${params.pipeline.variant_filter.min_gq}) & (AD>=${params.pipeline.variant_filter.min_ad}) & \
	(VAF>=${params.pipeline.variant_filter.min_vaf}) & (VAF<=${params.pipeline.variant_filter.max_vaf}))\" -f '[%CHROM\\t%POS0\\t%END\n]' lr.asm.vcf.gz > lr.asm.hq.bed;

	COLLAPSED_COV_MIN=\$(tail -n 1 ${lr_cov} | awk '{printf(\"%.0f\\n\", ${params.pipeline.ratio_cov.lower_bound.strict}*\$3)}');
	COLLAPSED_COV_MAX=\$(tail -n 1 ${lr_cov} | awk '{printf(\"%.0f\\n\", ${params.pipeline.ratio_cov.upper_bound.strict}*\$3)}');
	COV=\$(\${samtools} depth -@ ${task.cpus} -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} -b lr.asm.hq.bed lr.asm.bam | \
	awk -v min_cov=\${COLLAPSED_COV_MIN} -v max_cov=\${COLLAPSED_COV_MAX} '{if ((\$3>=min_cov) && (\$3<=max_cov)) {SUM+=\$3; COUNT+=1}} END {if (COUNT==0){print \"0\"} else {printf(\"%.0f\\n\", SUM/COUNT)}}')
	COLLAPSED_COV_MIN=\$(echo -e \"\${COV}\" | awk '{printf(\"%.0f\\n\", ${params.pipeline.ratio_cov.lower_bound.strict}*\$0)}');
	COLLAPSED_COV_MAX=\$(echo -e \"\${COV}\" | awk '{printf(\"%.0f\\n\", ${params.pipeline.ratio_cov.upper_bound.strict}*\$0)}');

	\${samtools} depth -@ ${task.cpus} -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} lr.asm.bam | \
	awk -v min_cov=\${COLLAPSED_COV_MIN} -v max_cov=\${COLLAPSED_COV_MAX} 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"; CONTIG=\"\"; POS_S=-1; POS_E=-1} \
	{if ((\$3>=min_cov) && (\$3<=max_cov)) {if ((\$1==CONTIG) && (\$2==POS_E)) {POS_E+=1} \
	else { if (POS_S!=-1) {print CONTIG, (POS_S-1), (POS_E-1)}; CONTIG=\$1; POS_S=\$2; POS_E=\$2+1;}}} \
	END {if (POS_S!=-1) {print CONTIG, (POS_S-1), (POS_E-1)}}' > lr.asm.collapsed.bed;

	\${bcftools} view -i \"(type=='snp') && (FILTER=='PASS') && ((GT=='AR') & (GQ>=${params.pipeline.variant_filter.min_gq}) & (AD>=${params.pipeline.variant_filter.min_ad}) & \
	(VAF>=${params.pipeline.variant_filter.min_vaf}) & (VAF<=${params.pipeline.variant_filter.max_vaf}))\" -R lr.asm.collapsed.bed -Oz lr.asm.vcf.gz | \
	\${bcftools} sort -Oz -o lr.asm.filtered.vcf.gz
	tabix -p vcf lr.asm.filtered.vcf.gz
	"""
}

process hapResAsmPolish_phase_varPhasedFilter {

	label 'medium_node'

	input:

		tuple path('lr.asm.unfilt.vcf.gz'), path('lr.asm.unfilt.vcf.gz.tbi')
		tuple path('lr.asm.hap.bam'), path('lr.asm.hap.bam.bai')
		tuple path('lr.asm.hap.filt.vcf.gz'), path('lr.asm.hap.filt.vcf.gz.tbi')

	output:
		tuple path('lr.asm.filt.2.vcf.gz'), path('lr.asm.filt.2.vcf.gz.tbi')
		

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	bcftools=\${BCFTOOLS:-${params.tools.bcftools.bin}}

	\${bcftools} query -i \"(type=='snp') & (PS>=0)\" -f '%CHROM\\t%POS0\\t%END\n' lr.asm.hap.filt.vcf.gz > lr.asm.hap.filt.bed
	\${samtools} depth -@ ${task.cpus} -a -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} -b lr.asm.hap.filt.bed lr.asm.hap.bam | \
	awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {print \$1, \$2, (\$2+1), \$3}' > lr.asm.hap.filt.bed.tmp
	mv -f lr.asm.hap.filt.bed.tmp lr.asm.hap.filt.bed
	COV_STDEV=\$(awk 'BEGIN {i=1} {SUM+=\$4; COUNT+=1; SAMP[i]=\$4; i+=1} END {MEAN=SUM/COUNT; SUMSQ=0; for (i=1; i<=COUNT; i++) {SUMSQ+=(SAMP[i]-MEAN)^2}; print MEAN \"\\t\" sqrt(SUMSQ/COUNT)}' lr.asm.hap.filt.bed)
	COV=\$(echo -e \"\${COV_STDEV}\" | awk '{print \$1}')
	AAD_COV=\$(echo -e \"\${COV_STDEV}\" | awk '{print \$2}')

	awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {if (\$1!=CONTIG){if (CONTIG!=\"\") {MEAN=SUM/COUNT; SUMSQ=0; for (i=1; i<=COUNT; i++) {SUMSQ+=(SAMP[i]-MEAN)^2}; \
	print CONTIG, MEAN, sqrt(SUMSQ/COUNT)}; CONTIG=\$1; SUM=\$4; COUNT=1; i=1} else {SUM+=\$4; COUNT+=1}; SAMP[i]=\$4; i+=1} \
	END {if (CONTIG!=\"\") {MEAN=SUM/COUNT; SUMSQ=0; for (i=1; i<=COUNT; i++) {SUMSQ+=(SAMP[i]-MEAN)^2}; print CONTIG, MEAN, sqrt(SUMSQ/COUNT)};}' lr.asm.hap.filt.bed > lr.asm.hap.filt.cov

	cat lr.asm.hap.filt.cov <(\${samtools} depth -@ ${task.cpus} -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} lr.asm.hap.bam) | \
	awk -v COV_GLOB=\$(echo -e \"\${COV}\" | awk '{printf(\"%.0f\\n\", \$0)}') -v AAD_GLOB=\$(echo -e \"\${AAD_COV}\" | awk '{printf(\"%.0f\\n\", \$0)}') -v NB_PHASED_HAP=\$(wc -l lr.asm.hap.filt.cov | cut -d ' ' -f1) \
	'BEGIN {CONTIG=\"\"; POS_S=-1; POS_E=-1; COV_GLOB_SUP_LOW=(COV_GLOB-AAD_GLOB)*${params.pipeline.ratio_cov.lower_bound.lenient}; COV_GLOB_SUP_UP=(COV_GLOB+AAD_GLOB)*${params.pipeline.ratio_cov.upper_bound.lenient}} \
	{if (NR<=NB_PHASED_HAP){COV_LOC[\$1]=\$2; AAD_LOC[\$1]=\$3} else {COV=COV_GLOB; AAD=AAD_GLOB; if (\$1 in COV_LOC){COV=COV_LOC[\$1]; AAD=AAD_LOC[\$1]}; \
	if ((\$3>=(COV-AAD)*${params.pipeline.ratio_cov.lower_bound.strict}) && (\$3<=(COV+AAD)*${params.pipeline.ratio_cov.upper_bound.strict}) && (\$3>=COV_GLOB_SUP_LOW) && (\$3<=COV_GLOB_SUP_UP)) \
	{if ((\$1==CONTIG) && (\$2==POS_E)) {POS_E+=1} else { if (POS_S!=-1) {print CONTIG \"\\t\" (POS_S-1) \"\\t\" (POS_E-1)}; CONTIG=\$1; POS_S=\$2; POS_E=\$2+1}}}} \
	END {if (POS_S!=-1) {print CONTIG \"\\t\" (POS_S-1) \"\\t\" (POS_E-1)}}' > lr.asm.hap.filt.cov.bed

	\${bcftools} view -i \"(type=='snp') && (FILTER=='PASS') && ((GT=='AR') & (GQ>=${params.pipeline.variant_filter.min_gq}) & (AD>=${params.pipeline.variant_filter.min_ad}) & \
	(VAF>=${params.pipeline.variant_filter.min_vaf}) & (VAF<=${params.pipeline.variant_filter.max_vaf}))\" -R lr.asm.hap.filt.cov.bed -Oz lr.asm.unfilt.vcf.gz | \
	\${bcftools} sort -Oz -o lr.asm.filt.vcf.gz; tabix -p vcf lr.asm.filt.vcf.gz;

	\${bcftools} query -i \"type=='snp'\" -f \"%CHROM\\t%POS0\\t%END\n\" lr.asm.filt.vcf.gz > lr.asm.filt.bed;
	\${samtools} mpileup -B -q ${params.pipeline.min_mapq.strict} -Q 0 -l lr.asm.filt.bed lr.asm.hap.bam | cut -f1,2,6 | \
	awk 'BEGIN {for(n=0;n<256;n++){ord[sprintf(\"%c\",n)]=n}} \
	{l=split(\$3,a,\"\"); SUM=0; for(i=1;i<=l;i++){SUM+=ord[a[i]]-33}; MEAN=0; if (l!=0){MEAN=SUM/l}; if (MEAN>=20){print \$1 \"\\t\" (\$2-1) \"\\t\" \$2}}' > lr.asm.filt.hq.bed

	\${bcftools} view -i \"(type=='snp') && (FILTER=='PASS') && ((GT=='AR') & (GQ>=${params.pipeline.variant_filter.min_gq}) & (AD>=${params.pipeline.variant_filter.min_ad}) & \
	(VAF>=${params.pipeline.variant_filter.min_vaf}) & (VAF<=${params.pipeline.variant_filter.max_vaf}))\" -R lr.asm.filt.hq.bed -Oz lr.asm.filt.vcf.gz | \
	\${bcftools} sort -Oz -o lr.asm.filt.2.vcf.gz
	tabix -p vcf lr.asm.filt.2.vcf.gz
	"""
}

process hapResAsmPolish_phase_selectBestPhasedAlignRec {

	label 'medium_node'

	input:
		tuple path('lr.asm.hap.bam'), path('lr.asm.hap.bam.bai')

	output:
		tuple path('lr_asm_hap.phasedSupp2prim.bam'), path('lr_asm_hap.phasedSupp2prim.bam.bai')

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	phased_supp2prim_py=\${PHASED_SUPP2PRIM_PY:-${params.python.script.phased_supp2prim}}

	python \${phased_supp2prim_py} -t ${task.cpus} -b lr.asm.hap.bam -B lr_asm_hap.phasedSupp2prim.bam;
	\${samtools} index -@ ${task.cpus} lr_asm_hap.phasedSupp2prim.bam;
	"""	
}

process hapResAsmPolish_polishAlt {

	label 'medium_node'

	input:
		tuple path(asm_fa), path(asm_fai), path(asm_bed)
		tuple path('hap_lr_polish.bam'), path('hap_lr_polish.bam.bai')
		tuple path('hap_lr_polish.vcf.gz'), path('hap_lr_polish.vcf.gz.tbi')

	output:
		path 'flye.polish.alt/assembly.all.fasta'

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	flye=\${FLYE:-${params.tools.flye.bin}}
	seqkit=\${SEQKIT:-${params.tools.seqkit.bin}}

	extract_ref_reads_py=\${EXTRACT_REF_READS_PY:-${params.python.script.extract_ref_reads}}

	TASK_MEM=\$(echo -e \"${task.memory}\" | cut -d \" \" -f1);
	MEM_PER_THREADS_SORT=\$(bc -l <<< \"(\${TASK_MEM} / ${task.cpus}) * ${params.tools.samtools.sort.mem_safety_ratio} * 1000\" | awk '{printf(\"%.0f\", \$0)}');

	python \${extract_ref_reads_py} -t ${task.cpus} -b hap_lr_polish.bam -v hap_lr_polish.vcf.gz -B hap_lr_polish.alt.bam -p --alt
	\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M hap_lr_polish.alt.bam > hap_lr_polish.alt.sorted.bam
	if [ ! -s hap_lr_polish.alt.sorted.bam ]; then echo \"File hap_lr_polish.alt.sorted.bam does not exist or is empty\" 1>&2; exit 1; fi
	\${samtools} quickcheck hap_lr_polish.alt.sorted.bam
	if [ ! \$? -eq 0 ]; then echo \"File hap_lr_polish.alt.sorted.bam is malformed\" 1>&2; exit 1; fi
	mv -f hap_lr_polish.alt.sorted.bam hap_lr_polish.alt.bam
	\${samtools} index -@ ${task.cpus} hap_lr_polish.alt.bam
	\${flye} --polish-target ${asm_fa} --nano-hq hap_lr_polish.alt.bam -t ${task.cpus} -o flye.polish.alt
	cd flye.polish.alt
	\${seqkit} seq -w0 -m${params.pipeline.min_len_collapsed_segment} polished_1.fasta | awk '{if (substr(\$1,1,1)==\">\"){print \$1 \"_alt\"} else {print \$0}}' > assembly.alt.fasta
	cat ../${asm_fa} assembly.alt.fasta > assembly.all.fasta
	"""
}

process hapResAsmPolish_trioBin {

	label 'medium_node'

	input:
		tuple path('mixhap_asm.fasta'), path('mixhap_asm.fasta.fai'), path('mixhap_asm.bed')
		tuple path('hapres_asm.fasta'), path('hapres_asm.fasta.fai'), path('hapres_asm.bed')

		path hapres_asm_including_alt_fa

		tuple path('hap_lr_polish.bam'), path('hap_lr_polish.bam.bai')
		tuple path('hap_lr_polish.vcf.gz'), path('hap_lr_polish.vcf.gz.tbi')

		tuple path('trio_binning_parentA.fq.gz'), path('trio_binning_parentA.gfa.gz'), path('trio_binning_parentA.bfi')
		tuple path('trio_binning_parentB.fq.gz'), path('trio_binning_parentB.gfa.gz'), path('trio_binning_parentB.bfi')

	output:
		tuple path("trio_binning_h0.bed"), path("trio_binning_h1.bed"), path("trio_binning_h2.bed")

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	minimap2=\${MINIMAP2:-${params.tools.minimap2.bin}}
	triobin=\${TRIO_BINNING:-${params.tools.trio_binning.bin}}

	bin_ambiguous_py=\${BIN_AMBIGUOUS_HAP_PY:-${params.python.script.bin_ambiguous}}
	detect_hap_collide_py=\${DETECT_HAP_COLLIDE_PY:-${params.python.script.detect_hap_collide}}

	TASK_MEM=\$(echo -e \"${task.memory}\" | cut -d \" \" -f1);
	MEM_PER_THREADS_SORT=\$(bc -l <<< \"(\${TASK_MEM} / ${task.cpus}) * ${params.tools.samtools.sort.mem_safety_ratio} * 1000\" | awk '{printf(\"%.0f\", \$0)}');

	\${triobin} phase -t ${task.cpus} -v -f -S -U 1 -C ${hapres_asm_including_alt_fa} -r trio_binning -o trio_binning;
	cat trio_binning.tsv | tr '_' '\\t' | sort -k1,1 -k2,2 -k3,3 | awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} \
	{if (\$4==\"alt\"){print \$1 \"_\" \$2 \"_\" \$3 \"_\" \$4,\$5,\$6,\$7,\$8} else {print \$1 \"_\" \$2 \"_\" \$3,\$4,\$5,\$6,\$7}}' > trio_binning.sorted.tsv;
	mv -f trio_binning.sorted.tsv trio_binning.tsv;
	awk 'BEGIN {LN=1; PREV=\"\"} { if (LN==1){LN1=\$0; PREV=\$1; LN=0} else if (substr(\$1,1,length(\$1)-4)==PREV) {print LN1; print \$0; LN=1} else {LN1=\$0; PREV=\$1; LN=0}}' trio_binning.tsv > trio_binning.paired.tsv;
	awk 'BEGIN {LN=1; PREV=\"\"} { if (LN==1){LN1=\$0; PREV=\$1; LN=0} else if (substr(\$1,1,length(\$1)-4)==PREV) {LN=1} else {print LN1; LN1=\$0; PREV=\$1; LN=0}} END {if (LN==0){print LN1}}' trio_binning.tsv > trio_binning.unpaired.tsv;
	awk '{SUM=\$3+\$4; if (NR%2==1){CONTIGR=\$0; FORCED_PHASING_R=0; LOW_RATIO_UNIQ_R=0; if (SUM==0) {HR=0; COUNTR=0} else if (\$3>=\$4) {HR=1; COUNTR=\$3} else {HR=2; COUNTR=\$4}; \
	if (SUM!=0){ if (COUNTR/SUM<0.8){FORCED_PHASING_R=1}; if (COUNTR/(\$2+SUM)<0.0005){LOW_RATIO_UNIQ_R=1}}} else {CONTIGA=\$0; FORCED_PHASING_A=0; LOW_RATIO_UNIQ_A=0; if (SUM==0) {HA=0; COUNTA=0} \
	else if (\$3>=\$4) {HA=1; COUNTA=\$3} else {HA=2; COUNTA=\$4}; if (SUM!=0) {if (COUNTA/SUM<0.8){FORCED_PHASING_A=1}; if (COUNTA/(\$2+SUM)<0.0005){LOW_RATIO_UNIQ_A=1}}; \
	if (HR!=HA){if ((HR==0) && (HA!=0)) {print CONTIGA} else {print CONTIGR}} else if (((FORCED_PHASING_R==1) && (FORCED_PHASING_A==0)) || ((LOW_RATIO_UNIQ_R==1) && (COUNTA>=1.2*COUNTR)) || \
	((FORCED_PHASING_R==0) && (FORCED_PHASING_A==0) && (COUNTA>=2*COUNTR))) {print CONTIGA} else {print CONTIGR}}}' trio_binning.paired.tsv > trio_binning.best_binning.tsv;
	cat <({ grep -v alt trio_binning.best_binning.tsv || true; }) <({ grep alt trio_binning.best_binning.tsv || true; } | awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"} {print substr(\$1,1,length(\$1)-4),\$2,\$4,\$3,\$5}') trio_binning.unpaired.tsv | \
	sort > trio_binning.best_binning.sorted.tsv;
	mv -f trio_binning.best_binning.sorted.tsv trio_binning.best_binning.tsv;

	mkdir -p binAmbiguous; cd binAmbiguous;
	awk '{SUM=\$3+\$4; CONTIGR=\$1; FORCED_PHASING_R=0; LOW_RATIO_UNIQ_R=0; if (SUM==0) {HR=0; COUNTR=0} else if (\$3>=\$4) {HR=1; COUNTR=\$3} else {HR=2; COUNTR=\$4}; \
	if (SUM!=0){ if (COUNTR/SUM<=0.6){FORCED_PHASING_R=1}; if (COUNTR/(\$2+SUM)<0.00025){LOW_RATIO_UNIQ_R=1}}; \
	if ((HR!=0) && (FORCED_PHASING_R==0) && (LOW_RATIO_UNIQ_R==0)) {print CONTIGR \"\\t\" HR} else {print CONTIGR \"\\t0\"}}' ../trio_binning.best_binning.tsv > hap_bin.conf.list;

	awk '{SUM=\$3+\$4; CONTIGR=\$1; FORCED_PHASING_R=0; LOW_RATIO_UNIQ_R=0; if (SUM==0) {HR=0; COUNTR=0} else if (\$3>=\$4) {HR=1; COUNTR=\$3} else \
	{HR=2; COUNTR=\$4}; print CONTIGR \"\t\" HR}' ../trio_binning.best_binning.tsv > hap_bin.unconf.list;

	python \${bin_ambiguous_py} -t ${task.cpus} -b ../hap_lr_polish.bam -v ../hap_lr_polish.vcf.gz -B hap_bin.conf.list > hap_bin.conf.1.list;
	awk '{if (\$2!=0){print \$0}}' hap_bin.conf.1.list | sort -k1,1 > hap_bin.conf.1.binned.list;
	awk '{if (\$2==0){print \$0}}' hap_bin.conf.1.list | sort -k1,1 > hap_bin.conf.1.unbinned.list;
	cat hap_bin.conf.1.binned.list <(join hap_bin.conf.1.unbinned.list hap_bin.unconf.list | awk '{print \$1 \"\\t\" \$4 \"\\t\" \$3}') | sort -k1,1 > hap_bin.conf.2.all.list;

	mkdir -p hap2mixed; cd hap2mixed;
	\${minimap2} -t ${task.cpus} ${params.tools.minimap2.param.asm} ../../mixhap_asm.fasta ../../hapres_asm.fasta > asm_hap2mixed.sam;
	\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M asm_hap2mixed.sam > asm_hap2mixed.bam;
	if [ ! -s asm_hap2mixed.bam ]; then echo \"File asm_hap2mixed.bam does not exist or is empty\" 1>&2; exit 1; fi;
	\${samtools} quickcheck asm_hap2mixed.bam;
	if [ ! \$? -eq 0 ]; then echo \"File asm_hap2mixed.bam is malformed\" 1>&2; exit 1; fi;
	\${samtools} index -@ ${task.cpus} asm_hap2mixed.bam; rm -rf asm_hap2mixed.sam;
	python \${detect_hap_collide_py} -t ${task.cpus} -b asm_hap2mixed.bam -p ../hap_bin.conf.2.all.list | sort > trio_binning.hap_collision.list;
	cut -f1 trio_binning.hap_collision.list > trio_binning.hap_collision.nohap.list;

	cd ../..

	for HAP in {0..2}; do awk -v hap=\${HAP} '{if (\$2==hap){print \$1}}' binAmbiguous/hap_bin.conf.2.all.list > h\${HAP}.1.list; done;
	for HAP in {0..2}; do comm -23 h\${HAP}.1.list binAmbiguous/hap2mixed/trio_binning.hap_collision.nohap.list > h\${HAP}.list; done;

	for HAP in {0..2}
	do
		awk -v hap=\${HAP} '{if (\$3==hap){print \$1}}' binAmbiguous/hap2mixed/trio_binning.hap_collision.list >> h\${HAP}.list
		sort h\${HAP}.list > h\${HAP}.tmp.list; mv -f h\${HAP}.tmp.list h\${HAP}.list
	done

	for HAP in {0..2}; do join <(sort -k1,1 hapres_asm.bed) h\${HAP}.list | sort -k1,1 -k2,2n > trio_binning_h\${HAP}.bed; done;
	"""
}

process hapResAsmPolish_extractPhased {

	label 'medium_node'

	input:
		tuple path('h0.bed'), path('h1.bed'), path('h2.bed')
		tuple path('lr.hap.bam'), path('lr.hap.bam.bai')
		tuple path('lr.hap.vcf.gz'), path('lr.hap.vcf.gz.tbi')

	output:
		tuple val('H1'), path('h1.hap.fastq.gz')
		tuple val('H2'), path('h2.hap.fastq.gz')

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	pigz=\${PIGZ_BIN:-${params.tools.pigz.bin}}

	extract_ref_reads_py=\${EXTRACT_REF_READS_PY:-${params.python.script.extract_ref_reads}}

	\${samtools} view -@ ${task.cpus} -b -M -L h1.bed -o h1.bam lr.hap.bam
	\${samtools} index -@ ${task.cpus} h1.bam
	\${samtools} view -@ ${task.cpus} -b -M -L h2.bed -o h2.bam lr.hap.bam
	\${samtools} index -@ ${task.cpus} h2.bam

	python \${extract_ref_reads_py} -t ${task.cpus} -b h1.bam -v lr.hap.vcf.gz -p | \${pigz} -p ${task.cpus} -c > h1.hap.fastq.gz
	python \${extract_ref_reads_py} -t ${task.cpus} -b h2.bam -v lr.hap.vcf.gz -p | \${pigz} -p ${task.cpus} -c > h2.hap.fastq.gz
	python \${extract_ref_reads_py} -t ${task.cpus} -b h1.bam -v lr.hap.vcf.gz -p --alt | \${pigz} -p ${task.cpus} -c >> h2.hap.fastq.gz
	python \${extract_ref_reads_py} -t ${task.cpus} -b h2.bam -v lr.hap.vcf.gz -p --alt | \${pigz} -p ${task.cpus} -c >> h1.hap.fastq.gz

	\${pigz} -p ${task.cpus} -t h1.hap.fastq.gz
	if [ ! \$? -eq 0 ]; then echo \"File h1.haplotigs.fastq.gz is malformed\" 1>&2; exit 1; fi
	\${pigz} -p ${task.cpus} -t h2.hap.fastq.gz
	if [ ! \$? -eq 0 ]; then echo \"File h2.haplotigs.fastq.gz is malformed\" 1>&2; exit 1; fi

	rm -rf h*.bam*
	"""
}

process hapResAsmPolish_getCollapsedHom {

	label 'medium_node'

	input:
		tuple path('h0.bed'), path('h1.bed'), path('h2.bed')
		tuple path('lr.hap.phasedSupp2prim.bam'), path('lr.hap.phasedSupp2prim.bam.bai')
		tuple path('lr.hap.vcf.gz'), path('lr.hap.vcf.gz.tbi')

	output:
		tuple val('H1'), path('h1.het_or_collapsed.fastq.gz'), path('h1.homozygous.fastq.gz'), emit: fastq_h1
		tuple val('H2'), path('h2.het_or_collapsed.fastq.gz'), path('h2.homozygous.fastq.gz'), emit: fastq_h2
		path('cov_stdev.tsv'), emit: cov

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	bcftools=\${BCFTOOLS:-${params.tools.bcftools.bin}}
	bedtools=\${BEDTOOLS:-${params.tools.bedtools.bin}}
	pigz=\${PIGZ_BIN:-${params.tools.pigz.bin}}

	extract_reads_py=\${EXTRACT_READS_PY:-${params.python.script.extract_reads}}
	extract_ref_reads_py=\${EXTRACT_REF_READS_PY:-${params.python.script.extract_ref_reads}}
	find_coll_hom_py=\${FIND_COLL_HOM_PY:-${params.python.script.find_coll_hom}}

	python \${extract_reads_py} -t ${task.cpus} -b lr.hap.phasedSupp2prim.bam -i HP:i:0 -I HP -B lr.hap.unphased.bam
	\${samtools} index -@ ${task.cpus} lr.hap.unphased.bam
	bcftools query -i \"(type=='snp') & (PS>=0)\" -f '%CHROM\t%POS0\t%END\n' lr.hap.vcf.gz > lr.hap.phased.bed
	\${samtools} depth -@ ${task.cpus} -a -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} -b lr.hap.phased.bed lr.hap.phasedSupp2prim.bam | \
	awk '{print \$1 \"\\t\" \$2 \"\\t\" (\$2+1) \"\\t\" \$3}' > lr.hap.phased.bed.tmp
	mv -f lr.hap.phased.bed.tmp lr.hap.phased.bed
	awk 'BEGIN {i=1} {SUM+=\$4; COUNT+=1; SAMP[i]=\$4; i+=1} END {MEAN=SUM/COUNT; SUMSQ=0; for (i=1; i<=COUNT; i++) {SUMSQ+=(SAMP[i]-MEAN)^2}; print MEAN \"\\t\" sqrt(SUMSQ/COUNT)}' lr.hap.phased.bed > cov_stdev.tsv
	COV=\$(awk '{print \$1}' cov_stdev.tsv); AAD_COV=\$(awk '{print \$2}' cov_stdev.tsv)
	python ${params.python.script.extract_ref_reads} -t ${task.cpus} -b lr.hap.phasedSupp2prim.bam -v lr.hap.vcf.gz -B lr.hap.ref.bam --only_include_phased --include_supplementary;
	\${samtools} index -@ ${task.cpus} lr.hap.ref.bam
	join <(\${samtools} depth -@ ${task.cpus} -aa -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} -b lr.hap.phased.bed lr.hap.ref.bam | \
	awk '{print \$1 \"\\t\" \$2 \"\\t\" \$3}' | sort -k1,1 -k2,2n) \
	<(\${samtools} depth -@ ${task.cpus} -aa -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} -b lr.hap.phased.bed lr.hap.unphased.bam | \
	awk '{print \$1 \"\\t\" \$2 \"\\t\" (\$3/2)}' | sort -k1,1 -k2,2n) | \
	tr ' ' '\\t' | awk '{if (\$1!=CONTIG){if (CONTIG!=\"\") {MEAN=SUM/COUNT; SUMSQ=0; for (i=1; i<=COUNT; i++) {SUMSQ+=(SAMP[i]-MEAN)^2}; print CONTIG \"\\t\" MEAN \"\\t\" sqrt(SUMSQ/COUNT)}; \
	CONTIG=\$1; SUM=\$3+\$5; COUNT=1; i=1} else {SUM+=\$3+\$5; COUNT+=1}; SAMP[i]=\$3+\$5; i+=1} END {if (CONTIG!=\"\") {MEAN=SUM/COUNT; SUMSQ=0; for (i=1; i<=COUNT; i++) {SUMSQ+=(SAMP[i]-MEAN)^2}; \
	print CONTIG \"\\t\" MEAN \"\\t\" sqrt(SUMSQ/COUNT)};}' > lr.hap.phased.ref.cov
	rm -rf lr.hap.ref.bam*

	python \${extract_ref_reads_py} -t ${task.cpus} -b lr.hap.phasedSupp2prim.bam -v lr.hap.vcf.gz -B lr.hap.alt.bam --alt --only_include_phased --include_supplementary
	\${samtools} index -@ ${task.cpus} lr.hap.alt.bam

	join <(\${samtools} depth -@ ${task.cpus} -aa -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} -b lr.hap.phased.bed lr.hap.alt.bam | \
	awk '{print \$1 \"\\t\" \$2 \"\\t\" \$3}' | sort -k1,1 -k2,2n) \
	<(\${samtools} depth -@ ${task.cpus} -aa -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} -b lr.hap.phased.bed lr.hap.unphased.bam | \
	awk '{print \$1 \"\\t\" \$2 \"\\t\" (\$3/2)}' | sort -k1,1 -k2,2n) | \
	tr ' ' '\t' | awk '{if (\$1!=CONTIG){if (CONTIG!=\"\") {MEAN=SUM/COUNT; SUMSQ=0; for (i=1; i<=COUNT; i++) {SUMSQ+=(SAMP[i]-MEAN)^2}; print CONTIG \"\\t\" MEAN \"\\t\" sqrt(SUMSQ/COUNT)}; \
	CONTIG=\$1; SUM=\$3+\$5; COUNT=1; i=1} else {SUM+=\$3+\$5; COUNT+=1}; SAMP[i]=\$3+\$5; i+=1} END {if (CONTIG!=\"\") {MEAN=SUM/COUNT; SUMSQ=0; for (i=1; i<=COUNT; i++) {SUMSQ+=(SAMP[i]-MEAN)^2}; \
	print CONTIG \"\\t\" MEAN \"\\t\" sqrt(SUMSQ/COUNT)};}' > lr.hap.phased.alt.cov
	rm -rf lr.hap.alt.bam*

	awk '{if (\$1!=CONTIG){if (CONTIG!=\"\") {MEAN=SUM/COUNT; SUMSQ=0; for (i=1; i<=COUNT; i++) {SUMSQ+=(SAMP[i]-MEAN)^2}; print CONTIG \"\\t\" MEAN \"\\t\" sqrt(SUMSQ/COUNT)}; CONTIG=\$1; SUM=\$4; COUNT=1; i=1} \
	else {SUM+=\$4; COUNT+=1}; SAMP[i]=\$4; i+=1} END {if (CONTIG!=\"\") {MEAN=SUM/COUNT; SUMSQ=0; for (i=1; i<=COUNT; i++) {SUMSQ+=(SAMP[i]-MEAN)^2}; print CONTIG \"\\t\" MEAN \"\\t\" sqrt(SUMSQ/COUNT)};}' \
	lr.hap.phased.bed > lr.hap.phased.all.cov

	join <(sort -k1 lr.hap.phased.all.cov) <(join <(sort -k1 lr.hap.phased.ref.cov) <(sort -k1 lr.hap.phased.alt.cov) | tr ' ' '\\t') | tr ' ' '\\t' > lr.hap.phased.cov
	cat lr.hap.phased.cov <(\${samtools} depth -@ ${task.cpus} -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} lr.hap.phasedSupp2prim.bam) | \
	awk -v MIN_COV_GLOB=\$(bc -l <<< \"\${COV}*${params.pipeline.ratio_cov.lower_bound.strict}\" | awk '{printf(\"%.0f\", \$0)}') -v NB_PHASED_HAP=\$(wc -l lr.hap.phased.cov | cut -d ' ' -f1) \
	'BEGIN {CONTIG=\"\"; POS_S=-1; POS_E=-1} {if (NR<=NB_PHASED_HAP){MIN_COV_LOC[\$1]=\$2-(\$2-\$4)/2} else {MIN_COV=MIN_COV_GLOB; if (\$1 in MIN_COV_LOC){MIN_COV=MIN_COV_LOC[\$1]}; \
	if (\$3>=MIN_COV) {if ((\$1==CONTIG) && (\$2==POS_E)) {POS_E+=1} else { if ((POS_S!=-1) && (POS_E-POS_S>=${params.pipeline.min_len_collapsed_segment})) {print CONTIG \"\\t\" (POS_S-1) \"\\t\" (POS_E-1)}; CONTIG=\$1; POS_S=\$2; POS_E=\$2+1}}}} \
	END {if ((POS_S!=-1) && (POS_E-POS_S>=${params.pipeline.min_len_collapsed_segment})) {print CONTIG \"\\t\" (POS_S-1) \"\\t\" (POS_E-1)}}' > lr.hap.phased.cov.bed

	python \${extract_reads_py} -t ${task.cpus} -b lr.hap.phasedSupp2prim.bam -i HP:i:1 -i HP:i:2 -B lr.hap.phased.bam
	\${samtools} index -@ ${task.cpus} lr.hap.phased.bam

	\${samtools} depth -@ ${task.cpus} -aa -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} lr.hap.phased.bam | \
	awk 'BEGIN {CONTIG=\"\"; POS_S=-1; POS_E=-1} {if (\$3>=${params.pipeline.variant_filter.min_depth}) {if ((\$1==CONTIG) && (\$2==POS_E)) {POS_E+=1} \
	else { if ((POS_S!=-1) && (POS_E-POS_S>=${params.pipeline.min_len_collapsed_segment})) {print CONTIG \"\\t\" (POS_S-1) \"\\t\" (POS_E-1)}; CONTIG=\$1; POS_S=\$2; POS_E=\$2+1;}}} \
	END {if ((POS_S!=-1) && (POS_E-POS_S>=${params.pipeline.min_len_collapsed_segment})) {print CONTIG \"\\t\" (POS_S-1) \"\\t\" (POS_E-1)}}' > lr.hap.phased.${params.pipeline.variant_filter.min_depth}x.bed
	rm -rf lr.hap.phased.bam*

	\${bedtools} merge -i lr.hap.phased.${params.pipeline.variant_filter.min_depth}x.bed -d ${params.pipeline.min_len_collapsed_segment} | \
	sort -k1,1 -k2,2n > lr.hap.phased.${params.pipeline.variant_filter.min_depth}x.merge.bed

	\${bedtools} intersect -wa -a lr.hap.phased.cov.bed -b lr.hap.phased.${params.pipeline.variant_filter.min_depth}x.merge.bed | sort -k1,1 -k2,2n | uniq > lr.hap.phased.cov.1.bed
	\${bedtools} subtract -A -a lr.hap.phased.cov.bed -b lr.hap.phased.${params.pipeline.variant_filter.min_depth}x.merge.bed | sort -k1,1 -k2,2n | uniq > lr.hap.phased.cov.2.bed

	python \${find_coll_hom_py} -t ${task.cpus} -b lr.hap.phasedSupp2prim.bam -r lr.hap.phased.cov.2.bed >> lr.hap.phased.cov.1.bed
	cat lr.hap.phased.cov <(\${samtools} depth -@ ${task.cpus} -J -Q ${params.pipeline.min_mapq.strict} -G ${params.pipeline.samtools.depth.filter_sec} lr.hap.phasedSupp2prim.bam) | \
	awk -v MAX_COV_GLOB=\$(bc -l <<< \"\${COV}*${params.pipeline.ratio_cov.upper_bound.strict}\" | awk '{printf(\"%.0f\", \$0)}') -v NB_PHASED_HAP=\$(wc -l lr.hap.phased.cov | cut -d ' ' -f1) \
	'BEGIN {CONTIG=\"\"; POS_S=-1; POS_E=-1} {if (NR<=NB_PHASED_HAP){if (\$4>\$6){MAX_COV_LOC[\$1]=\$2+\$4/2} else {MAX_COV_LOC[\$1]=\$2+\$6/2}} else \
	{MAX_COV=MAX_COV_GLOB; if (\$1 in MAX_COV_LOC){MAX_COV=MAX_COV_LOC[\$1]}; if (\$3>=MAX_COV) {if ((\$1==CONTIG) && (\$2==POS_E)) {POS_E+=1} else \
	{ if ((POS_S!=-1) && (POS_E-POS_S>=${params.pipeline.min_len_collapsed_segment})) {print CONTIG \"\\t\" (POS_S-1) \"\\t\" (POS_E-1)}; CONTIG=\$1; POS_S=\$2; POS_E=\$2+1}}}} \
	END {if ((POS_S!=-1) && (POS_E-POS_S>=${params.pipeline.min_len_collapsed_segment})) {print CONTIG \"\\t\" (POS_S-1) \"\\t\" (POS_E-1)}}' > lr.hap.collapsed.bed

	\${samtools} view -@ ${task.cpus} -b -M -L lr.hap.collapsed.bed -o lr.hap.unphased.collapsed.bam lr.hap.unphased.bam
	\${samtools} index -@ ${task.cpus} lr.hap.unphased.collapsed.bam

	\${samtools} bam2fq -@ ${task.cpus} -n lr.hap.unphased.collapsed.bam | awk '{if (NR%4==1){print substr(\$0,2,length(\$0)-1)}}' > lr.hap.collapsed.list
	rm -rf lr.hap.unphased.collapsed.bam*

	\${samtools} view -@ ${task.cpus} -b -M -L lr.hap.phased.cov.1.bed lr.hap.unphased.bam > lr.hap.unphased.hom_candidates.bam
	\${samtools} index -@ ${task.cpus} lr.hap.unphased.hom_candidates.bam

	python \${extract_reads_py} -t ${task.cpus} -b lr.hap.unphased.hom_candidates.bam -L lr.hap.collapsed.list | \${pigz} -p ${task.cpus} -c > lr.hap.homozygous.fastq.gz

	\${pigz}  -d -c -p ${task.cpus} lr.hap.homozygous.fastq.gz | awk '{if (((NR-1)%8)<4){print \$0}}' | \${pigz} -p ${task.cpus} -c > h1.homozygous.fastq.gz
	\${pigz}  -d -c -p ${task.cpus} lr.hap.homozygous.fastq.gz | awk '{if (((NR-1)%8)>=4){print \$0}}' | \${pigz} -p ${task.cpus} -c > h2.homozygous.fastq.gz
	\${pigz}  -d -c -p ${task.cpus} lr.hap.homozygous.fastq.gz | awk '{if (NR%4==1){print substr(\$0,2,length(\$0)-1)}}' > lr.hap.homozygous.list

	python \${extract_reads_py}  -t ${task.cpus} -b lr.hap.unphased.bam -B lr.hap.unphased.het_or_collapsed.bam -L lr.hap.homozygous.list

	\${samtools} index -@ ${task.cpus} lr.hap.unphased.het_or_collapsed.bam

	\${samtools} view -@ ${task.cpus} -b -M -L h1.bed lr.hap.unphased.het_or_collapsed.bam | \
	\${samtools} bam2fq -@ ${task.cpus} -n - | \${pigz} -p ${task.cpus} -c > h1.het_or_collapsed.fastq.gz

	\${samtools} view -@ ${task.cpus} -b -M -L h2.bed lr.hap.unphased.het_or_collapsed.bam | \
	\${samtools} bam2fq -@ ${task.cpus} -n - | \${pigz} -p ${task.cpus} -c > h2.het_or_collapsed.fastq.gz
	"""			
}

process hapResAsmPolish_dualAsm {

	label 'large_node'

	input:
		tuple val(HAP), path('haplotigs.fastq.gz'), path('het_or_collapsed.fastq.gz'), path('homozygous.fastq.gz')

	output:
		tuple val("${HAP}"), path('assembly.filtered.fasta'), path('assembly.filtered.bed'), path("hap.nodup.fastq.gz")

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	flye=\${FLYE:-${params.tools.flye.bin}}
	seqkit=\${SEQKIT:-${params.tools.seqkit.bin}}
	seqtk=\${SEQTK:-${params.tools.seqtk.bin}}
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	pigz=\${PIGZ_BIN:-${params.tools.pigz.bin}}

	cat haplotigs.fastq.gz homozygous.fastq.gz het_or_collapsed.fastq.gz > hap.fastq.gz
	\${seqkit} rmdup -j ${task.cpus} -n -o hap.nodup.fastq hap.fastq.gz
	rm -rf hap.fastq.gz
	\${pigz} -p ${task.cpus} hap.nodup.fastq;

	FLYE_SUBSAMPLE_PARAM=\"\";
	if [ ${params.genome_size} -gt 0 ]; then FLYE_SUBSAMPLE_PARAM=\"--asm-coverage ${params.pipeline.assembly.subsampling_cov} --genome-size ${params.genome_size}\"; fi;
	\${flye} --nano-hq hap.nodup.fastq.gz -o . -t ${task.cpus} --read-error ${params.pipeline.assembly.max_error_rate.strict} --scaffold \${FLYE_SUBSAMPLE_PARAM}
	if [ ! \$? -eq 0 ]; then echo \"Flye has exited prematurely\" 1>&2; exit 1; fi
	ASM_COMPLETED=\$({ grep \"INFO: Final assembly\" <(tail -n 50 flye.log) || true; } | wc -l)
	if [ \${ASM_COMPLETED} -eq 0 ] || [ ! -s assembly.fasta ]; then echo \"Flye log incomplete or empty assembly\" 1>&2; exit 1; fi

	\${seqkit} seq -w0 assembly.fasta > assembly.fasta.tmp
	mv -f assembly.fasta.tmp assembly.fasta
	rm -rf 00-assembly 10-consensus 20-repeat 30-contigger 40-polishing
	tail -n +2 assembly_info.txt | awk '{if (\$7==\"*\"){print \$1}}' > contigs.list
	\${seqtk} subseq assembly.fasta contigs.list | awk '{if (NR%2==1){print \">${HAP}#\" substr(\$0,2,length(\$0)-1)} else {print \$0}}' > assembly.filtered.fasta
	\${samtools} faidx assembly.filtered.fasta
	awk '{print \$1 \"\\t0\\t\" \$2}' assembly.filtered.fasta.fai | sort -k1,1 -k2,2n > assembly.filtered.bed
	"""
}

process hapResAsmPolish_polishDualAsm_1 {

	label 'medium_node'

	input:
		tuple path('h1.assembly.filtered.fasta'), path('h1.assembly.filtered.bed'), path("h1.nodup.fastq.gz")
		tuple path('h2.assembly.filtered.fasta'), path('h2.assembly.filtered.bed'), path("h2.nodup.fastq.gz")
		path('cov_stdev.tsv')
		path('lr.filtered.fq.gz')

	output:
		tuple path('assembly.filtered.fasta'), path('assembly.filtered.fasta.fai'), path('assembly.filtered.bed'), emit: asm
		tuple path('lr.asm.bam'), path('lr.asm.bam.bai'), emit: bam

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	minimap2=\${MINIMAP2:-${params.tools.minimap2.bin}}
	seqtk=\${SEQTK:-${params.tools.seqtk.bin}}

	TASK_MEM=\$(echo -e \"${task.memory}\" | cut -d \" \" -f1);
	MEM_PER_THREADS_SORT=\$(bc -l <<< \"(\${TASK_MEM} / ${task.cpus}) * ${params.tools.samtools.sort.mem_safety_ratio} * 1000\" | awk '{printf(\"%.0f\", \$0)}');

	HAP_COV=\$(awk '{print \$1/2}' cov_stdev.tsv);

	for HAP_A in {1..2}
	do
		HAP_B=2

		if [ \${HAP_A} -eq 2 ]; then HAP_B=1; fi

		\${minimap2} -t ${task.cpus} ${params.tools.minimap2.param.asm} h\${HAP_B}.assembly.filtered.fasta h\${HAP_A}.assembly.filtered.fasta > h\${HAP_A}_map2_h\${HAP_B}.sam
		\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M h\${HAP_A}_map2_h\${HAP_B}.sam > h\${HAP_A}_map2_h\${HAP_B}.bam
		if [ ! -s h\${HAP_A}_map2_h\${HAP_B}.bam ]; then echo \"File h\${HAP_A}_map2_h\${HAP_B}.bam does not exist or is empty\" 1>&2; exit 1; fi
		\${samtools} quickcheck h\${HAP_A}_map2_h\${HAP_B}.bam
		if [ ! \$? -eq 0 ]; then echo \"File h\${HAP_A}_map2_h\${HAP_B}.bam is malformed\" 1>&2; exit 1; fi
		\${samtools} index -@ ${task.cpus} h\${HAP_A}_map2_h\${HAP_B}.bam; rm -rf h\${HAP_A}_map2_h\${HAP_B}.sam

		\${samtools} view -F 2048 -@ ${task.cpus} --tag de h\${HAP_A}_map2_h\${HAP_B}.bam | \
		awk '{ERR=\"1.0\"; for (i=6; i<=NF; i+=1){ if (substr(\$i,1,5)==\"de:f:\"){ERR=substr(\$i,6,length(\$i)-5); if (ERR<=0.0005) {print \$1}; break;}}}' | \
		sort | uniq > h\${HAP_A}_map2_h\${HAP_B}.high_sim.list
		join h\${HAP_A}_map2_h\${HAP_B}.high_sim.list h\${HAP_A}.assembly.filtered.bed > h\${HAP_A}_map2_h\${HAP_B}.high_sim.bed
		rm -rf h\${HAP_A}_map2_h\${HAP_B}.bam*

		\${minimap2} -t ${task.cpus} -ax map-hifi -Y h\${HAP_A}.assembly.filtered.fasta h\${HAP_A}.nodup.fastq.gz > lr.h\${HAP_A}.sam
		\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M lr.h\${HAP_A}.sam > lr.h\${HAP_A}.bam
		if [ ! -s lr.h\${HAP_A}.bam ]; then echo \"File lr.h\${HAP_A}.bam does not exist or is empty\" 1>&2; exit 1; fi
		\${samtools} quickcheck lr.h\${HAP_A}.bam
		if [ ! \$? -eq 0 ]; then echo \"File lr.h\${HAP_A}.bam is malformed\" 1>&2; exit 1; fi
		\${samtools} index -@ ${task.cpus} lr.h\${HAP_A}.bam; rm -rf lr.h\${HAP_A}.sam

		\${samtools} depth -@ ${task.cpus} -aa -J -Q ${params.pipeline.min_mapq.strict} -G 1540 -b h\${HAP_A}_map2_h\${HAP_B}.high_sim.bed lr.h\${HAP_A}.bam | \
		awk -v cov=\"\${HAP_COV}\" 'BEGIN {CONTIG=\"\"; TOT_BP=0; TOT_COV_BP=0} {if (\$1!=CONTIG){if ((CONTIG!=\"\") && (TOT_COV_BP/TOT_BP>=0.9)) {print CONTIG}; CONTIG=\$1; TOT_BP=0; TOT_COV_BP=0}; \
		if (\$3<cov*0.6){TOT_COV_BP+=1}; TOT_BP+=1} END {if ((CONTIG!=\"\") && (TOT_COV_BP/TOT_BP>=0.9)) {print CONTIG}}' > h\${HAP_A}_map2_h\${HAP_B}.high_sim.low_cov.list
		comm -23 <(cut -f1 h\${HAP_A}.assembly.filtered.bed | sort) <(sort h\${HAP_A}_map2_h\${HAP_B}.high_sim.low_cov.list) > h\${HAP_A}_map2_h\${HAP_B}.high_sim.high_cov.list
		\${seqtk} subseq h\${HAP_A}.assembly.filtered.fasta h\${HAP_A}_map2_h\${HAP_B}.high_sim.high_cov.list > assembly.filtered.h\${HAP_A}.fasta
		rm -rf lr.h\${HAP_A}.bam*;
	done

	cat assembly.filtered.h1.fasta assembly.filtered.h2.fasta > assembly.filtered.fasta
	\${samtools} faidx assembly.filtered.fasta
	awk '{print \$1 \"\\t0\\t\" \$2}' assembly.filtered.fasta.fai | sort -k1,1 -k2,2n > assembly.filtered.bed

	\${minimap2} -t ${task.cpus} -ax map-hifi -Y -I 6G assembly.filtered.fasta lr.filtered.fq.gz > lr.asm.sam
	\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M lr.asm.sam > lr.asm.bam
	if [ ! -s lr.asm.bam ]; then echo \"File lr.asm.bam does not exist or is empty\" 1>&2; exit 1; fi
	\${samtools} quickcheck lr.asm.bam
	if [ ! \$? -eq 0 ]; then echo \"File lr.asm.bam is malformed\" 1>&2; exit 1; fi
	\${samtools} index -@ ${task.cpus} lr.asm.bam
	rm -rf lr.asm.sam
	"""
}

process hapResAsmPolish_polishDualAsm_2 {

	label 'medium_node'

	input:
		tuple path(asm_fa), path(asm_fai), path(asm_bed)
		tuple path('lr.asm.vcf.gz'), path('lr.asm.vcf.gz.tbi')
		path('cov_stdev.tsv')

	output:
		tuple path('assembly.filtered.fasta'), path('assembly.filtered.fasta.fai'), path('assembly.filtered.bed')

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	seqtk=\${SEQTK:-${params.tools.seqtk.bin}}

	polish_snp_py=\${POLISH_SNP_PY:-${params.python.polish_snp}}

	python \${polish_snp_py} -a ${asm_fa} -v lr.asm.vcf.gz -g 20 -d \$(awk '{printf(\"%.0f\", \$1)}' cov_stdev.tsv) > assembly.filtered.tmp.fasta
	\${seqtk} seq -L ${params.pipeline.min_len_contig} assembly.filtered.tmp.fasta > assembly.filtered.fasta
	rm -rf assembly.filtered.tmp.fasta

	\${samtools} faidx assembly.filtered.fasta
	awk '{print \$1 \"\\t0\\t\" \$2}' assembly.filtered.fasta.fai | sort -k1,1 -k2,2n > assembly.filtered.bed
	"""	
}

process hapResAsmPolish_polishDualAsm_3 {

	label 'medium_node'

	publishDir "${params.out_dir}"

	input:
		tuple path(asm_fa), path(asm_fai), path(asm_bed)
		tuple path('sr.bam'), path('sr.bam.bai')

	output:
		tuple path('H1.polished.fasta'), path('H2.polished.fasta')

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	bcftools=\${BCFTOOLS:-${params.tools.bcftools.bin}}
	seqtk=\${SEQTK:-${params.tools.seqtk.bin}}

	mkdir -p polishing

	cut -f1 ${asm_fai} | \
	parallel -j ${task.cpus} \"\
		\${samtools} view -b sr.bam {} > polishing/{#}.bam; \
		\${samtools} index polishing/{#}.bam; \
		\${seqtk} subseq ${asm_fa} <(echo -e {}) > polishing/{#}.fasta;
		\${bcftools} mpileup -Ou -f ${asm_fa} polishing/{#}.bam | \${bcftools} call -mv -Ou | \${bcftools} view --exclude 'GT==\\"het\\"' -Oz -o polishing/{#}.vcf.gz; \
		tabix -p vcf polishing/{#}.vcf.gz; \
		cat polishing/{#}.fasta | \${bcftools} consensus polishing/{#}.vcf.gz > polishing/{}.fasta; \
		rm -rf polishing/{#}.*; \
	\"

	cat polishing/H1#*.fasta > H1.polished.fasta
	cat polishing/H2#*.fasta > H2.polished.fasta

	rm -rf polishing/
	"""
}

workflow {

	if (!params.proband_lr_bam_in && !params.proband_lr_fq_in) error("PROBAND: No corrected long reads in input.")
	if (params.proband_lr_bam_in && params.proband_lr_fq_in) error("PROBAND: Input corrected long reads can be in BAM or FASTQ format but not both.")

	if (!params.proband_sr_bam_in && !params.proband_sr_fq_in) error("PROBAND: No short reads in input.")
	if (params.proband_sr_bam_in && params.proband_sr_fq_in) error("PROBAND: Input short reads can be provided in BAM or FASTQ in input but not both.")

	if (!params.father_sr_bam_in && !params.father_sr_fq_in) error("FATHER: No short reads in input.")
	if (params.father_sr_bam_in && params.father_sr_fq_in) error("FATHER: Input short reads can be provided in BAM or FASTQ in input but not both.")

	if (!params.mother_sr_bam_in && !params.mother_sr_fq_in) error("MOTHER: No short reads in input.")
	if (params.mother_sr_bam_in && params.mother_sr_fq_in) error("MOTHER: Input short reads can be provided in BAM or FASTQ in input but not both.")

	proband_lr_bam = params.proband_lr_bam_in ? [params.proband_lr_bam_in] : []
	proband_lr_fq = params.proband_lr_fq_in ? [params.proband_lr_fq_in] : []

    proband_sr_fq = params.proband_sr_fq_in ? Channel.fromPath(params.proband_sr_fq_in) : extractPairedIllumina(Channel.fromPath(params.proband_sr_bam_in))

    parents_sr = Channel.of(
    				['father', params.father_sr_bam_in ? [params.father_sr_bam_in] : [], params.father_sr_fq_in ? [params.father_sr_fq_in] : []],
    				['mother', params.mother_sr_bam_in ? [params.mother_sr_bam_in] : [], params.mother_sr_fq_in ? [params.mother_sr_fq_in] : []]
    			)

    parents_sr_dbg = extractIllumina_buildDBG(parents_sr)

	parents_sr_dbg.flatten().map { file -> tuple(file.name.substring(0,6), file) }.groupTuple().branch {
        paternal: it[0] == "father"
        		return it[1]
        maternal: it[0] == "mother"
        		return it[1]
    }.set { parents_sr_dbg_bin }

	sr_chunks_fq = splitPairedIllumina(proband_sr_fq).flatten() // Split a paired Illumina FASTQ into chunks

	lr_filt_fq = filterONT(proband_lr_bam, proband_lr_fq) // Filter input corrected long reads based on corrected base quality
	mixhap_asm = mixHapAsm(lr_filt_fq) // Mixed haplotype assembly from corrected filtered long reads

	lr_map2mixhap_asm = mapLRtoAsm_1(lr_filt_fq, mixhap_asm) // Remap filtered corrected long reads to mix haplotype assembly

	filt_mixhap_asm_call = mixHapAsm_filter_call(lr_map2mixhap_asm.bam, mixhap_asm)
	filt_mixhap_asm_phase = mixHapAsm_filter_phase(mixhap_asm, lr_map2mixhap_asm.bam, filt_mixhap_asm_call, false)
	filt_mixhap_asm = mixHapAsm_filter(filt_mixhap_asm_phase.hap_bam, mixhap_asm, filt_mixhap_asm_call, lr_map2mixhap_asm.cov) // Filter haplotigs from mixed haplotype assembly

	sr_chunks_fq.combine(filt_mixhap_asm).set { sr_chunks_filt_mixhap_asm }

	sr_filt_mixhap_asm_chunks_bam = mapPairedIllumina_correct(sr_chunks_filt_mixhap_asm, true) // Map the splitted Illumina chunks to the filtered mixed haplotype assembly

	sr_filt_mixhap_asm_chunks_bam.multiMap { it ->
		bam: it[0]
		bam_bai: it[1]
	}
	.set { sr_filt_mixhap_asm_chunks_bam_split } // List of (bam, bam.bai) -> [[all bam], [all bam.bai]]

	sr_filt_mixhap_asm_merge_bam = mergePairedIlluminaBAM_correct(sr_filt_mixhap_asm_chunks_bam_split.bam.collect(), sr_filt_mixhap_asm_chunks_bam_split.bam_bai.collect()) // Merge the Illumina chunks to one single BAM file

	lr_filt_map2mixhap_asm = mapLRtoAsm_2(lr_filt_fq, filt_mixhap_asm) // Map the filtered long reads to the filtered mixed haplotype assembly

	local_lr_corr_chunks_prep = prepLocalCorrection(lr_filt_map2mixhap_asm.bam, sr_filt_mixhap_asm_merge_bam) // Prep filtered ONT reads for local correction with short reads

	local_lr_corr_chunks_prep.flatten().map { file -> tuple(file.name.substring(13,22), file) }.groupTuple().set { local_lr_corr_chunks_prep_bin } // What a shitshow

	local_lr_corr_chunks_fq = localCorrection(local_lr_corr_chunks_prep_bin, lr_filt_map2mixhap_asm.bam, sr_filt_mixhap_asm_merge_bam).collect() // Perform local correction for each chunk
	local_lr_corr_merge_fq = mergeLocalCorrection(local_lr_corr_chunks_fq, lr_filt_map2mixhap_asm.bam) // Merge the corrected chunks

	hapres_asm = hapResAsm(local_lr_corr_merge_fq) // Haplotype-resolved assembly of the locally corrected reads
	hapres_filt_asm = hapResAsm_filter(local_lr_corr_merge_fq, hapres_asm) // Basic haplotig filtering based on min coverage and length

	split_ps_1 = hapResAsm_polish1_1(local_lr_corr_merge_fq, hapres_filt_asm, lr_filt_map2mixhap_asm.bam, lr_filt_map2mixhap_asm.cov)
	split_ps_2 = var_CallFilterPhase_1(hapres_filt_asm, split_ps_1.bam, split_ps_1.bed, true)
	split_ps_3 = hapResAsm_polish1_2(hapres_filt_asm, split_ps_1.bam, split_ps_2.ps_bed, lr_filt_map2mixhap_asm.cov)
	split_ps = var_CallFilterPhase_2(split_ps_3.asm, split_ps_3.bam, split_ps_3.bed, false)

	hapres_filt_polish_asm = hapResAsm_polish2(split_ps_3.asm, split_ps.hap_bam, split_ps.hap_vcf)

	hapres_filt_polish_asm_phasing_filt_bam = hapResAsmPolish_phase_mapFilter(local_lr_corr_merge_fq, hapres_filt_polish_asm)
	hapres_filt_polish_asm_phasing_filt_vcf = hapResAsmPolish_phase_varCall(hapres_filt_polish_asm, hapres_filt_polish_asm_phasing_filt_bam)
	hapres_filt_polish_asm_phasing_filt_vcf_filt = hapResAsmPolish_phase_varFilter(hapres_filt_polish_asm_phasing_filt_bam, hapres_filt_polish_asm_phasing_filt_vcf, lr_filt_map2mixhap_asm.cov)
	hapres_filt_polish_asm_phasing_filt_phased1 = hapResAsmPolish_phase_varPhase1(hapres_filt_polish_asm, hapres_filt_polish_asm_phasing_filt_bam, hapres_filt_polish_asm_phasing_filt_vcf_filt, false)
	hapres_filt_polish_asm_phasing_filt_phased1_filt = hapResAsmPolish_phase_varPhasedFilter(hapres_filt_polish_asm_phasing_filt_vcf, hapres_filt_polish_asm_phasing_filt_phased1.hap_bam, hapres_filt_polish_asm_phasing_filt_phased1.hap_vcf)
	hapres_filt_polish_asm_phasing_filt_phased2 = hapResAsmPolish_phase_varPhase2(hapres_filt_polish_asm, hapres_filt_polish_asm_phasing_filt_bam, hapres_filt_polish_asm_phasing_filt_phased1_filt, false)
	hapres_filt_polish_asm_phasing_filt_phased2_bam = hapResAsmPolish_phase_selectBestPhasedAlignRec(hapres_filt_polish_asm_phasing_filt_phased2.hap_bam)

	hapres_filt_polish_asm_withAlts = hapResAsmPolish_polishAlt(hapres_filt_polish_asm, hapres_filt_polish_asm_phasing_filt_phased2_bam, hapres_filt_polish_asm_phasing_filt_phased2.hap_vcf)

	triobin = hapResAsmPolish_trioBin(	filt_mixhap_asm, hapres_filt_polish_asm, hapres_filt_polish_asm_withAlts, hapres_filt_polish_asm_phasing_filt_phased2.hap_bam,
										hapres_filt_polish_asm_phasing_filt_phased2.hap_vcf, parents_sr_dbg_bin.paternal, parents_sr_dbg_bin.maternal)

	lr_hap_fq = hapResAsmPolish_extractPhased(triobin, hapres_filt_polish_asm_phasing_filt_phased2_bam, hapres_filt_polish_asm_phasing_filt_phased2.hap_vcf)
	lr_het_hom_coll_fq = hapResAsmPolish_getCollapsedHom(triobin, hapres_filt_polish_asm_phasing_filt_phased2_bam, hapres_filt_polish_asm_phasing_filt_phased2.hap_vcf)

	dual_asm = hapResAsmPolish_dualAsm(lr_hap_fq.mix().join(lr_het_hom_coll_fq.fastq_h1.mix(lr_het_hom_coll_fq.fastq_h2)))

	dual_asm.groupTuple().branch {
        h1: it[0] == "H1"
        	return it.tail()
        h2: it[0] == "H2"
        	return it.tail()
    }.set { dual_asm_bin }

	dual_asm_polishCovErr = hapResAsmPolish_polishDualAsm_1(dual_asm_bin.h1, dual_asm_bin.h2, lr_het_hom_coll_fq.cov, lr_filt_fq)
	dual_asm_polishVarCall = hapResAsmPolish_polishDualAsm_call(dual_asm_polishCovErr.bam, dual_asm_polishCovErr.asm)
	dual_asm_polishSNP = hapResAsmPolish_polishDualAsm_2(dual_asm_polishCovErr.asm, dual_asm_polishVarCall, lr_het_hom_coll_fq.cov)

	sr_chunks_fq.combine(dual_asm_polishSNP).set { sr_chunks_dual_asm } // Combine each [sr.bam, sr.bam.bai] chunk with the created dual assembly

	sr_dual_asm_chunks_bam = mapPairedIllumina_polish(sr_chunks_dual_asm, false) // Map the splitted Illumina chunks to the dual assembly

	sr_dual_asm_chunks_bam.multiMap { it ->
		bam: it[0]
		bam_bai: it[1]
	}
	.set { sr_dual_asm_chunks_bam_split } // List of (bam, bam.bai) -> [[all bam], [all bam.bai]]

	// Merge the Illumina chunks to one single BAM file
	sr_dual_asm_chunks_bam_merged = mergePairedIlluminaBAM_polish(sr_dual_asm_chunks_bam_split.bam.collect(), sr_dual_asm_chunks_bam_split.bam_bai.collect())

	final_dual_asm = hapResAsmPolish_polishDualAsm_3(dual_asm_polishSNP, sr_dual_asm_chunks_bam_merged)
}
