process mapLRtoAsm {

	label 'medium_node'

	input:
		path lr_fq
		tuple path(asm_fa), path(asm_fai), path(asm_bed)

	output:
		tuple path('lr.asm.bam'), path('lr.asm.bam.bai'), emit: bam
		path('lr.asm.cov'), emit: cov

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	minimap2=\${MINIMAP2:-${params.tools.minimap2.bin}}

	TASK_MEM=\$(echo -e \"${task.memory}\" | cut -d \" \" -f1)
	MEM_PER_THREADS_SORT=\$(bc -l <<< \"(\${TASK_MEM} / ${task.cpus}) * ${params.tools.samtools.sort.mem_safety_ratio} * 1000\" | awk '{printf(\"%.0f\", \$0)}');

	\${minimap2} -t ${task.cpus} ${params.tools.minimap2.param.ont_reads} -Y -I \$(du -L -BG ${asm_fa} | cut -f1) ${asm_fa} ${lr_fq} > lr.asm.sam;
	\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M lr.asm.sam > lr.asm.bam;
	if [ ! -s lr.asm.bam ]; then echo \"File lr.asm.bam does not exist or is empty\" 1>&2; exit 1; fi;
	\${samtools} quickcheck lr.asm.bam;
	if [ ! \$? -eq 0 ]; then echo \"File lr.asm.bam is malformed\" 1>&2; exit 1; fi;
	\${samtools} index -@ ${task.cpus} lr.asm.bam; rm -rf lr.asm.sam;

	\${samtools} depth -@ ${task.cpus} -aa -J -Q ${params.pipeline.min_mapq_strict} -G ${params.pipeline.samtools.depth.filter_sec} lr.asm.bam | \
	awk 'BEGIN {FS=\"\\t\"; OFS=\"\\t\"; CONTIG=\"\"; SUM=0; COUNT=0; SUM_ALL=0; COUNT_ALL=0} \
	{if (\$1!=CONTIG) {if (CONTIG!=\"\") {MEAN=0; if (COUNT>=1) {MEAN=SUM/COUNT}; print CONTIG, COUNT, MEAN}; CONTIG=\$1; SUM=0; COUNT=0}; SUM+=\$3; COUNT+=1; SUM_ALL+=\$3; COUNT_ALL+=1} \
	END {if (CONTIG!=\"\") {MEAN=0; if (COUNT>=1) {MEAN=SUM/COUNT}; print CONTIG, COUNT,MEAN}; MEAN_ALL=0; if (COUNT_ALL>=1) {MEAN_ALL=SUM_ALL/COUNT_ALL}; print \"all\", COUNT_ALL, MEAN_ALL}' > lr.asm.cov;
	"""
}

process var_CallFilterPhase {

	label 'medium_node'

	input:

		tuple path(asm_fa), path(asm_fai), path(asm_bed)
		tuple path('lr.bam'), path('lr.bam.bai')
		path var_regions_bed // Optional, must be empty list ([]) if not used. Workaround till Nextflow implement this.
		val skipHaplotypeBAM // true or false

	output:

		tuple path('margin/MARGIN_PHASED.haplotagged.bam'), path('margin/MARGIN_PHASED.haplotagged.bam.bai'), optional: true, emit: hap_bam
		tuple path('margin/MARGIN_PHASED.phased.vcf.gz'), path('margin/MARGIN_PHASED.phased.vcf.gz.tbi'), emit: hap_vcf
		path 'margin/MARGIN_PHASED.phaseset.bed', emit: ps_bed

	shell '/bin/bash', '-euo', 'pipefail'

	script:

		def var_regions_bed_bash = var_regions_bed != [] ? "${var_regions_bed}" : ''
		def skipHaplotypeBAM_opt = skipHaplotypeBAM == true ? "--skipHaplotypeBAM" : ''

		"""
		bcftools=\${BCFTOOLS:-${params.tools.bcftools.bin}}

		mkdir -p margin

		run_pepper_margin_deepvariant call_variant -b lr.bam -f ${asm_fa} -o pepper -t ${task.cpus} -s Sample \
		--hifi --pepper_include_supplementary --pepper_min_mapq ${params.pipeline.min_mapq.strict} --dv_min_mapping_quality ${params.pipeline.min_mapq.strict}

		if [ -s ${var_regions_bed_bash} ]
		then

			\${bcftools} view -i \"(type=='snp') && (FILTER=='PASS') && ((GT=='AR') & (GQ>=${params.pipeline.variant_filter.min_gq}) & \
			(AD>=${params.pipeline.variant_filter.min_ad}) & (VAF>=${params.pipeline.variant_filter.min_vaf}) & (VAF<=${params.pipeline.variant_filter.max_vaf}))\" \
			-R ${var_regions_bed_bash} -Oz pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz | \${bcftools} sort -Oz -o pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.filtered.vcf.gz
		else

			\${bcftools} view -i \"(type=='snp') && (FILTER=='PASS') && ((GT=='AR') & (GQ>=${params.pipeline.variant_filter.min_gq}) & \
			(AD>=${params.pipeline.variant_filter.min_ad}) & (VAF>=${params.pipeline.variant_filter.min_vaf}) & (VAF<=${params.pipeline.variant_filter.max_vaf}))\" \
			-Oz pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz | \${bcftools} sort -Oz -o pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.filtered.vcf.gz
		fi

		tabix -p vcf pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.filtered.vcf.gz

		margin phase lr.bam ${asm_fa} pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.filtered.vcf.gz ${projectDir}/${params.pmdv.r08.margin.config} \
		-t ${task.cpus} -o margin/MARGIN_PHASED ${skipHaplotypeBAM_opt}

		\${bcftools} view -Oz -o margin/MARGIN_PHASED.phased.vcf.gz margin/MARGIN_PHASED.phased.vcf
		tabix -p vcf margin/MARGIN_PHASED.phased.vcf.gz;

		rm -rf margin/MARGIN_PHASED.phased.vcf

		if [ -s margin/MARGIN_PHASED.haplotagged.bam ]
		then
			samtools index -@ ${task.cpus} margin/MARGIN_PHASED.haplotagged.bam
		fi
		"""
}

process varCall_hifi {

	label 'medium_node'

	input:

		tuple path(asm_fa), path(asm_fai), path(asm_bed)
		tuple path('lr.bam'), path('lr.bam.bai')

	output:

		tuple path('PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz'), path('PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz.tbi')

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	run_pepper_margin_deepvariant call_variant -b lr.bam -f ${asm_fa} -o . -t ${task.cpus} -s Sample \
	--hifi --pepper_include_supplementary --pepper_min_mapq ${params.pipeline.min_mapq.strict} --dv_min_mapping_quality ${params.pipeline.min_mapq.strict}
	"""
}

process varCall_ontR9_trainedModels {

	label 'medium_node'

	input:
		tuple path('lr.bam'), path('lr.bam.bai')
		tuple path(asm_fa), path(asm_fai), path(asm_bed)

	output:
		tuple path('PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz'), path('PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz.tbi')

	shell '/bin/bash', '-euo', 'pipefail'

	script:

		def model_snp = ("${params.pmdv_r07_model_r9_guppy5_sup_snp}" == "") ? "${projectDir}/${params.pmdv.r07.models.snp}" : "${params.pmdv_r07_model_r9_guppy5_sup_snp}"
		def model_hp = ("${params.pmdv_r07_model_r9_guppy5_sup_hp}" == "") ? "${projectDir}/${params.pmdv.r07.models.hp}" : "${params.pmdv_r07_model_r9_guppy5_sup_hp}"
		def model_dv = ("${params.pmdv_r07_model_r9_guppy5_sup_dv}" == "") ? "${projectDir}/${params.pmdv.r07.models.dv}" : "${params.pmdv_r07_model_r9_guppy5_sup_dv}"

		"""
		run_pepper_margin_deepvariant call_variant -b lr.bam -f ${asm_fa} -o . -t ${task.cpus} -s Sample --ont_r9_guppy5_sup \
		--pepper_model ${model_snp} --pepper_hp_model ${model_hp} --dv_model ${model_dv} \
		--pepper_min_mapq ${params.pipeline.min_mapq.strict} --dv_min_mapping_quality ${params.pipeline.min_mapq.strict} \
		--pepper_include_supplementary --no_pepper_quantized
		"""
}

process varPhase {

	label 'medium_node'

	input:

		tuple path(asm_fa), path(asm_fai), path(asm_bed)
		tuple path('lr.bam'), path('lr.bam.bai')
		tuple path('lr.asm.vcf.gz'), path('lr.asm.vcf.gz.tbi')
		val skipHaplotypeBAM // true or false (whether phased BAM is output)

	output:

		tuple path('margin/MARGIN_PHASED.haplotagged.bam'), path('margin/MARGIN_PHASED.haplotagged.bam.bai'), optional: true, emit: hap_bam
		tuple path('margin/MARGIN_PHASED.phased.vcf.gz'), path('margin/MARGIN_PHASED.phased.vcf.gz.tbi'), emit: hap_vcf
		path 'margin/MARGIN_PHASED.phaseset.bed', emit: ps_bed

	shell '/bin/bash', '-euo', 'pipefail'

	script:

		def skipHaplotypeBAM_opt = (skipHaplotypeBAM == true) ? "--skipHaplotypeBAM" : ''

		"""
		bcftools=\${BCFTOOLS:-${params.tools.bcftools.bin}}

		mkdir -p margin

		margin phase lr.bam ${asm_fa} lr.asm.vcf.gz ${projectDir}/${params.pmdv.r08.margin.config} -t ${task.cpus} -o margin/MARGIN_PHASED ${skipHaplotypeBAM_opt}

		\${bcftools} view -Oz -o margin/MARGIN_PHASED.phased.vcf.gz margin/MARGIN_PHASED.phased.vcf
		tabix -p vcf margin/MARGIN_PHASED.phased.vcf.gz;

		rm -rf margin/MARGIN_PHASED.phased.vcf

		if [ -s margin/MARGIN_PHASED.haplotagged.bam ]
		then
			samtools index -@ ${task.cpus} margin/MARGIN_PHASED.haplotagged.bam
		fi
		"""
}


/*
Map Illumina reads to assembly, sort and index output BAM file.
*/
process mapPairedIllumina {

	label 'medium_node'

	input:
		tuple path(sr_fq), path(asm_fa), path(asm_fai), path(asm_bed)
		val keepSecondary // true or false (whether we want the alignment to keep all secondary alignments)

	output:
		tuple path("sr.${task.process}.bam"), path("sr.${task.process}.bam.bai")

	shell '/bin/bash', '-euo', 'pipefail'

	script:

		def keepAllSecondary_opt = (keepSecondary == true) ? "--secondary=yes" : ''

		"""
		samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
		minimap2=\${MINIMAP2:-${params.tools.minimap2.bin}}

		TASK_MEM=\$(echo -e \"${task.memory}\" | cut -d \" \" -f1) # Get total RAM allocated to job in GB
		MEM_PER_THREADS_SORT=\$(bc -l <<< \"(\${TASK_MEM} / ${task.cpus}) * ${params.tools.samtools.sort.mem_safety_ratio} * 1000\" | awk '{printf(\"%.0f\", \$0)}') # Compute RAM available per core in MB

		\${minimap2} -t ${task.cpus} -ax sr -Y -I \$(du -L -BG ${asm_fa} | cut -f1) ${keepAllSecondary_opt} ${asm_fa} ${sr_fq} > sr.${task.process}.sam # Map long reads back to assembly
		\${samtools} sort -@ ${task.cpus} -m \${MEM_PER_THREADS_SORT}M -o sr.${task.process}.bam sr.${task.process}.sam # Create BAM file from sorted SAM
		if [ ! -s sr.${task.process}.bam ]; then echo \"File sr.${task.process}.bam does not exist or is empty\" 1>&2; exit 1; fi
		\${samtools} quickcheck sr.${task.process}.bam # Run file truncation check on resulting BAM file
		if [ ! \$? -eq 0 ]; then echo \"File sr.${task.process}.bam is malformed\" 1>&2; exit 1; fi
		\${samtools} index -@ ${task.cpus} sr.${task.process}.bam
		rm -rf sr.${task.process}.sam
		"""
}

process mergePairedIlluminaBAM {

	label 'medium_node'

	input:
		path("?????????.bam")
		path("?????????.bam.bai")

	output:
		tuple path("sr.${task.process}.bam"), path("sr.${task.process}.bam.bai")

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}

	\${samtools} merge -@ ${task.cpus} -o sr.${task.process}.bam *.bam # Merge the BAM files
	if [ ! -s sr.${task.process}.bam ]; then echo \"File sr.${task.process}.bam does not exist or is empty\" 1>&2; exit 1; fi
	\${samtools} quickcheck sr.${task.process}.bam # Run file truncation check on resulting BAM file
	if [ ! \$? -eq 0 ]; then echo \"File sr.${task.process}.bam is malformed\" 1>&2; exit 1; fi
	\${samtools} index -@ ${task.cpus} sr.${task.process}.bam
	"""
}

