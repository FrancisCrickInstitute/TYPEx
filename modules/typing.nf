
process TYPE {

	label 'max'
	cache false

	input:
		tuple val(method)
		val subset
		tuple val (markers)
		tuple val (major_markers)
		val stratify
		path nfDir
		val files
		

	output:	
		val params.output_dir
	
	script:
	"""
		export BASE_DIR=$baseDir
		export PARAMS_CONF=${params.params_config}
		export ANN_CONF=${params.annotation_config}
		export COL_CONF=${params.color_config}
		
		typing.R --wDir ${params.output_dir} \
			--nfDir ${nfDir} \
			--method ${method}  \
			--subset ${subset} \
			--markers ${markers} \
			--run ${params.release} \
			--panel ${params.panel} \
			--regFile ${params.sample_file} \
			--tissAreaDir "${params.output_dir}/tissue_seg/" \
			--celltypeModelFile ${params.output_dir}/review/${params.panel}_${params.major_markers}.RData \
			--stratify ${stratify} \
			--mostFreqCellType "${params.exclude_cell_lineage}" \
			--major_markers "${major_markers}"

	"""
}

process build_strata_model {
	
	label 'medium_mem'
	cache false
	
	input:
		val files
		val ref
		val csm
		val tissue
	output:
		val params.output_dir
	script:
	"""
		export BASE_DIR=$baseDir
		export PARAMS_CONF=${params.params_config}
		export ANN_CONF=${params.annotation_config}
		export COL_CONF=${params.color_config}
		
		cell_typing_review_summary.R \
			--wDir ${params.output_dir} \
			--cellReviewDir ${params.output_dir}/review \
			--major_method ${params.major_method}  \
            --subset major \
			--panel ${params.panel} \
			--regFile ${params.sample_file} \
			--mostFreqCellType "${params.exclude_cell_lineage}" \
			--major_markers ${params.major_markers}

	"""
}

