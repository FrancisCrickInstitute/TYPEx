
process TYPE {

	label 'max'

	input:
		tuple val(method)
		val subset
		tuple val (markers)
		tuple val (major_markers)
		val stratify
		path nfDir
		val files
		

	output:	
		val params.outDir
	
	script:
	"""
		export BASE_DIR=$baseDir
		export PARAMS_CONF=${params.paramsConfig}
		export ANN_CONF=${params.annotationConfig}
		export COL_CONF=${params.colorConfig}
		
		typing.R --wDir ${params.outDir} \
			--nfDir ${nfDir} \
			--method ${method}  \
			--subset ${subset} \
			--markers ${markers} \
			--run ${params.release} \
			--panel ${params.panel} \
			--regFile ${params.sampleFile} \
			--tissAreaDir "${params.outDir}/tissue_seg/" \
			--celltypeModelFile ${params.outDir}/review/${params.panel}_${params.major_markers}.RData \
			--stratify ${stratify} \
			--mostFreqCellType ${params.mostFreqCellType} \
			--major_markers "${major_markers}"

	"""
}

process build_strata_model {
	
	label 'xs'
	
	input:
		val files
		val ref
		val csm
		val tissue
	output:
		val params.outDir
	script:
	"""
		export BASE_DIR=$baseDir
		export PARAMS_CONF=${params.paramsConfig}
		export ANN_CONF=${params.annotationConfig}
		export COL_CONF=${params.colorConfig}
		
		cell_typing_review_summary.R \
			--wDir ${params.outDir} \
			--cellReviewDir ${params.outDir}/review \
			--major_method ${params.major_method}  \
            --subset major \
			--panel ${params.panel} \
			--regFile ${params.sampleFile} \
			--mostFreqCellType ${params.mostFreqCellType} \
			--major_markers ${params.major_markers}

	"""
}

