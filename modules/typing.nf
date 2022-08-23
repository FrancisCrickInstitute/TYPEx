
process TYPE {
//	label  'max'
	 label 'xs' 

	input:
		tuple val(method)
		val subset
		val markers
		val stratify
		val nfDir
		val files
		

	output:	
		val params.outDir
	
	script:
	"""
		export BASE_DIR=$baseDir
		typing.R --wDir ${params.outDir} \
			--nfDir ${nfDir} \
			--method ${method}  \
			--subset ${subset} \
			--markers ${markers} \
			--run ${params.run} \
			--panel ${params.panel} \
			--celltypeReviewFile ${params.celltypeReviewFile} \
			--regFile ${params.regFile} \
			--tissAreaDir "${params.outDir}/tissue_seg/" \
			--cellAssignFile ${params.cellAssignFile} \
			--celltypeModelFile ${params.outDir}/review/${params.panel}_${params.major_markers}.RData \
			--stratify ${stratify} \
			--stratify_label "${params.stratify_label}" \
			--major_markers "${params.major_markers}"
			
	"""
}

process review_major_types {
	
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
		cell_typing_review_summary.R \
			--wDir ${params.outDir} \
			--cellReviewDir ${params.outDir}/review \
			--major_method ${params.major_method}  \
            --subset major --panel ${params.panel} \
			--subtype_markers ${params.subtype_markers} \
			--celltypeReviewFile ${params.celltypeReviewFile} \
			--stratify_label ${params.stratify_label}

	"""
}

