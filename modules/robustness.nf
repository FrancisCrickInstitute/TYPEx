
process call_subsampled {
        
		label 'medium_mem'

        input:
                tuple val(iteration), val(method), val(markers)
				val nfDir
				val files

		output:
			val iteration

        script:
        """
                export BASE_DIR=$baseDir
                typing.R --wDir ${params.outDir}  --nfDir ${nfDir} \
						--method ${method}  \
                        --subset sampled \
						--markers $markers \
                        --run ${params.run} --panel ${params.panel} \
                        --celltypeReviewFile ${params.celltypeReviewFile} --regFile ${params.regFile} \
						--iter ${iteration} \
						--tissAreaDir "${params.outDir}/tissue_seg" \
						--cellAssignFile ${params.cellAssignFile} \
						--stratify_label ${params.stratify_label}
        """
}


process match_clusters {
	tag "match"

	input:
		val inDir
		tuple val(method)
		val files
		val call_method
	output:
		val method
	
	script:
    """
        export BASE_DIR=$baseDir
		match_clusters.R --inDir ${inDir} --method ${method} \
			--run ${params.run} --panel ${params.panel} --reference_method ${method} 
    """
}

process plot_subsampled {
	
	input:
		val methods
	script:
    	"""
            export BASE_DIR=$baseDir
			matchedClusterStats.R --wDir "${params.outDir}/sampled/robustness" --outDir "${params.outDir}/sampled/robustness/plots"
	    """
}
