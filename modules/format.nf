
/*
 * STEP 1: based on output from imcyto -> format into intensity matrices
 */

process format_input {

	label 'xs'

	publishDir "${params.outDir}/nfData/", mode: params.publish_dir_mode, overwrite: true

	input:
		tuple val(imageID), path(inFile), val(feature)

	output:
		path("$feature/*")

	script:
		"""
			split_input_by_feature_per_file.pl ${inFile} ${imageID} ${feature} ${params.imcyto}
			
		"""
}

/*
 * STEP 2: merge formatted input per feature
 */
process collate_features {

	label 'medium_mem'
	
	publishDir "${params.outDir}/features/", mode: params.publish_dir_mode, overwrite: true
	
	input:
		tuple val(feature)
		val files

	output:
		path("${feature}/*")
	
    script:
    """
		
			export BASE_DIR=${baseDir}
			collate_by_feature.R --inDir "${params.outDir}/nfData/" \
				--panel ${params.panel} \
				--run ${params.run} \
				--feature ${feature}
				
    """
}
