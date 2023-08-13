
/*
 * STEP 1: based on output from imcyto -> format into intensity matrices
 */

process format_input {

	label 's'
//	executor "local"

	publishDir "${params.output_dir}/nfData/", mode: params.publish_dir_mode, overwrite: true

	input:
		tuple val(imageID), path(inFile), val(feature)

	output:
		path("$feature/*")

	script:
		"""
			split_input_by_feature_per_file.pl ${inFile} ${imageID} ${feature} ${params.deep_imcyto}
			
		"""
}

/*
 * STEP 2: merge formatted input per feature
 */
process collate_features {

	label 's'
	
	publishDir "${params.output_dir}/features/", mode: params.publish_dir_mode, overwrite: true
	
	input:
		tuple val(feature)
		val files

	output:
		path("${feature}/*")
	
    script:
    """
		
			export BASE_DIR=${baseDir}
			export PARAMS_CONF=${params.params_config}
			export ANN_CONF=${params.annotation_config}
			export COL_CONF=${params.color_config}
			
			collate_by_feature.R --inDir "${params.output_dir}/nfData/" \
				--panel ${params.panel} \
				--run ${params.release} \
				--feature ${feature}
				
    """
}
