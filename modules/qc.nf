process qc_select_images {
	tag "match"

	input:
		tuple val (ref_markers)
	output:
		val method
	
	script:
    """
        export BASE_DIR=$baseDir
		export PARAMS_CONF=${params.params_config}
		export ANN_CONF=${params.annotation_config}
		export COL_CONF=${params.color_config}
		
        qc_select_images.R \
			--inDir "${params.output_dir}/summary/" \
			--outDir "${params.output_dir}/qc/" \
			--ref_markers ${ref_markers}

    """
}

process qc_create_single_channel_images {


		label 'x'

		maxRetries 1
        
		publishDir path: "${params.output_dir}/qc/", 
				   mode: params.publish_dir_mode, 
				   overwrite: true
				   
		input:
			val selection
			

		output:
			path ("*")

		script:
		"""
			## Using params.input_dir as an absolute path
			ImageJ-linux64 --ij2 --headless --run "${baseDir}/bin/raw_image_by_file.ijm" \
				'inDir=\"${params.image_dir}/\", posFile=\"${params.output_dir}/qc/overlay_examples.txt"'
		"""
}

process qc_overlay {
	
	input:
		val methods
	script:
    	"""
			# SET NEW DOCKER FILE
            export BASE_DIR=$baseDir
			export PARAMS_CONF=${params.params_config}
			export ANN_CONF=${params.annotation_config}
			export COL_CONF=${params.color_config}
			
			matchedClusterStats.R --wDir "${params.output_dir}/sampled/robustness" \
				--outDir "${params.output_dir}/sampled/robustness/plots"
	    """
}
