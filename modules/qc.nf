process qc_select_images {

	label 'x'

	input:
		tuple val (ref_markers)
		val subset
		tuple val (method)
		val typing

	output:
		val method
	
	script:
    """
        export BASE_DIR=$baseDir
		export PARAMS_CONF=${params.params_config}
		export ANN_CONF=${params.annotation_config}
		export COL_CONF=${params.color_config}
		
        qc_select_images.R \
			--inDir "${params.output_dir}/summary/${subset}_${ref_markers}_${method}/" \
			--outDir "${params.output_dir}/summary/${subset}_${ref_markers}_${method}/qc/" \
			--ref_markers ${params.major_markers}

    """
}

process qc_create_single_channel_images {


		label 'x'

		maxRetries 1
        
		input:
			tuple val (ref_markers)
        	val subset
        	tuple val (method)
			val selection

		output:
			val method

		script:
		"""
			[[ ! -d ${subset}_${ref_markers}_${method}/qc ]]
				mkdir -p ${subset}_${ref_markers}_${method}/qc
			## Using params.input_dir as an absolute path
			ImageJ-linux64 --ij2 --headless --run "${baseDir}/bin/raw_image_by_file.ijm" \
				'inDir=\"${params.image_dir}/\", posFile=\"${params.output_dir}/summary/${subset}_${ref_markers}_${method}/qc/overlay_examples.txt", outDir=\"${params.output_dir}/summary/${subset}_${ref_markers}_${method}/qc\"'
		"""
}

process qc_overlay {
	
	input:
		tuple val (ref_markers)
		val subset
		tuple val (method)
		val mccs
		val images

	script:
    	"""
            export BASE_DIR=$baseDir
			export PARAMS_CONF=${params.params_config}
			export ANN_CONF=${params.annotation_config}
			export COL_CONF=${params.color_config}
			
			plot_overlays.R \
				--rawDir ${params.output_dir}/summary/${subset}_${ref_markers}_${method}/qc/ \
				--maskDir ${params.input_dir} --mccs ${mccs}	 \
				--inDir ${params.output_dir}/summary/${subset}_${ref_markers}_${method}/ \
				--outDir "${params.output_dir}/summary/${subset}_${ref_markers}_${method}/" \
				--posFile ${params.output_dir}/summary/${subset}_${ref_markers}_${method}/qc/overlay_examples.txt \
				--panel ${params.panel} \
				--run ${params.release}
	    """
}

process qc_intensity_heatmap {
	
	input:
		tuple val (ref_markers)
		val subset
		tuple val (method)
		val typing
		

	script:
    	"""
            export BASE_DIR=$baseDir
			export PARAMS_CONF=${params.params_config}
			export ANN_CONF=${params.annotation_config}
			export COL_CONF=${params.color_config}
			
			plot_intensity_heatmap.R --inDir ${params.output_dir}/ \
				--panel ${params.panel} \
				--subset ${subset} \
				--method ${method} \
				--markers ${ref_markers}

	    """
}
