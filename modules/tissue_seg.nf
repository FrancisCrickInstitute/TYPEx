
process preprocess_panck {

	maxRetries 1

	publishDir path: "${params.output_dir}/composites/", 
			   mode: params.publish_dir_mode, 
			   overwrite: false

	input:
		tuple val(tumour), val(immune), val(stroma), val(auxStroma), val(dna)

	output:
		path ("*")

	script:
	"""
		## Using params.input_dir as an absolute path
		ImageJ-linux64 --ij2 --headless --run "${baseDir}/bin/panck_filters.ijm" \
			'rootDir=\"${params.image_dir}/\", tumour=\"${tumour}\"'
	
	"""
}


process create_composites {

    label 'medium_mem'
	
	publishDir path:"${params.output_dir}/composites/", 
			   mode: params.publish_dir_mode, 
			   overwrite: false

	input:
		val files
		tuple val(tumour), val(immune), val(stroma), val(auxStroma), val(dna)

	output:
		path ("*")
		

    script:
    """
		# compositeDir has the output directory from the previous
		ImageJ-linux64 --ij2 --headless \
			--run "${baseDir}/bin/create_composites.ijm" \
				'rootDir="${params.image_dir}/", panel=\"${params.panel}\", compositeDir=\"${params.output_dir}/composites/\", tumour=\"${tumour}\", immune=\"${immune}\", stroma=\"${stroma}\", auxStroma=\"${auxStroma}\", dna=\"${dna}\"'

    """

}

process run_classifier {

	label "s"
	
	// publishDir path: "${params.output_dir}/composites/probs/", 
			//	  mode: params.publish_dir_mode, 
				//  overwrite: true

	input:
		val files
	output:
		val params.output_dir

	script:
	
	"""
		files=` find ${params.output_dir}/composites/ -name "*tiff" `
		echo $files
		if [ -n "` find ${params.output_dir}/composites/ -name "*tiff" `" ]; then
			/ilastik-release/run_ilastik.sh --headless \
				--project=${params.tissue_seg_model} \
				${params.output_dir}/composites/*tiff
		else
			echo 'No composite images for tissue segmentation created'
		fi
	
	"""
}

process process_probs {

	
	maxRetries 1

	publishDir path: "${params.output_dir}/tissue_seg/", 
			   mode: params.publish_dir_mode,
			   overwrite: true

	input:
        val files
    output:
        path("*")

	script:
    """
	
		ImageJ-linux64 --ij2 --headless \
			--run $baseDir/bin/tissueseg_postprocess.ijm \
			'wDir=\"${params.output_dir}/\"'
		
    """	
}

process ts_exporter {
	
	maxRetries 1

	input:
		val files

	output:
		val params.output_dir
		
	script:
    """
	
		export BASE_DIR=$baseDir
		export PARAMS_CONF=${params.params_config}
		export ANN_CONF=${params.annotation_config}
		export COL_CONF=${params.color_config}
		
		summarise_tissue_seg.R \
			--tissAreaDir "${params.output_dir}/tissue_seg" \
			--panel ${params.panel}
			
    """
}

process mask_overlay {

	maxRetries 1
    tag "post-process"
	
	publishDir "${params.output_dir}/tissue_seg/",
				mode: params.publish_dir_mode,
				overwrite: false

	input:
		tuple path(maskDir), val(maskRegEx), val(regionType), val(imageRegEx), val(regionRegEx)
		val cellObjFile
		val files
		val cells

	output:
		path("*")

	script:
    """
			
		cell_objects_regional_info.py \
			--maskDir ${maskDir} \
			--tissueAreaDir "${params.output_dir}/tissue_seg" \
			--panel ${params.panel} \
			--cellObjFile ${cellObjFile} \
			--maskRegEx "${maskRegEx}" \
			--regionType "${regionType}" \
			--regionRegEx "${regionRegEx}" \
			--imageRegEx "${imageRegEx}"

    """

}


