
process preprocess_panck {

	maxRetries 1

	publishDir path: "${params.outDir}/composites/", 
			   mode: params.publish_dir_mode, 
			   overwrite: true

	output:
		path ("*")

	script:
	"""
		## Using params.inDir as an absolute path
		ImageJ-linux64 --ij2 --headless --run "${baseDir}/bin/panck_filters.ijm" \
			'rootDir=\"${params.imageDir}/\"'
	
	"""
}


process create_composites {

    label 'medium_mem'
	
	publishDir path:"${params.outDir}/composites/", 
			   mode: params.publish_dir_mode, 
			   overwrite: true

	input:
		val files

	output:
		path ("*")
		

    script:
    """
		# compositeDir has the output directory from the previous
		ImageJ-linux64 --ij2 --headless \
			--run "${baseDir}/bin/create_composites.ijm" \
				'rootDir="${params.imageDir}/", panel=\"${params.panel}\", mcd=true, compositeDir=\"${params.outDir}/composites/\"'

    """

}

process run_classifier {
	
	// publishDir path: "${params.outDir}/composites/probs/", 
			//	  mode: params.publish_dir_mode, 
				//  overwrite: true

	input:
		val files
	output:
		val params.outDir

	script:
	
	"""
		files=` find ${params.outDir}/composites/ -name *tiff `
		echo $files
		if [ -n "` find ${params.outDir}/composites/ -name *tiff`" ]; then
			/ilastik-release/run_ilastik.sh --headless \
				--project=${params.tissue_seg_model} \
				${params.outDir}/composites/*tiff
		else
			echo 'No composite images for tissue segmentation created'
		fi
	
	"""
}

process process_probs {

	
	maxRetries 1

	publishDir path: "${params.outDir}/tissue_seg/", 
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
			'wDir=\"${params.outDir}/\"'
		
    """	
}

process ts_exporter {
	
	maxRetries 1

	input:
		val files

	output:
		val params.outDir
		
	script:
    """
	
		export BASE_DIR=$baseDir
		export PARAMS_CONF=${params.paramsConfig}
		export ANN_CONF=${params.annotationConfig}
		export COL_CONF=${params.colorConfig}
		
		summarise_tissue_seg.R \
			--tissAreaDir "${params.outDir}/tissue_seg" \
			--panel ${params.panel}
			
    """
}

process mask_overlay {

	maxRetries 1
    tag "post-process"
	
	publishDir "${params.outDir}/tissue_seg/",
				mode: params.publish_dir_mode,
				overwrite: true

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
			--tissueAreaDir "${params.outDir}/tissue_seg" \
			--panel ${params.panel} \
			--cellObjFile ${cellObjFile} \
			--maskRegEx "${maskRegEx}" \
			--regionType "${regionType}" \
			--regionRegEx "${regionRegEx}" \
			--imageRegEx "${imageRegEx}"

    """

}


