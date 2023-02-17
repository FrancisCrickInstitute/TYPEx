

process csm_submit {

	tag "CSM"
	label 'xs'

	input:
		tuple val(imageID), path(tifFile), path(inFile)
		val(csmDir)

	output:
		val (csmDir)

	script:
	"""
		echo ${inFile} ${csmDir}
		cell_specific_masks_per_image.py ${inFile} ${tifFile} ${imageID} ${csmDir} 
	"""
} 


process csm_export {
	
	tag "CSM export"
	label 'mem_medium'

	script:

	input:
		val (csmDir)
		val files
	output:
		val params.output_dir

	"""
		export BASE_DIR=$baseDir
		export PARAMS_CONF=${params.params_config}
		export ANN_CONF=${params.annotation_config}
		export COL_CONF=${params.color_config}
		
		csm_exporter.R --inDir ${csmDir} \
			--imcyto_run ${params.release} --panel ${params.panel} \
			--markers ${params.major_markers} --nndist 5 \
			--regFile ${params.sample_file} \
			--tissAreaDir "${params.output_dir}/tissue_seg/" 
	"""
}
