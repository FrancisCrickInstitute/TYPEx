

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
		val params.outDir

	"""
		export BASE_DIR=$baseDir
		export PARAMS_CONF=${params.paramsConfig}
		export ANN_CONF=${params.annotationConfig}
		export COL_CONF=${params.colorConfig}
		
		csm_exporter.R --inDir ${csmDir} \
			--imcyto_run ${params.release} --panel ${params.panel} \
			--markers ${params.major_markers} --nndist 5 \
			--regFile ${params.sampleFile} \
			--tissAreaDir "${params.outDir}/tissue_seg/" 
	"""
}
