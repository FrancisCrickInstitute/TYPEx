

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
		csm_exporter.R --inDir ${csmDir} \
			--imcyto_run ${params.run} --panel ${params.panel} \
			--celltypeReviewFile ${params.celltypeReviewFile} \
			--markers ${params.major_markers} --nndist 5 \
			--tissAreaDir "${params.outDir}/tissue_seg/" 
	"""
}
