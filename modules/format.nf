
/*
 * STEP 1: based on output from imctools -> format into intenstity matrices
 */

process format_input {

	tag "Formating imctools output as input for typing"
	label 'xs'
	executor 'local'
	
	input:
		tuple val(imageID), path(inFile)
		path nfOut

	output:
		path("${nfOut}")
		// export PATH="/camp/lab/swantonc/working/angelom/tools/miniconda3/bin":$PATH

	script:
		"""
			echo ${inFile} ${nfOut} ${imageID}
			# [[ ! -d "${nfOut}" ]] && mkdir -p ${nfOut}
			split_files_by_feature_per_file.pl ${inFile} ${imageID} ${nfOut} ${params.run}
				# "${params.inDir}/${params.panel}/${params.run}"
		"""
}

process collate_features {

    tag "Collate"
	label 'medium_mem'

	input:
		path inDir
		path outDir
		tuple val(feature)
		val files

	output:
		val(outDir)
	
    script: // Script in typex/bin/format/
    """
			export BASE_DIR=${baseDir}

			# echo $inDir $outDir
			[[ ! -d "${outDir}" ]] && mkdir -p ${outDir}
			collate_by_feature.R --inDir "${inDir}/" --outDir "${outDir}/${feature}" --panel ${params.panel} --run ${params.run} --feature ${feature}
                #                --inDir ${params.inDir}/nfrun --outDir ${params.outDir}/features \
    """
}
