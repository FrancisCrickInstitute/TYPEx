
process preprocess_panck {

	maxRetries 1
    tag "composites"

	input:
		path inDir
	output:
		val params.outDir

	script:
	"""
		# rootDir="/camp/project/proj-tracerx-lung/tctProjects/rubicon/" + study + "/" + cohort + "/imc/outputs/nextflow/" + panel + "/" + runPattern + "/";
		# outDir="/camp/home/angelom/labwd/analyses/imc//tissue_classification/composites/reanalysed_" + runPattern + "/";
		ImageJ-linux64 --ij2 --headless --run $baseDir/bin/panck_filters.ijm 'rootDir=\"${params.inDir}/${params.panel}/${params.run}/\", outDir=\"${params.outDir}/composites/\", panel=\"${params.panel}\", runPattern="\${params.run}\", mcd=true '
	
	"""
}


process create_composites {

    maxRetries 1
    tag "composites"
//	label 'medium_mem'

	input:
		path inDir
		val files

	output:
		val params.outDir

    script:
    """
		ImageJ-linux64 --ij2 --headless --run $baseDir/bin/create_composites.ijm 'rootDir=\"${params.inDir}/${params.panel}/${params.run}/\", outDir=\"${params.outDir}/composites/\", panel=\"${params.panel}\", mcd=true '

    """
}

process run_classifier {

	maxRetries 1
    tag "post-process"
//	label 's'

	input:
		val files
	output:
		val params.outDir

	script:
	"""
	   /ilastik-release/run_ilastik.sh --headless --project=$baseDir/model/Classifier_final_ts_simpler_3.ilp ${params.outDir}/composites/P1_TMA_REC*roi_5*tiff
	"""
}

process process_probabilities {

	
	maxRetries 1
    tag "post-process"

	input:
        val files
    output:
        val params.outDir

	script:
    """
		ImageJ-linux64 --ij2 --headless --run $baseDir/bin/tissueseg_postprocess.ijm 'wDir=\"${params.outDir}/\", outDir=\"${params.outDir}/tissue_seg/\"'
    """	
}

process overlay {

	maxRetries 1
    tag "post-process"

	input: 
		val files

	output:
		val params.outDir

	script:
    """
		cell_objects_regional_info.py --analysisDir ${params.outDir} --run ${params.run} --panel ${params.panel}
    """

}

process export_tissue_seg {
	
	maxRetries 1
    tag "post-process"

	input:
		val files

	output:
		val params.outDir
	script:
    """
		export BASE_DIR=$baseDir
		summarise_tissue_seg.R --tissAreaDir "${params.outDir}/tissue_seg" \
			 --panel ${params.panel} \
			--cohort ${params.cohort} --study ${params.study}
    """
}


