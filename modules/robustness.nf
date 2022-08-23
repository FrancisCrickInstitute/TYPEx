process call_method {
    tag "major"
    label 'medium_mem'
    maxRetries 1

    input:
        tuple val(method)
        val nfDir
        val files

    output:
        val params.outDir

    //sbatch -n 1 --part=cpu --time=3-00:00:00 --mail-type=END --mail-user=mihaela.angelova@crick.ac.uk --mem-per-cpu=250G  \
    //                  --wrap="Rscript $SCRIPTDIR/typing.R --cohort $cohort --study $study --method $method --run $run --panel $panel --markers $marker --subset $subset "
    script:
    """
        export BASE_DIR=$baseDir
        typing.R --wDir ${params.outDir} --nfDir ${nfDir} \
            --cohort ${params.cohort} --study ${params.study} --method ${method}  \
            --subset major --markers ${params.major_markers} \
            --run ${params.run} --panel ${params.panel} \
            --celltypeReviewFile ${params.celltypeReviewFile} \
            --regFile ${params.regFile} \
			--tissAreaDir "${params.outDir}/tissue_seg" \
            --cellAssignFile ${params.cellAssignFile}
    """
}

process call_subsampled {
        tag "rob"
        label 'medium_mem'

        input:
                tuple val(iteration), val(method), val(markers)
				val nfDir
				val files

		output:
			val iteration

        script:
        """
                export BASE_DIR=$baseDir
                typing.R --wDir ${params.outDir}  --nfDir ${nfDir} \
                        --cohort ${params.cohort} --study ${params.study} --method ${method}  \
                        --subset sampled --markers $markers \
                        --run ${params.run} --panel ${params.panel} \
                        --celltypeReviewFile ${params.celltypeReviewFile} --regFile ${params.regFile} \
						--iter ${iteration} \
						--tissAreaDir "${params.outDir}/tissue_seg" \
						--cellAssignFile ${params.cellAssignFile}
        """
}


process match_clusters {
	tag "match"

	input:
		val inDir
		tuple val(method)
		val files
		val call_method
	output:
		val method
	
	script:
    """
        export BASE_DIR=$baseDir
		match_clusters.R --inDir ${inDir} --method ${method} \
			--run ${params.run} --panel ${params.panel} --reference_method ${method} 
    """
}

process plot_subsampled {
	
	input:
		val methods
	script:
    	"""
            export BASE_DIR=$baseDir
			# [[ ! -d "${params.outDir}/sampled/robustness/plots" ]] && mkdir "${params.outDir}/sampled/robustness/plots"
			matchedClusterStats.R --wDir "${params.outDir}/sampled/robustness" --outDir "${params.outDir}/sampled/robustness/plots"
	    """
}
