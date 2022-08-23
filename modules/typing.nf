
// ml load Anaconda3
// #ml load R/3.6.0-foss-2019a
// #export NUMTHREADS=1

// #SCRIPTDIR=~/lung/scripts/typing
// #export PATH="/camp/home/angelom/.virtualenvs/r-tensorflow/bin":$PATH
// #export PYTHONPATH="/camp/home/angelom/.virtualenvs/r-tensorflow/bin":"/camp/home/angelom/.virtualenvs/r-tensorflow/lib/python3.7/":"/camp/home/angelom/.virtualenvs/r-tensorflow/lib/python3.7/site-packages/":$PYTHONPATH
// #export TENSORFLOW_PYTHON="/camp/home/angelom/.virtualenvs/r-tensorflow/bin/python3.7"


process call_major {
	tag "major"
	label 'xs' //'medium_mem' //'max'
	maxRetries 1
	
	input:
		tuple val(method)
		val nfDir
		val files

	output:	
		val params.outDir
	
	//sbatch -n 1 --part=cpu --time=3-00:00:00 --mail-type=END --mail-user=mihaela.angelova@crick.ac.uk --mem-per-cpu=250G  \
	//					--wrap="Rscript $SCRIPTDIR/typing.R --cohort $cohort --study $study --method $method --run $run --panel $panel --markers $marker --subset $subset "
	script:
	"""
		export BASE_DIR=$baseDir
		typing.R --wDir ${params.outDir} --nfDir ${nfDir} \
			--cohort ${params.cohort} --study ${params.study} --method ${method}  \
			--subset major --markers ${params.major_markers} \
			--run ${params.run} --panel ${params.panel} \
			--celltypeReviewFile ${params.celltypeReviewFile} \
			--regFile ${params.regFile} \
			--tissAreaDir "${params.outDir}/tissue_seg/" \
			--cellAssignFile ${params.cellAssignFile}
	"""
} 

process call_major_ref {
    tag "major"
    label 'xs' // 'medium_mem' //'max'
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
            --subset major --markers ${params.major_markers_ref} \
            --run ${params.run} --panel ${params.panel} \
            --celltypeReviewFile ${params.celltypeReviewFile} \
            --regFile ${params.regFile} \
			--tissAreaDir "${params.outDir}/tissue_seg/" \
            --cellAssignFile ${params.cellAssignFile}
    """
}


process review_major_types {

	tag "review"
    label 'xs'
    maxRetries 1
	
	input:
		val files
		val ref
		val csm
		val tissue
	output:
		val params.outDir
	script:
	"""
		export BASE_DIR=$baseDir
		cell_typing_review_summary.R --wDir ${params.outDir} --cellReviewDir ${params.outDir}/review \
            --cohort ${params.cohort} --study ${params.study} \
			--major_method ${params.major_method}  \
            --subset major --run ${params.run} --panel ${params.panel} \
			--subtype_markers all \
			--celltypeReviewFile  ${params.celltypeReviewFile}

	"""
}

process call_subtypes {

        tag "sub"
        label 'xs' // 'medium_mem'
        maxRetries 1

        input:
                tuple val(method)
				val nfDir
                val files
				val tissue_seg

		output:
			val params.outDir

        script:
        """
                export BASE_DIR=$baseDir
                typing.R --wDir ${params.outDir} --nfDir ${nfDir} \
                        --cohort ${params.cohort} --study ${params.study} --method $method  \
                        --subset subtypes --markers ${params.subtype_markers} \
                        --run ${params.run} --panel ${params.panel} \
                        --celltypeReviewFile ${params.celltypeReviewFile} --regFile ${params.regFile} \
						--celltypeModelFile ${params.outDir}/review/${params.cohort}_${params.run}_${params.panel}_${params.major_markers}.RData \
						--tissAreaDir "${params.outDir}/tissue_seg/" \
						--cellAssignFile ${params.cellAssignFile}
        """
}

process call_subtypes_nostrat {

        tag "sub"
        label 'xs' // 'medium_mem'
        maxRetries 1

        input:
                tuple val(method)
                val nfDir
                val files
				val tissue_seg

		output:
			val params.outDir

        script:
        """
                export BASE_DIR=$baseDir
                typing.R --wDir ${params.outDir} --nfDir ${nfDir} \
                        --cohort ${params.cohort} --study ${params.study} --method $method  \
                        --subset subtypes --markers ${params.subtype_markers} \
                        --run ${params.run} --panel ${params.panel} \
                        --celltypeReviewFile ${params.celltypeReviewFile} --regFile ${params.regFile} \
						--tissAreaDir "${params.outDir}/tissue_seg/" \
                        --celltypeModelFile ${params.outDir}/review/${params.cohort}_${params.run}_${params.panel}_${params.major_markers}.RData \
						--stratify false \
						 --cellAssignFile ${params.cellAssignFile}
        """
}

