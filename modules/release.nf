
process exporter {

        tag 'export'
        label 'xs'
        maxRetries 1

    input:
        tuple val(method)
        val typing
		val tisue_seg
		val nostrat

    script:
    """
        export BASE_DIR=$baseDir
        cell_density_exporter.R --subset subtypes \
            --method $method --run ${params.run} \
            --panel ${params.panel} --markers ${params.subtype_markers} \
			--regFile ${params.regFile} --inDir ${params.outDir} \
			--outDir ${params.outDir}/summary \
			--cohort ${params.cohort} \
			--cellAssignFile ${params.cellAssignFile} \
			--tissAreaDir "${params.outDir}/tissue_seg/"
    """
}
 

process exporter_major {

        tag 'export'
        label 'xs'
        maxRetries 1

    input:
        tuple val(method)
        val typing

    script:
    """
        export BASE_DIR=$baseDir
        cell_density_exporter.R --subset major \
            --method $method --run ${params.run} \
            --panel ${params.panel} --markers ${params.major_markers} \
			--regFile ${params.regFile} --inDir ${params.outDir} \
			--cohort ${params.cohort} \
			--outDir ${params.outDir}/summary \
			--cellAssignFile ${params.cellAssignFile} \
			--tissAreaDir "${params.outDir}/tissue_seg/"
    """
}


process plot_dr {

        tag "sub"
        label  'median_mem' //'max'
        maxRetries 1

        input:
                val method
                val nfDir
                val files

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
                        --cellAssignFile ${params.cellAssignFile} --stratify false
        """
}
