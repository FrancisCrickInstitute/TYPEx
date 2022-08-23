
process exporter {

	label 'xs'
	
	publishDir path: "${params.outDir}/summary/", 
			   mode: params.publish_dir_mode, 
			   overwrite: true
	
    input:
        tuple val(method)
		val	subset
        val typing
		
	output:
		path ("*")

    script:
    """
        export BASE_DIR=$baseDir
        cell_density_exporter.R --subset ${subset} \
            --method $method --run ${params.run} \
            --panel ${params.panel} --markers ${params.subtype_markers} \
			--regFile ${params.regFile} --inDir ${params.outDir} \
			--cellAssignFile ${params.cellAssignFile} \
			--tissAreaDir "${params.outDir}/tissue_seg/"
    """
}


process plot_dr {
	
		label 'median_mem'  // 'max'
        
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
						--method $method  \
                        --subset subtypes --markers ${params.subtype_markers} \
                        --run ${params.run} --panel ${params.panel} \
                        --celltypeReviewFile ${params.celltypeReviewFile} --regFile ${params.regFile} \
                        --celltypeModelFile ${params.outDir}/review/${params.panel}_${params.major_markers}.RData \
                        --tissAreaDir "${params.outDir}/tissue_seg/" \
                        --cellAssignFile ${params.cellAssignFile} --stratify false
        """
}
