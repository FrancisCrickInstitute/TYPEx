
process exporter {

	label 'xs'
	
	publishDir path: "${params.outDir}/summary/", 
			   mode: params.publish_dir_mode, 
			   overwrite: true
	
    input:
        tuple val(method)
		val	subset
		tuple val (markers)
		val ref_markers
        val typing
		
	output:
		path ("*")

    script:
    """
        export BASE_DIR=$baseDir
		export PARAMS_CONF=${params.paramsConfig}
		export ANN_CONF=${params.annotationConfig}
		export COL_CONF=${params.colorConfig}
		
        cell_density_exporter.R \
			--subset ${subset} \
            --method $method --run ${params.release} \
            --panel ${params.panel} --markers ${markers} \
			--regFile ${params.sampleFile} --inDir ${params.outDir} \
			--tissAreaDir "${params.outDir}/tissue_seg/" \
			--ref_markers ${ref_markers}
    """
}

process plot_dr {
	
		label 'max'
        
        input:
                val method
                val nfDir
                val files

        output:
            val params.outDir

        script:
        """
                export BASE_DIR=$baseDir
				export PARAMS_CONF=${params.paramsConfig}
				export ANN_CONF=${params.annotationConfig}
				export COL_CONF=${params.colorConfig}
				
                typing.R --wDir ${params.outDir} --nfDir ${nfDir} \
						--method $method  \
                        --subset subtypes --markers ${params.subtype_markers} \
                        --run ${params.release} --panel ${params.panel} \
                        --regFile ${params.regFile} \
                        --celltypeModelFile ${params.outDir}/review/${params.panel}_${params.major_markers}.RData \
                        --tissAreaDir "${params.outDir}/tissue_seg/" \
						--stratify false
        """
}
