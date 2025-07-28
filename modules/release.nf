
process exporter {

	label 'xs'
	cache false
	
	publishDir path: "${params.output_dir}/summary/", 
			   mode: params.publish_dir_mode, 
			   overwrite: false
	
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
		export PARAMS_CONF=${params.params_config}
		export ANN_CONF=${params.annotation_config}
		export COL_CONF=${params.color_config}
		
        cell_density_exporter.R \
			--subset ${subset} \
            --method $method --run ${params.release} \
            --panel ${params.panel} --markers ${markers} \
			--regFile ${params.sample_file} --inDir ${params.output_dir} \
			--tissAreaDir "${params.output_dir}/tissue_seg/" \
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
            val params.output_dir

        script:
        """
                export BASE_DIR=$baseDir
				export PARAMS_CONF=${params.params_config}
				export ANN_CONF=${params.annotation_config}
				export COL_CONF=${params.color_config}
				
                typing.R --wDir ${params.output_dir} --nfDir ${nfDir} \
						--method $method  \
                        --subset subtypes --markers ${params.subtype_markers} \
                        --run ${params.release} --panel ${params.panel} \
                        --regFile ${params.regFile} \
                        --celltypeModelFile ${params.output_dir}/review/${params.panel}_${params.major_markers}.RData \
                        --tissAreaDir "${params.output_dir}/tissue_seg/" \
						--stratify false
        """
}
