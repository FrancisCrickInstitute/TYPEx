
process call_subsampled {


		label 'max'
	

        input:
                tuple val(iteration), val(method), val(markers)
				path nfDir
				val files

		output:
			val iteration

        script:
        """
                export BASE_DIR=$baseDir
				export PARAMS_CONF=${params.params_config}
				export ANN_CONF=${params.annotation_config}
				export COL_CONF=${params.color_config}
				
                typing.R --wDir ${params.output_dir}  \
						--nfDir ${nfDir} \
						--method ${method}  \
                        --subset sampled \
						--markers $markers \
                        --run ${params.release} \
						--panel ${params.panel} \
						--regFile ${params.regFile} \
						--iter ${iteration} \
						--tissAreaDir "${params.output_dir}/tissue_seg" \
					    --major_markers "${params.major_markers}"
        """
}

process match_clusters {
	tag "match"

	input:
		val inDir
		tuple val(method)
		tuple val (source_markers)
		tuple val (ref_markers)
		val files
		val call_method
	output:
		val method
	
	script:
    """
        export BASE_DIR=$baseDir
		export PARAMS_CONF=${params.params_config}
		export ANN_CONF=${params.annotation_config}
		export COL_CONF=${params.color_config}
		
		match_clusters.R --inDir ${inDir} \
			--method ${method} \
			--run ${params.release} \
			--panel ${params.panel} \
			--reference_method ${method} \
			--markers ${source_markers} \
			--reference_markers ${ref_markers}

    """
}

process plot_subsampled {
	
	input:
		val methods
	script:
    	"""
            export BASE_DIR=$baseDir
			export PARAMS_CONF=${params.params_config}
			export ANN_CONF=${params.annotation_config}
			export COL_CONF=${params.color_config}
			
			matchedClusterStats.R --wDir "${params.output_dir}/sampled/robustness" \
				--outDir "${params.output_dir}/sampled/robustness/plots"
	    """
}
