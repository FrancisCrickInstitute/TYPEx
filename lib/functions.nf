
def helpMessage() {

    out = 'TYPEx - Cell PhenoTYPing with MultiplEX imaging data'
    out += '''
    Usage:
        nextflow run -profile local --input <input_folder> --outdir <path> main.nf
    Mandatory arguments:
      -profile                      Configuration profile to use
      --input                       The folder with the required input files
      --outdir                      Directory to publish results
    '''.stripIndent()
	
    log.info(out)
}

def get_cell_files(imcyto, mccs) {

	if(imcyto) {
		println 'Processing files from deep-imcyto run'
		// tuple (imageID, file)
		if(mccs) {
			cellFilePattern="${params.input_dir}/consensus_cell_segmentation/**/cells.csv"
		} else {
			cellFilePattern="${params.input_dir}/simple_segmentation/**/cells.csv"
		}
		samplePattern=cellFilePattern
			.replaceAll('\\*', '.*')
			.replaceAll('.\\*.\\*', '(.*)')
			.replaceAll('\\/+', '.')
	
		cellFiles=Channel
			.fromPath("${cellFilePattern}", relative:false)
			.map{ file ->
						tuple(
							file.toString()
								.replaceAll(samplePattern, '$1')
								.replaceAll('\\/', '-'),
							file) }
	} else {
		// tuple ("-", file) when sample pattern not provided
		println 'Processing files independently from deep-imcyto'
		cellFiles=Channel
				.fromPath("${params.input_table}", relative:false)
				.map{ file -> tuple("-", file) }
	}
	return(cellFiles)
}

def get_imcyto_raw_masks() {

	rawMasksPattern="${params.input_dir}/imctools/**/full_stack/*.tiff"
	rawSamplePattern=rawMasksPattern.
		replaceAll('\\*', '.*').
		replaceAll('.\\*.\\*', '(.*)').
		replaceAll('\\/+', '.')
	println(rawSamplePattern)
	rawMasks=Channel
        .fromPath("${rawMasksPattern}", relative:false)
    	.map{ file -> tuple(file.toString()
							.replaceAll(rawSamplePattern, '$1')
							.replaceAll('\\/', '-'), file)
		}
	return(rawMasks)
}


def parse_json_file(overlayConfigFile) {

	def overlayConfig=file(overlayConfigFile)
	def jSlurp=new groovy.json.JsonSlurper()
	def ts_params=jSlurp.parse(overlayConfig)
	return(ts_params)
}


def get_tissue_masks_config(overlayConfigFile) {

	/*  Overlay cell objects with segmented/annotated tissue regions. 
	   Segmented tissue regions are either user-defined in tissue_segmentation.json 
	   and/or produced with the workflow TISSEG */
	   
	Channel.fromPath(overlayConfigFile) | 
			flatMap { parse_json_file(it) } |
			map {
			  entry -> tuple(
			  	(! file(entry.value.tissueDir).exists() ? "${params.output_dir}/tissue_seg" : entry.value.tissueDir),
				entry.value.maskRegEx, 
				entry.value.annotationName,
				entry.value.imgNameRegExIndex,
				entry.value.regionRegExIndex)
			}
}



