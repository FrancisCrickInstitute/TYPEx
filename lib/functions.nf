
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

def get_cell_files(imcyto, cellprofiler) {

	if(imcyto) {
		println 'Processing files from deep-imcyto run'
		// tuple (imageID, file)
		if(cellprofiler) {
			cellFilePattern="${params.input_dir}/consensus_cell_segmentation/**/cells.csv"
		} else {
			cellFilePattern="${params.input_dir}/simple_segmentation/**/cells.csv"
		}
		samplePattern=cellFilePattern
			.replaceAll('\\*', '.*')
			.replaceAll('.\\*.\\*\\/cells.csv', '(.*)/cells.csv')
			.replaceAll('.\\*.\\*', '.*')
			.replaceAll('\\/+', '.')
		samplePattern=".*" + samplePattern
	
		cellFiles=Channel
			.fromPath("${cellFilePattern}", relative:false, checkIfExists: true)
			.map{ file ->
						tuple(
							file.toString()
								.replaceAll(samplePattern, '$1')
								.replaceAll('\\/', '-'),
							file) }
			.ifEmpty(exit 1, "ERROR: Did not find deep-imcyto output in ${params.input_dir})
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
		replaceAll('.\\*.\\*\\/full_stack', '(.*)/full_stack').
		replaceAll('.\\*.\\*', '.*').
		replaceAll('\\/+', '.')
	rawSamplePattern=".*" + rawSamplePattern
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
			flatMap { parse_json_file(it).masks } |
			map {
			  entry -> tuple(
				(! file(entry.value.tissueDir).exists() ? "${params.output_dir}/tissue_seg" : entry.value.tissueDir),
				entry.value.maskRegEx, 
				entry.value.annotationName,
				entry.value.imgNameRegExIndex,
				entry.value.regionRegExIndex)
			}
}

def get_tissue_seg_markeres(overlayConfigFile) {

   entry =	Channel.fromPath(overlayConfigFile) |
      			flatMap { parse_json_file(it).markers } | collect(flat: true)
	println(entry)
	entry.view()
}


