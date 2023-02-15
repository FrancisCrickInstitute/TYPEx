
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

