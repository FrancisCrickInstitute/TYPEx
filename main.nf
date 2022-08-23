#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// to define help function & hlp argument
include { helpMessage } from './lib/core_functions'
include { MAJOR_ONLY } from './workflows/major_only.nf' 

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


workflow {

	MAJOR_ONLY()
		
}


