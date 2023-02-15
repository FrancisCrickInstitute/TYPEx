#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// to define help function & hlp argument
include { helpMessage } from './lib/functions.nf'
include { TYPEx } from './workflows/typex.nf' 

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


workflow {

	TYPEx()
		
}


