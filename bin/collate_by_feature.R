#!/usr/bin/env Rscript

# Collate files by feature
library(tidyr)
source(glue::glue(Sys.getenv("BASE_DIR"), '/conf/settings.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/utilities.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/imc_utils.R'))


arg_parser=argparser::arg_parser("Collate image info by feature")
add=argparser::add_argument
arg_parser=add(arg_parser, "--inDir", 
	default="tracerx", help="Panel of markers")
arg_parser=add(arg_parser, "--panel",
	default="p1", help="Panel of markers")
arg_parser=add(arg_parser, "--run",
	default="final", help="NextFlow run")	
arg_parser=add(arg_parser, "--feature",
	default="LocationCenter", help="Panel of markers")

args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))

for(feature in args$feature) {

	outDir=feature
	if(!dir.exists(outDir))
		dir.create(outDir)
	
	inData=load_files(
		inDir=file.path(args$inDir, feature),
		run=args$run, 
		col.exclude = pars$channels_exclude)
	
	fileOut=f("{feature}/{args$panel}_{feature}_{args$run}")
	inData$imagename=basename(inData$imagename) %>%  gsub(".txt", "", .)
	colnames(inData)=gsub("^file$", "imagename", colnames(inData))
	
	idCols=c('ObjectNumber', 'imagename')
	columnsOrder= c(idCols, colnames(inData)[! colnames(inData) %in% idCols])
	inData=inData[ , ..columnsOrder ]
	write.table(inData, file  = f("{fileOut}.csv"), 
		sep = ",", quote = F, row.names = F)
	fst::write.fst(inData, path = f("{fileOut}.fst"), compress = 75)

}

