#!/usr/bin/env Rscript
options(warn=-1)
library(ggplot2) 
library(tidyr)
library(RColorBrewer)
suppressPackageStartupMessages(library(ComplexHeatmap))

baseDir=Sys.getenv("BASE_DIR")
source(glue::glue("{baseDir}/conf/settings.R"))
source(glue::glue("{baseDir}/lib/utilities.R"))
source(glue::glue("{baseDir}/lib/imc_utils.R"))
source(glue::glue("{baseDir}/lib/celltype_tree_utils.R"))
source(glue::glue("{baseDir}/lib/launcher.R"))
source(glue::glue("{baseDir}/lib/assigner.R"))
source(glue::glue("{baseDir}/lib/summary.R"))
source(glue::glue("{baseDir}/lib/plotter.R"))


arg_parser=argparser::arg_parser("Run a cell typing method")
add=argparser::add_argument
arg_parser=add(arg_parser, "--wDir",
	help="Working directory")
arg_parser=add(arg_parser, "--nfDir", 
	help="Directory with image-level matrices")
arg_parser=add(arg_parser, "--subset", 
	help="sampled|major|subtypes")
arg_parser=add(arg_parser, arg="--method",
	default="FastPG", 
	help="Rphenograph|flowSOM|kmeans|rtsne|umap") 
arg_parser=add(arg_parser, "--run", 
	default="final", help="NextFlow run")
arg_parser=add(arg_parser, "--panel",
	default="p2", help="Panel of markers")

## METADATA
arg_parser=add(arg_parser, "--regFile", help="Metadata")

## ASSIGNMENT
arg_parser=add(arg_parser, "--cellAssignFile", 
	help="List with reassignments")
arg_parser=add(arg_parser, "--markers",
	help="Marker list name in cell_type_annotation.json")
arg_parser=add(arg_parser, "--major_markers",
	default="major_markers", 
	help="Marker list name in cell_type_annotation.json")
arg_parser=add(arg_parser, "--tissAreaDir", 
	default=NULL, help="file or NULL")

## SUBSAMPLE
arg_parser=add(arg_parser, "--iter", 
	default='1', help="Iteration for subsampling")
arg_parser=add(arg_parser, "--resampleBy", 
	default=NULL, help="file or NULL")

# STRATIFICATION BY CONFIDENCE
arg_parser=add(arg_parser, "--celltypeReviewFile",
	help='File with image annotations')
arg_parser=add(arg_parser, arg="--celltypeModelFile", 
	help="Review model")
arg_parser=add(arg_parser, "--stratify", 
	default=T, help="file or NULL")	
arg_parser=add(arg_parser, "--stratify_label", help="file or NULL")

args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))


pars=c(pars, args, pars_stratify[[args$stratify_label]])
# Run-specific markers
cellTypeReview=read.csv(args$celltypeReviewFile, 
	sep = "\t", stringsAsFactors = F)

if(! pars$markers %in% names(marker_gene_list))
  stop("Marker list -", pars$markers, "- not found")

if(! is.null(pars_stratify$stratify_label) & pars$stratify & pars$subset == majorDir) {
	abrv=strsplit(pars$excludedCellType, split = '')[[1]][1:3] %>% 
		tolower %>% paste0(., collapse = "")
	ref_markers_list=paste1(pars$markers, abrv)
	truncated_list=remove_node(
		marker_gene_list[[pars$markers]], pars$stratify_label)
	marker_gene_list[[ref_markers_list]] = truncated_list 
	pars$markers=ref_markers_list
}
# Set output directory
outDir=file.path(args$wDir, with(pars, f(analysisPath)))
if(args$subset == sampledDir)
	outDir=file.path(outDir, args$iter)

if(! args$stratify & args$subset == subtypesDir)
	outDir=paste(outDir, args$stratify, sep = "_")

if(!dir.exists(outDir))
  dir.create(outDir, recursive=T)

# Path and prefix used for the files, uniquely identifying that run
runID=file.path(outDir, paste1(pars[[pars$method]]))
write.table(x = data.frame(unlist(pars[[pars$method]])),
	file=f("{runID}.log"), append=T)
cat("Output dir:", outDir, "\n")
print(pars[[pars$method]])

for(feature in pars$features) {
	
	inDir=f("{args$nfDir}/{feature}")
	resample_frac=pars[[pars$method]]$resample_frac
	resample=F
	if('resample' %in% names(pars$method))
		resample = pars[[pars$method]]$resample
	toResample=(args$subset == sampledDir | resample) &
		! file.exists(f("{runID}.pars.yaml"))
	
	inData=load_files(
		inDir=inDir,
		resample=toResample,
		resample_frac=pars[[pars$method]]$resample_frac,
		run=pars$run, 
		col.exclude = pars$channels_exclude, 
		resampleBy = args$resampleBy)
		# run=pars$run; col.exclude = pars$channels_exclude
		
	start=Sys.time()
	cat("Running", feature, pars$method, "\n")
	
	## Filter out excluded samples
	imagenames=gsub(".txt", "", basename(inData$imagename))
	excludeIDs=cellTypeReview$imagename[
		which(cellTypeReview$control == "excluded" &
			cellTypeReview$imagename %in% imagenames)]
	inData=inData[! imagenames %in% excludeIDs, ]
	cat("Kept cells", sum(! imagenames %in% excludeIDs), 
		"out of", length(imagenames), '\n')
	print(nrow(inData))
	
	if(pars$method %in% c("MC", dimred_methods) | pars$subset == subtypesDir) {
		refMethod=ifelse(
					pars$subset == subtypesDir & ! pars$method %in% dimred_methods, 
					"cellassign", 
					pars[[pars$method]]$refMethod)
			
			ref <- args; ref$method=refMethod
			if(pars$method %in% dimred_methods) {
				ref$subset=pars$subset
				ref$markers=pars$markers
			} else {
				ref$subset=majorDir
				ref$markers=pars$major_markers
			}
		clustDir=file.path(args$wDir, with(ref, f(analysisPath)))
		clustFile=paste0(pars[[refMethod]], collapse="_")

		clustFile=f("{clustDir}/{clustFile}.clusters.txt")
		if(! file.exists(clustFile)) 
			stop("File does not exist", clustFile, "\n")
		cat("INFO: File with clusters:", clustFile, "\n",
			file=f("{runID}.log"), append=T)
		clustData=data.table::fread(clustFile, sep = "\t")		
		if(pars$subset == subtypesDir) {
			
			rowSel=get_filtered_objects(clustData, 
				celltypeModelFile=pars$celltypeModelFile)

			print(which(is.na(rowSel)))
			if(pars$stratify) {
				
				if(is.null(rowSel))	{
					print('No stratification for the following cell types')
					print(sum(is.na(rowSel)))
					print(table(clustData$cluster[is.na(rowSel)]))
				} else {
					filteredInd=which(! rowSel | is.na(rowSel))
					print(clustData[filteredInd, ])
					clustData$cluster[filteredInd]=f("Excluded {clustData$cluster[filteredInd]}")
				}
			}
		}
		
		clustMatch=match(
			with(inData, paste(ObjectNumber, gsub("(.*).txt", "\\1", basename(imagename)))),
			with(clustData, paste(object, imagename)))
			print(sum(is.na(clustMatch)))
			
			
		if(pars$method %in% dimred_methods & pars$subset == 'subtypes') {
			clusters = clustData$majorType[clustMatch]
		} else {
			clusters=clustData$cluster[clustMatch]
		}
	} else {
			clusters="black"
	}
	print(table(clusters))
	# Running analysis
	results=run_method(
		inData = inData,
		method=pars$method,
		feature=feature, 
		pars=pars, runID=runID, 
		colors=clusters,
		wDir=args$wDir,
		regFile=args$regFile, 
		nfDir=args$nfDir)
		
	if(! pars$method %in% dimred_methods) {
		out=summarise_output(inData=inData, 
			method=pars$method, pars=pars,
			runID=runID,
			runOutput=results$runOutput,
			columnNames=results$columnNames,
			colors=clusters,
			nfDir=args$nfDir)
	} else {
		out=NA
	}
	if(file.exists(args$celltypeReviewFile) & !is.na(out)) {
		controls=unique(cellTypeReview$control)
		for(subset in controls) {
			IDs=cellTypeReview$imagename[cellTypeReview$control == subset]
			if(! sum(out$imagename %in% IDs))
				next
			fst::write.fst(out[out$imagename %in% IDs, ], 
				path=f("{runID}.{subset}.fst"), compress = 75)
			if(subset == "review") 
				write.tab(out[out$imagename %in% IDs, ], 
					file=f("{runID}.{subset}.txt"))
		}
	}
	end=Sys.time()
	
	cat('Printing parameters in ', f("{runID}.pars.yaml"), '\n')
	
	if(! file.exists(f("{runID}.pars.yaml"))) {
		paramsList=list(
			channels_excluded  =pars$channels_exclude,
			additional_features=pars$additional_features,
			features           =pars$features,
			area               =pars$area_exclude,
			method             =pars$method,
			parameters         =pars[[pars$method]],
			channels_exclude   =pars$channels_exclude,
			additional_features=pars$additional_features,
			min_expression     =pars$min_expresssion,
			excludedImage      =excludeIDs,
			defaultMarkers     =pars$defaultMarkersList,
			refMarkers         =pars$subtype_markers,
			markers            =marker_gene_list[[args$markers]],
			phenotypeMarkers   =marker_gene_list[[pars$major_markers]],
			magnitude          =pars$magnitude,
			ilastik            =pars$ilastik,
			ubiquitous         =pars$ubiquitous)
			# rlist::list.save(paramsList,
				#	file           =f("${runID}.pars.yaml"))
	}
	cat(f("Analysis run {runID} finished in"),
		end - start, "\n", file=f("{runID}."), append=T)
	cat("Analysis for", runID, "finished in ", end - start, "\n")
	cat(f("Settings and parameters are saved in the log file: {runID}.log"), "\n")
}
