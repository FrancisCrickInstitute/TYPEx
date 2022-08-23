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
source(glue::glue("{baseDir}/lib/launcher.R"))
source(glue::glue("{baseDir}/lib/assigner.R"))
source(glue::glue("{baseDir}/lib/plotter.R"))


arg_parser=argparser::arg_parser("Run a cell typing method")
add=argparser::add_argument
arg_parser=add(arg_parser, "--wDir", help="Working directory")
arg_parser=add(arg_parser, "--nfDir", help="Directory with image-level matrices")
arg_parser=add(arg_parser, "--subset", help="sampled/major/subtypes")
arg_parser=add(arg_parser, "--cohort", default="tx100", help="NextFlow run")
arg_parser=add(arg_parser, "--study", default="tracerx", help="NextFlow run")
arg_parser=add(arg_parser, arg="--method", default="flowSOM",
               help=paste("Rphenograph/flowSOM/kmeans/rtsne/umap"))
# "immunoClust", "grannet", "x-shift", "Rclusterpp", "autoencoder", "SAUCI"
arg_parser=add(arg_parser, "--run", default="final", help="NextFlow run")
arg_parser=add(arg_parser, "--panel", default="p2", help="Panel of markers")
arg_parser=add(arg_parser, "--markers", default="all",
               help="Marker lists defined in cell_type_annotation.json")

arg_parser=add(arg_parser, "--cellAssignFile", help="List with reassignments")
arg_parser=add(arg_parser, "--celltypeReviewFile", help='File with image annotations')
arg_parser=add(arg_parser, "--regFile", help="Metadata")
arg_parser=add(arg_parser, "--iter", default='1', help="Iteration for subsampling")
arg_parser=add(arg_parser, "--resampleBy", default=NULL, help="file or NULL")
arg_parser=add(arg_parser, arg="--celltypeModelFile", default="output/review", help="Review model")
arg_parser=add(arg_parser, "--major_markers", default="mcsa",
               help="Marker lists defined in cell_type_annotation.json")
arg_parser=add(arg_parser, "--stratify", default=T, help="file or NULL")
arg_parser=add(arg_parser, "--tissAreaDir", default=NULL, help="file or NULL")


args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))
pars=c(pars, args)

# Run-specific markers
cellTypeReview=read.csv(args$celltypeReviewFile, sep = "\t", stringsAsFactors = F)

if(! pars$markers %in% c(names(marker_gene_list), "all"))
  stop("Marker list -", pars$markers, "- not found")

# Set output directory
outDir=file.path(args$wDir, with(pars, f(analysisPath)))
if(args$subset == sampledDir) {
  # find iteration name
   outDir=file.path(outDir, args$iter)
#  resampleBy='file'
#  if(pars$subset == sampledDir) resampleBy = NULL
}
if(! args$stratify & args$subset == subtypesDir) {
	outDir=paste(outDir, args$stratify, sep = "_")
}

if(!dir.exists(outDir))
  dir.create(outDir, recursive=T)
logFile = file.path(outDir, "execution_time.txt")

for(feature in pars$features) {
  
  runID=file.path(outDir, paste1(pars[[pars$method]]))
  inDir=file.path(args$nfDir, feature)
  resample_frac=pars[[pars$method]]$resample_frac
  print(resample_frac)
  print(pars[[pars$method]]$resample)
  resample=F
  if('resample' %in% names(pars$method))
		resample = pars[[pars$method]]$resample
  print(f("{runID}.pars.yaml"))
  print((args$subset == sampledDir | resample) & !file.exists(f("{runID}.pars.yaml")))
  inData=load_files(inDir=inDir,
                    resample=(args$subset == sampledDir | resample) & !file.exists(f("{runID}.pars.yaml")),
                    resample_frac=pars[[pars$method]]$resample_frac,
                    run=pars$run, col.exclude = pars$channels_exclude,
                    resampleBy = args$resampleBy)
  run=pars$run; col.exclude = pars$channels_exclude
  start=Sys.time()
  cat("Running", feature, pars$method, "\n")
  print(pars[[pars$method]])
  write.table(x = data.frame(unlist(pars[[pars$method]])),
        file=file.path(paste0(runID, ".log")), append=T)
  cat("Output dir:", outDir, "\n")
  # Filter out excluded samples
  excludeIDs=cellTypeReview$imagename[which(cellTypeReview$control == "excluded")]
  imagenames=gsub(".txt", "", basename(inData$file))
  inData=inData[! imagenames %in% excludeIDs, ]
  cat("Excluded images", excludeIDs, '\n')
	print(nrow(inData)) 

  if(pars$method %in% c("MC", 'umap', 'rtsne') |  pars$subset == subtypesDir) {

    refMethod=ifelse(pars$subset == subtypesDir & !pars$method %in% c('umap', 'rtsne'), "cellassign", pars[[pars$method]]$refMethod)
	ref=args
	ref$method=refMethod
	if(pars$method %in% c('umap', 'rtsne')) {
		ref$subset=pars$subset
		ref$markers=pars$markers
	} else {
		ref$subset=majorDir
		ref$markers=pars$major_markers
	}
    clustDir=file.path(args$wDir, with(ref, f(analysisPath)))
    clustFile=paste0(pars[[refMethod]], collapse="_")
    clustFile=f("{clustDir}/{clustFile}.clusters.txt")
	print(clustFile)
    if(!file.exists(clustFile)) stop("File does not exist", clustFile, "\n")
    cat("File with clusters:", clustFile, "\n", file=file.path(paste0(runID, ".log")), append=T)
    clustData=data.table::fread(clustFile, sep = "\t")
	clustData$cluster=gsub('CD([48])$', 'CD\\1 T cells', clustData$cluster)
    if(pars$subset == subtypesDir) {
      rowSel=get_filtered_objects(clustData, subset = majorDir, run = pars$run, 
                                  markers=pars$subtype_markers, panel = pars$panel,
                                  cohort=pars$cohort, celltypeModelFile=pars$celltypeModelFile)
	  if(pars$stratify) {
		  print(table(rowSel))
		  if(any(is.na(rowSel))) {
				print('No stratification for the following cell types')
				print(sum(is.na(rowSel)))
				print(table(clustData$cluster[is.na(rowSel)]))
		  }
      	  filteredInd=which(! rowSel | is.na(rowSel))
	      clustData$cluster[filteredInd] = paste("Excluded", clustData$cluster[filteredInd])
    	  print(table(clustData$cluster))
		  print(sum(is.na(clustData$cluster)))
		print(nrow(clustData))
	  }
    }
    clustMatch=match(with(inData, paste(ObjectNumber, gsub("(.*).txt", "\\1", basename(file)))),
                     with(clustData, paste(object, imagename)))
	print(sum(is.na(clustMatch)))
	if(pars$method %in% c('umap', 'rtsne') & pars$subset == 'subtypes') {
		clusters = clustData$majorType[clustMatch]
	} else {
	    clusters=clustData$cluster[clustMatch]
	}
  } else {
    clusters="black"
  }
  print(table(clusters))

  # Running analysis
  results=run_method(inData = inData, method=pars$method, feature=feature, 
                   pars=pars, runID=runID, colors=clusters, wDir=args$wDir,
		   regFile=args$regFile, nfDir=args$nfDir)
  # method=pars$method; colors=clusters; runOutput=results$runOutput; columnNames=results$columnNames
  if(! pars$method %in% c('umap', 'tsne')) {
	  out=summarise_output(inData=inData, method=pars$method, pars=pars,
                       runID=runID,
                       runOutput=results$runOutput, 
                       columnNames=results$columnNames,
                       colors=clusters,  nfDir=args$nfDir)
  } else {
	out=NA
  }
  if(file.exists(args$celltypeReviewFile) & !is.na(out)) {
    controls=unique(cellTypeReview$control)
    for(subset in controls) {
      IDs=cellTypeReview$imagename[cellTypeReview$control == subset]
      if(! sum(out$imagename %in% IDs)) next
      fst::write.fst(out[out$imagename %in% IDs, ], path=f("{runID}.{subset}.fst"), compress = 75)
      if(subset=="review") write.tab(out[out$imagename %in% IDs, ], file=f("{runID}.{subset}.txt"))
    }
  }
  end=Sys.time()
  
  if(!file.exists(paste0(runID, ".pars.yaml")))
   rlist::list.save(list(channels_excluded  =pars$channels_exclude,
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
                   phenotypeMarkrs    =marker_gene_list[[pars$ref_markers]],
                   magnitude          =pars$magnitude,
                   ilastik            =pars$ilastik,
                   ubiquitous         =pars$ubiquitous),
              file=paste0(runID, ".pars.yaml"))
  # imgOut=paste0(runID, ".full.RData.gz")
  # save.image(file=imgOut, compress = "gzip")
  cat("Analysis run for", pars$method, feature, "in", end - start, "\n",
      file=file.path(paste0(runID, ".log")), append=T)
  cat("Analysis run for", runID, "in", end - start, "\n")
}
