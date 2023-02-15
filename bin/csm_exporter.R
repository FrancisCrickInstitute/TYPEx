#!/usr/bin/env Rscript
library(plyr)

#if(!require(fst))
#	install.packages('fst')
source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/utilities.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/imc_utils.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/celltype_tree_utils.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/conf/settings.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/assigner.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/summary.R'))

arg_parser=argparser::arg_parser("Summarize celltype-specific mask results")
add=argparser::add_argument
arg_parser=add(arg_parser, "--inDir", help="Input directort")
arg_parser=add(arg_parser, "--regFile",
	help="Input directory")
arg_parser=add(arg_parser, arg="--subset", 
	default="major", help=paste("major", "subtypes", "sampled"))
arg_parser=add(arg_parser, "--imcyto_run", 
	default="publication", help="NextFlow run")
arg_parser=add(arg_parser, "--panel",
	default="p2", help="Panel of markers")
arg_parser=add(arg_parser, "--method",
	default="csm", help="Method")
arg_parser=add(arg_parser, "--nndist",
	default=5, help="Distance from centroid")
arg_parser=add(arg_parser, "--markers",
	default="major_markers", help="Subset of markers")
arg_parser=add(arg_parser, "--tissAreaDir",
	default=NULL, help="Subset of markers")

args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))
nndist=args$nndist


fileNamePattern="^([^-]+)-(roi_[^_]+).*\\.([^.]+).ratios.csv"
files=list.frec(args$inDir, pattern=fileNamePattern)
cat("Number of files", length(files), "\n")
# Typing output format: 
# "file", "object", "cluster", "runID", "names", "centerX", "centerY", "panel"

data=vector(mode="list")
stats=vector(mode="list")
percell=vector(mode="list")

cat("Loading files", paste(args), "\n")
for(file in files)  {
	
  if(! which(files == file) %% 1000) print(which(files == file))
  marker=gsub(fileNamePattern, "\\3", basename(file))
  if(!marker %in% unlist(marker_gene_list[[args$markers]])) next  
  sample=gsub(fileNamePattern, "\\1-\\2", basename(file))
  runID=with(args, paste1(imcyto_run, panel, method, nndist))
  
  if(file.info(file)[["size"]] == 0) {
    data[[paste1(runID, sample)]]=c(imagename=sample, object=NA, 
                                    panel = args$panel, 
									cluster=marker, 
									cellType=NA, 
									runID=runID, 
                                    names=marker, majorType=NA, positive=NA, 
                                    positiveFiltered=NA, positiveSpecific=NA, 
                                    centerX=NA, centerY=NA,
                                    NN=NA, rawCenterX=NA, rawCenterY=NA)
    next
  }
  tmp=read.csv(file, sep=",")
  out=data.frame(
	  imagename=sample, object=tmp$ObjectNumber, panel=args$panel, 
      cluster=marker, runID=runID, names=marker,  positive=marker,
      centerX=tmp$Location_Center_X, centerY=tmp$Location_Center_Y, 
      NN=tmp$NN, rawCenterX=tmp$X, rawCenterY=tmp$Y)
  data[[args$panel]][[args$imcyto_run]][[paste(sample, marker)]]=as.data.frame(out)
  stats[[args$panel]][[args$imcyto_run]][[paste(sample, marker)]]=
    ddply(tmp, .(sample, marker),summarize,
          count=length(marker),
          overlap=length(marker[NN < nndist]))
}

print("Exporting data")
for(panel in names(data))	{
	
  for(run in names(data[[panel]]))  {
		
    runID=file.path(args$inDir, "nearest")
    dfMrg=data.table::rbindlist(data[[panel]][[run]], use.names = T)
    data.table::fwrite(dfMrg, file = f("{runID}.mrg.txt"), sep = "\t")
    fst::write_fst(dfMrg, path = f("{runID}.mrg.fst"))
		
    dfCell=dfMrg[, list(cluster = paste1(sort(cluster[NN < args$nndist])),
                        names=paste1(sort(names[NN < args$nndist])),
                        NN=paste1(NN)),
          by = c("imagename", "object", "panel", "runID",  "centerX", "centerY")]
    dfCell[, positive := names]
    
    print("Getting tissue category")
    tissue_categs=get_tissue_category(
			cellIDs = with(dfCell, paste(object, imagename)),
			panel = args$panel,
			tissAreaDir=args$tissAreaDir,
			class="regional")
    celltypes=assign_celltype(
		names=dfCell$names, 
		markers = marker_gene_list[[args$markers]])
    dfCell=data.frame(
		dfCell[, c("imagename", "object", "panel", "cluster")],
        cellType         =celltypes,
        dfCell[, c("runID", "names")],
        majorType        =NA,
        positive         =NA,
        dfCell[, c("centerX", "centerY")],
		area             =NA,
        meanIntensity    =NA,
        probability      =NA,
        tissue_categs,
        NN=dfCell$NN
    )
    
    fst::write_fst(dfCell, path = f("{runID}.clusters.fst"), compress=75)
    data.table::fwrite(dfCell, file = f("{runID}.clusters.txt"), sep = "\t")
    cat("Output saved in", runID, "\n")
		

	metaDf=read.csv(args$regFile, sep = "\t", stringsAsFactors = F)
	if("useImage" %in% colnames(metaDf)) {
			
      controls=unique(metaDf$useImage)
      for(subset in controls) {
				
        IDs=cellTypeReview$imagename[metaDf$useImage == subset]
        if(! sum(dfCell$imagename %in% IDs)) 
			next
		fst::write.fst(dfCell[ dfCell$imagename %in% IDs, ], 
			path=f("{runID}.{subset}.fst"), compress = 75)
        if(subset == "review")
			write.tab(dfCell[dfCell$imagename %in% IDs, ], 
			file=f("{runID}.{subset}.txt"))
      }
    } else {
			cat('WARNING: a celltypeReviewFile has not been provided.",
				"STRATIFICATION BY CONFIDENCE WILL NOT BE done.\n')
		}
  }
}
