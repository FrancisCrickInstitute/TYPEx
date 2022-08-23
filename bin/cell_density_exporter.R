#!/usr/bin/env Rscript

source(file.path(Sys.getenv('BASE_DIR'), "/conf/settings.R"))
source(file.path(Sys.getenv("BASE_DIR"), "/lib/imc_utils.R"))
source(file.path(Sys.getenv("BASE_DIR"), "/lib/utilities.R"))

library(plyr)

resultPattern="*\\.clusters.txt"
minNrCells=100

arg_parser=argparser::arg_parser("Summarize typing results")
add=argparser::add_argument
arg_parser=add(arg_parser, arg="--inDir", default="output", help="in")
arg_parser=add(arg_parser, arg="--outDir", default="output", help="in")
arg_parser=add(arg_parser, arg="--subset", default="subtypes", help=paste("all", "sampled"))
arg_parser=add(arg_parser, arg="--method", default="Rphenograph", help="Rphenograph/flowSOM/kmeans")
arg_parser=add(arg_parser, "--run", default="final", help="NextFlow run")
arg_parser=add(arg_parser, "--study", default="tracerx", help="NextFlow run")
arg_parser=add(arg_parser, "--panel", default="p2", help="Panel of markers")
arg_parser=add(arg_parser, "--cohort", default="tx100", help="NextFlow run")
arg_parser=add(arg_parser, "--markers", default="all", help="Marker lists.R")
arg_parser=add(arg_parser, "--ref_markers", default="mcsa", help="Marker lists")
arg_parser=add(arg_parser, "--regFile", default="../../../data/metadata.tracerx.txt", help="Marker lists")
arg_parser=add(arg_parser, "--cellAssignFile", help="List with reassignments")
arg_parser=add(arg_parser, "--tissAreaDir", help="List with reassignments")


args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))

inDir=file.path(args$inDir, with(args, f(analysisPath)))
outDir=file.path(args$outDir, gsub("\\/", "_", with(args, f(analysisPath))))
if(!dir.exists(outDir)) dir.create(outDir)
imgData=read.csv(args$regFile, sep = '\t', stringsAsFactors = F)
analysisID=gsub("/", "_", with(args,f(analysisPath)))
resultFiles=list.files(inDir, pattern=resultPattern, recursive=T)

if(!length(resultFiles)) cat("No files for", with(args, f(analysisPath)), " in ", inDir, "\n")

for(file in resultFiles) {
  
  dfIn=data.table::fread(f("{inDir}/{file}"))
  cat("Processing file", inDir, file, dim(dfIn), "\n")
  fileOut=gsub(".clusters", ".cell_stats", file)
  analysisID=file.path(outDir, gsub(".clusters.*", "", file))
  
  # OUTPUTS
  releaseOut=with(args, f("cell_objects_{cohort}_{run}_{panel}.txt"))
  releaseFst=gsub(".txt", ".fst", releaseOut)
  dpOut=f("{analysisID}.DP.txt")
  positiveOut=f("{analysisID}.positivity.txt")
  densityOut=f("{analysisID}.cell_density.txt")
  densityCategsOut=f("{analysisID}.cell_density.categs.txt")
  unassignedOut=f("{analysisID}.unassigned.txt")
  
  
  dfIn$OldCellType=dfIn$cellType
  dfIn$OldMajorType=dfIn$majorType
  dfIn$region=gsub('_.*', '', dfIn$region)
  
  # Review cell type assignments
  # cellTypes = dfIn$cellType; majorTypes = dfIn$majorType; positivity = dfIn$positive; panel=args$panel
  # cohort=args$cohort;study=args$study
  dfIn$cellType=review_cellType_by_major(cellTypes = dfIn$OldCellType, majorTypes = dfIn$OldMajorType,
                                         positivity = dfIn$positive, panel=args$panel,
                                         study=args$study, cohort=args$cohort, cellAssignFile=args$cellAssignFile)
  dfIn$majorType=review_major_by_cellType(cellTypes = dfIn$OldCellType, majorTypes = dfIn$OldMajorType,
                                         positivity = dfIn$positive, panel=args$panel,
                                         study=args$study, cohort=args$cohort, cellAssignFile=args$cellAssignFile)
  dfIn$majorType = gsub(' Region', '', dfIn$majorType)
  # write for each cohort
  cohorts=imgData$cohort[match(dfIn$imagename, imgData$imagename)]
  if(any(is.na(cohorts)) & args$cohort=='peace') {
    cohorts[is.na(cohorts)] = 'peace'
  } else if(any(is.na(cohorts)) & args$cohort=='origpanel') {
    cohorts[is.na(cohorts)] = 'ASB'
  }
  for(cohort in unique(cohorts)) {

	print(cohort)
    cohortSel=which(cohorts==cohort)
    dfSub=dfIn[cohortSel, ]
    dfSub$positive=gsub('pos:(.*) neg:', '\\1', dfSub$positive)
    
    write.table(dfSub[, c("imagename", 'object',  "centerX", "centerY", 'region', "majorType", "cellType", 'positive')],
                file = f("{inDir}/{cohort}_{releaseOut}"), sep = '\t', row.names = F)
    fst::write.fst(dfSub[, c("imagename", 'object', "centerX", "centerY", 'region', "majorType", "cellType", 'positive')],
              path = f("{inDir}/{cohort}_{releaseFst}"))  
    
    if(file.exists(args$cellAssignFile)) {

      guide=read.csv(args$cellAssignFile, sep='\t', stringsAsFactors = F)
      ind=match(paste(dfSub$OldCellType, dfSub$OldMajorType, dfSub$positive),
                with(guide, paste(cellType, majorType, positivity)))
      dfSub$Review=guide$Review[ind]
      dfSub$Comment=guide$Comment[ind]
      dfSub$clusterMajor=gsub(' [0-9]+$', '', dfSub$cluster)

		
      phenoStatsDf=ddply(dfSub, .(positive, OldCellType, OldMajorType, majorType, cellType, Review, Comment), summarise,
                         mostLikelyCelltype=paste(unique(clusterMajor), collapse='|', sep = ' '),
                         cellCount=length(imagename), clusters = paste(unique(clusterMajor), collapse='_', sep = ' '))
      phenoStatsDf=phenoStatsDf[order(phenoStatsDf$cellType, phenoStatsDf$cellCount, decreasing = F), ]
      phenoStatsDf$mostLikelyCelltype=gsub('Excluded', 'LowerConf', phenoStatsDf$mostLikelyCelltype)

      phenoStatsDf$cellType_markers=sapply(1:nrow(phenoStatsDf), function(x) {

        if(! phenoStatsDf$cellType[x] %in% names(marker_gene_list$export)) return("")
        celltype=phenoStatsDf$cellType[x]
        markers=marker_gene_list$export[[celltype]]
        # Only used for the model
        paste(markers, sep = ",", collapse=",")
      })

      write.tab(phenoStatsDf, file=f("{inDir}/check_{cohort}_{args$run}_{args$panel}.txt"))
	  print(colnames(phenoStatsDf))
	  phenoStatsSub = subset(phenoStatsDf, select =  c('positive', 'majorType', 'cellType', 'Comment', 'mostLikelyCelltype', 'cellType_markers'))
      write.tab(phenoStatsSub, file=f("{inDir}/phenotypes_{cohort}_{args$run}_{args$panel}.txt"))
      #stopifnot(!any(is.na(phenoStatsDf$Review)))

    }
  }
    
  # Tissue area
  imagenames=unique(dfIn$imagename)
  tissueArea=get_tissue_area(samples=imagenames, panel=args$panel,
						study=args$study, cohort=gsub('_mcsa', '', args$cohort),
						tissAreaDir=args$tissAreaDir)
  print(table(dfIn$region))
  regionArea=vector(mode='list')
  for(region in unique(dfIn$region)) {
	print(region)
    if(is.na(region)) next
    regionArea[[region]]=data.frame(Area=get_tissue_area(samples=imagenames, panel=args$panel,
                                                         study = args$study, cohort=args$cohort,
									category = region, tissAreaDir=args$tissAreaDir),
                                    Region=region, imagename=imagenames)
  }
  regionArea=do.call(rbind, regionArea)
  print(head(regionArea))
  cell_stats=vector(mode="list")
  fileOut=gsub(".clusters", ".cell_stats", file)
  fileOut=f("summary_{fileOut}")
  
  fileFst=gsub(".txt", ".fst", fileOut)
  
  id.vars=c('imagename', 'panel', 'runID')
  dfIn$positivity=gsub("pos:(.*) neg:.*", "\\1", dfIn$positive)
  stopifnot(! any(is.na(dfIn$cellDensity)))

    id.vars=c(id.vars, "cellType", "majorType", "cluster", "positivity")
    tissue.vars=c(id.vars, "cellType", "majorType", "cluster", "positivity", 'region')

  rowSel=1:nrow(dfIn)

  print("Summarize by number of CellTypes")
  dtIn=data.table::setDT(dfIn[rowSel, ])
  dfStats=dtIn[, .(cellCount=.N), by=id.vars]
  print("Add tissue area")
  dfStats[, cellDensity:=cellCount/tissueArea[match(imagename, names(tissueArea))]]
  
  print("Cell percentages")
  cellCountAssigned=dfStats[, .(totalCellCount=sum(cellCount[! cellType %in% c('Unassigned', 'Ambiguous') ])), by = imagename]
  cat("Removing the following images with less than", minNrCells, ' cells\n')
  selectedImages=cellCountAssigned$imagename[cellCountAssigned$totalCellCount > minNrCells]
  cellCountStats=dfStats[, .(totalCellCount=sum(cellCount)), by = imagename]
  dfStats[, cellPercentage:=cellCount/cellCountStats$totalCellCount[match(imagename, cellCountStats$imagename)]]
  cell_stats[['F']] = dfStats
  
  if(!is.null(regionArea)) {
	  print("Summarize by number of CellTypes")
	  dfStatsCateg=dfIn[, .(cellCount=.N), by=tissue.vars]
	  dfStatsCateg=subset(dfStatsCateg, !is.na(region))

	  print("Add tissue area")
	  dfStatsCateg$area=regionArea$Area[match(with(dfStatsCateg, paste(region, imagename)),
                                          with(regionArea, paste(Region, imagename)))]
	  dfStatsCateg[, cellDensity:=cellCount/area]
		print(colnames(dfStatsCateg))
	  dfStatsCateg = subset(dfStatsCateg, select = c( "imagename", "panel", "cellType", "majorType", "positivity", "region", "cellCount", "cellDensity"))
	  write.tab(dfStatsCateg, file=f("{inDir}/categs_{fileOut}"))
	  fst::write.fst(dfStatsCateg, path=f("{inDir}/categs_{fileFst}"))
 }
  
  
  dfStatsMrg=do.call(rbind, cell_stats)
  dfStatsMrg$CD3=get_marker_frequency(dfStatsMrg, "CD3")
  dfStatsMrg$CD4=get_marker_frequency(dfStatsMrg, "CD4")
  dfStatsMrg$CD8a=get_marker_frequency(dfStatsMrg, "CD8a")
  dfStatsMrg$CD4_CD8=with(dfStatsMrg, paste(CD3, CD8a, CD4, sep = ""))
  print(head(dfStatsMrg))
  tCellsPct=ddply(dfStatsMrg[grep("CD4\\+|CD8a\\+", dfStatsMrg$CD4_CD8),], .(CD4_CD8), summarise, 
                  count=sum(cellCount), total=sum(dfStatsMrg$cellCount))
  tCellsPct$pct=tCellsPct$count/sum(tCellsPct$count) * 100
  print(tCellsPct)
  write.table(tCellsPct, file=dpOut, sep = "\t", quote = F, row.names = F)
  
  
  unassignedStats=ddply(dfStatsMrg, .(imagename), summarise,
		unassigned=sum(cellCount[cellType %in% c('Unassigned')]),
		ambiguous=sum(cellCount[cellType %in% c('Ambiguous', 'T cells DP')]),
		total=sum(cellCount))
  unassignedStats$upct=with(unassignedStats, unassigned/total)
  unassignedStats$apct=with(unassignedStats, ambiguous/total)
  write.table(unassignedStats, file = unassignedOut, sep="\t", row.names = F, quote=F)
	
  exportAll=summary_table(dfStatsMrg)
  write.tab(exportAll, file=densityOut)  
  
  dfStatsMrg=subset(dfStatsMrg, imagename %in% selectedImages)
  print(colnames(dfStatsMrg))
  dfStatsMrg = subset(dfStatsMrg, select = c('imagename', 'panel', 'cellType', 'majorType', 'positivity', 'cellCount', 'cellDensity', 'cellPercentage'))
  write.tab(dfStatsMrg, file=f("{inDir}/{fileOut}"))

  markers=setdiff(unique(unlist(lapply(unique(dfStatsMrg$positivity),
                                       function(x) strsplit(toString(x), split="_")[[1]]))), c('NA', ''))
  dfMarker=sapply(markers, function(marker) {
    get_marker_frequency(dfStatsMrg, marker)
  })
  dfMarker=data.frame(imagename=dfStatsMrg$imagename, 
                      cellCount=dfStatsMrg$cellCount, 
                      cellDensity=dfStatsMrg$cellDensity,
                      cellType=dfStatsMrg$cellType,
                      dfMarker,
                      majorType=dfStatsMrg$majorType, stringsAsFactors = F)
  dfMarkMelt=reshape2::melt(dfMarker, id.vars=c("imagename", "majorType", "cellCount", "cellDensity"))
  dfMarkStats=ddply(dfMarkMelt, .(imagename, variable, value), summarise,
                    count=sum(cellCount, na.rm = T),
                    density = sum(cellDensity, na.rm = T))
  write.tab(dfMarkStats[grep('\\+$', dfMarkStats$value), ], file=positiveOut)
	
  if(nrow(dfStatsMrg) >0 ) fst::write.fst(dfStatsMrg, path=f("{inDir}/{fileFst}"))
  cat("Output saved in", f("{inDir}/{fileFst}"), "\n")
  
}
print(analysisID)

