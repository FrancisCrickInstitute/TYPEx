#!/usr/bin/env Rscript

source(file.path(Sys.getenv('BASE_DIR'), "/conf/settings.R"))
source(file.path(Sys.getenv("BASE_DIR"), "/lib/imc_utils.R"))
source(file.path(Sys.getenv("BASE_DIR"), "/lib/utilities.R"))
source(file.path(Sys.getenv("BASE_DIR"), "/lib/assigner.R"))
source(file.path(Sys.getenv("BASE_DIR"), "/lib/celltype_tree_utils.R"))

library(plyr)

resultPattern="*\\.clusters.txt"

arg_parser=argparser::arg_parser("Select images for QC")
add=argparser::add_argument
arg_parser=add(arg_parser, arg="--inDir", help="in")
arg_parser=add(arg_parser, arg="--outDir", default = 'qc', help="in")
arg_parser=add(arg_parser, "--panel", default="p1", help="Panel of markers")
arg_parser=add(arg_parser, "--ref_markers", default="major_markers", help="Marker lists")

args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))
pars=c(pars, args) 

if(! dir.exists(args$outDir))
	dir.create(args$outDir)
qcOut = f("{args$outDir}/overlay_examples.txt")


densityFile = f("{args$inDir}/cell_density_{args$panel}.txt")
dfIn = read.csv(densityFile, sep = '\t')
 #dfIn$majorType[is.na(dfIn$majorType)] = dfIn$cellType[is.na(dfIn$majorType)]

markerStats = marker_gene_list[[args$ref_markers]] %>% 
	unlist %>% table
print(markerStats)

selected = vector(mode = 'list')
for(majorCellType in unique(dfIn$majorType)) {
	
	if(is.na(majorCellType)) 
		next
	
	if(majorCellType %in% c('Unassigned', 'Ambiguous')) 
		next
	markers = get_celltype_markers(marker_gene_list[[args$ref_markers]], majorCellType)
	# get the most unique marker per cell type
	sel = markerStats[markers] == 1 &  markers %in% names(markerStats)
	if(! any(sel))
		sel = order(markerStats[markers])[1]
	markers = markers[sel]
	if(is.na(sel)) 
		next
	
	# for each of these markers with that cell type, select the image where this marker is most frequent
	dfSub = subset(dfIn, majorType == majorCellType)
	
	stats = ddply(dfSub, .(imagename), summarise, count = sum(cellCount))
	stats = stats[order(stats$count, decreasing = T), ]
	
	selected[[majorCellType]] = cbind(imagename = stats$imagename[1], 
									  marker=markers, majorType = majorCellType)
}

mrg = do.call(rbind, selected)
write.table(mrg, file = qcOut, quote = F,
	col.names = F, row.names = F, sep = '\t')
