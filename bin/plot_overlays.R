#!/usr/bin/env Rscript

baseDir=Sys.getenv("BASE_DIR")
source(glue::glue("{baseDir}/lib/imc_utils.R"))

library(knitr)
library(tibble)
library(magrittr)
library(glue)
library(raster)
library(tiff)
 
 
arg_parser=argparser::arg_parser("Select images for QC")
add=argparser::add_argument
arg_parser=add(arg_parser, arg="--rawDir", help="in")
arg_parser=add(arg_parser, arg="--inDir", help="in")
arg_parser=add(arg_parser, arg="--maskDir", help="in")
arg_parser=add(arg_parser, arg="--outDir", default = 'qc', help="in")
arg_parser=add(arg_parser, arg="--posFile", default = 'qc', help="in")
arg_parser=add(arg_parser, "--run", default="PHLEX_test", help="NextFlow run")
arg_parser=add(arg_parser, "--panel", default="p1", help="Panel of markers")
arg_parser=add(arg_parser, "--mccs", help="MCCS or simple segmentation")

args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))
f <- glue::glue

releaseOut=with(args, f("cell_objects_{run}_{panel}.fst"))
print(releaseOut)
data = fst::read.fst(file.path(args$inDir, releaseOut))

if(! dir.exists(args$outDir))
	dir.create(args$outDir,  recursive = T)

qcDf = read.csv(args$posFile, sep = '\t')
imagenames = unique(qcDf[, 1])
print(imagenames)


for(imagename in imagenames) {

	markersList = unique(qcDf[qcDf[, 1] == imagename, 2])
	for(marker in markersList) {

	  print(marker)
	  file=list.files(args$rawDir, pattern = f("^{marker}\\..{imagename}.*.png"), full.names = T)
	  print(file)
	  if(! length(file))
		stop('File not found', paste(file, marker, imagename, args$rawDir))
	  
	  outlineImg=list.files(
				file.path(args$maskDir, "simple_segmentation/",  gsub("-", "/", imagename)), 
			pattern = '.*nuclear_mask_nuclear_dilation.tiff', 
			full.name = T)
	  if(args$mccs)
	  	outlineImg=file.path(args$maskDir, "consensus_cell_segmentation", 
					gsub("-", "/", imagename), "total_cells_mask.tiff")
	  
	  pngOut = f("{args$outDir}/outlines_{basename(file)}")
	  # if(file.exists(pngOut)) next
	  data$type = get_marker_frequency(data, marker, "positive")
	  positive = intersect(grep("\\+", data$type), 
	  						which(data$imagename == imagename))
	  if(! length(positive))
		   next
	  raw = raster(file)
	  r <-raster(outlineImg)
	  pos = r
	  pos[! pos %in% data$object[positive]] = NA
	  # neg[neg==0] = NA
	  # neg[neg %in% data$object[positive]] = NA
	
	  png(pngOut, width=dim(raw)[1], height=dim(raw)[2], units = 'px')
	  par(omi=c(0,0,0,0), mgp=c(0,0,0),mar=c(0,0,0,0), family = 'D') 
	  # par(mar=c(0, 0, 0, 0))
	  plot(raw, axes = F, legend = F, 
		  col = colorRampPalette(c("white", "black"))(80),
		  maxpixels = dim(raw)[1] * dim(raw)[2] * 10, 
			xpd = T, frame = F, box = F, legend.mar = 0)
	  plot(pos, legend=F, axes = F,
	       maxpixels=dim(r)[2] * dim(r)[1] * 10, add = T, alpha=0.4,
	       col=hcl.colors(length(positive)), xpd = T, type = 'n')
	  dev.off()
	}
}
