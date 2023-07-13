#!/usr/bin/env Rscript


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
arg_parser=add(arg_parser, arg="--outDir", default = 'qc', help="in")
arg_parser=add(arg_parser, arg="--posFile", default = 'qc', help="in")
arg_parser=add(arg_parser, "--run", default="PHLEX_test", help="NextFlow run")
arg_parser=add(arg_parser, "--panel", default="p1", help="Panel of markers")


releaseOut=with(args, f("cell_objects_{run}_{panel}.txt"))
data = read.csv(file.path(args$inDir, releaseOut), sep = '\t')

if(! dir.exists(args$outDir))
	dir.create(args$outDir,  recursive = T)

qcDf = fst::read_fst(args$posFile, as.data.table = T)
imagenames = unique(qcDf[, 1])
markersList = unique(qcDf[, 2])
for(imagename in imagenames) {
	for(marker in markersList) {
	  print(marker)
	  file=list.files(rawDir, pattern = f("{marker}\\..{imagename}.*.png"), full.names = T)
	  pngOut = f("{outDir}/outlines_{basename(file)}")
	  # if(file.exists(pngOut)) next
	  outlineImg=file.path(maskDir, gsub("-", "/", imagename), "total_cells_mask.tiff")
	  data$type=get_marker_frequency(data, marker, "positive")
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
 
	  par(mar=c(0, 0, 0, 0))
	  png(pngOut, width=dim(raw)[2], height=dim(raw)[1], units = 'px')
	  plot(raw, axes = F, legend = F, 
		  col = colorRampPalette(c("white", "black"))(20),
		  maxpixels = dim(raw)[1] * dim(raw)[2]*10)
	  plot(pos, legend=F, axes = F,
	       maxpixels=dim(r)[1] * dim(r)[2], add = T, alpha=0.4,
	       col=hcl.colors(length(positive)))
	  dev.off()
	  browseURL(pngOut)
  
	}
}
