#!/usr/bin/env Rscript

baseDir=Sys.getenv("BASE_DIR")
source(glue::glue("{baseDir}/lib/imc_utils.R"))
source(glue::glue("{baseDir}/conf/settings.R"))

library(raster)
library(knitr)
library(tibble)
library(magrittr)
library(glue)
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

args = argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))
f <- glue::glue

releaseOut=with(args, f("cell_objects_{run}_{panel}.fst"))
print(file.path(args$inDir, releaseOut))
data = fst::read.fst(file.path(args$inDir, releaseOut))

if(! dir.exists(args$outDir))
	dir.create(args$outDir,  recursive = T)

imagenames = unique(data$imagename)

imgMapDir = f("{args$outDir}/maps")
if(! dir.exists(imgMapDir))
	dir.create(imgMapDir)


celltypes = unique(data$cellType)
celltypes = celltypes[order(match(celltypes, names(palette$cellTypeColors)))]

cellTypeColors = sapply(celltypes, function(type) {
  print(type)
  if(type %in% names(palette$cellTypeColors))
	  return(palette$cellTypeColors[[type]])
  return('grey')
})
print(cellTypeColors)
pdfOut = f("{imgMapDir}/legend.pdf")
#png(pngOut, width = 300, height = length(cellTypeColors) * 30, units = 'px')
pdf(pdfOut)
plot(1:length(cellTypeColors), rep(1, length(cellTypeColors)), axes = F, 
	xpd = F, type = "n", xaxt = 'n', yaxt='n')
legend("left", fill = cellTypeColors, legend = names(cellTypeColors),
			pch = 21, box.lty = 0, title = "Cell subtypes", cex = 1)
dev.off()

for(imagename in imagenames) {

	print(imagename)
	dataFlt = data[data$imagename == imagename, ]
	outlineImg = list.files(file.path(args$maskDir, "simple_segmentation/"), 
		pattern = '.*nuclear_mask_nuclear_dilation.tiff', 
		full.name = T, recursive = T)
	if(args$mccs)
		outlineImg = list.files(file.path(args$maskDir, "consensus_cell_segmentation"), 
  					pattern = "total_cells_mask.tiff", recursive = T, full.name = T)
	outlineImg = grep(gsub('-',  '.', f("{imagename}\\/")), outlineImg, value = T)
  	print(outlineImg)
  	print(file.exists(outlineImg))
  	r <- raster::raster(outlineImg)
  
	typeMat = r
	typeMat[typeMat == 0] = "black"
  
	for(cellType in celltypes) {
	  objects = dataFlt$object[dataFlt$cellType == cellType]
	  color = palette$cellTypeColors[[cellType]]
	  if(is.null(color))
		  color = 'grey'
	  typeMat[typeMat %in% objects] = color
	}

	pngOut = f("{imgMapDir}/types_{imagename}.png")
	png(pngOut, width = dim(typeMat)[1], height = dim(typeMat)[2], units = 'px')
	par(omi=c(0,0,0,0), mgp=c(0,0,0), mar = c(0, 0, 0, 0), family = 'D')
	plot(typeMat, legend=F, axes = F,
	   maxpixels=dim(r)[2] * dim(r)[1], alpha=1,
	   col=cellTypeColors, xpd = T)
	dev.off()
	
}
print(f('Cell type maps saved in {imgMapDir}'))

qcDf = read.csv(args$posFile, sep = '\t', header = F)
imagenames = unique(qcDf[, 1])

for(imagename in imagenames) {
	
	print(f("Outline overlays for {imagename}"))
	dataFlt = data[data$imagename == imagename, ]
	outlineImg = list.files(file.path(args$maskDir, "simple_segmentation/"), 
		pattern = '.*nuclear_mask_nuclear_dilation.tiff', 
		full.name = T, recursive = T)
	print(outlineImg)
	print(file.exists(outlineImg))
	if(args$mccs)
		outlineImg = list.files(file.path(args$maskDir, "consensus_cell_segmentation"), 
  					pattern = "total_cells_mask.tiff", recursive = T, full.name = T)
	outlineImg = grep(gsub('-',  '.', imagename), outlineImg, value = T)
  	
  	r <- raster::raster(outlineImg)
 
	markersList = unique(qcDf[qcDf[, 1] == imagename, 2])
	for(marker in markersList) {

	  print(marker)
	  file = list.files(args$rawDir, pattern = f("^{marker}\\..{imagename}.*.png"), full.names = T)
	  if(! length(file))
		stop('File not found', paste(file, marker, imagename, args$rawDir))
	
	  dataFlt$type = get_marker_frequency(dataFlt, marker, "positive")
	  positive = grep("\\+", dataFlt$type)
	
	  if(! length(positive))
		   next
	  print(file)
	  print(file.exists(file))
	  raw = raster::raster(file)
	 
	  pos = r
	  pos[! pos %in% data$object[positive]] = NA
	
	  pngOut = f("{args$rawDir}/outlines_{basename(file)}")
	  # if(file.exists(pngOut)) next
	  png(pngOut, width=dim(raw)[1], height=dim(raw)[2], units = 'px')
	  par(omi=c(0,0,0,0), mgp=c(0,0,0),mar=c(0,0,0,0), family = 'D') 
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
