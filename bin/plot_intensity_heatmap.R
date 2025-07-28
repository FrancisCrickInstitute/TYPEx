#!/usr/bin/env Rscript
maxDim = 2000
bgColor = "black"
reducetime = 2

options(warn = -1)
library(scales)
library(plyr)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(ggbeeswarm)
suppressPackageStartupMessages(library(ComplexHeatmap))

f <- glue::glue
split_by_other <- function(value) {
	gsub(" Tissue Segmentation.*", "\nTissue Segmentation", value, ignore.case = T)
}
source(file.path(Sys.getenv('BASE_DIR'), "/conf/settings.R"))
source(file.path(Sys.getenv('BASE_DIR'), "/lib/utilities.R"))
source(file.path(Sys.getenv('BASE_DIR'), "/lib/imc_utils.R"))
source(file.path(Sys.getenv('BASE_DIR'), "/lib/plotter.R"))
source(file.path(Sys.getenv('BASE_DIR'), "/lib/celltype_tree_utils.R"))

arg_parser=argparser::arg_parser("Summarise typing results")
add=argparser::add_argument
arg_parser=add(arg_parser, arg="--inDir", help="in")
arg_parser=add(arg_parser, arg="--subset", default="subtypes", help=paste("all", "sampled"))
arg_parser=add(arg_parser, arg="--method", default="FastPG", help="FastPG/Rphenograph/flowSOM/kmeans")
arg_parser=add(arg_parser, "--panel", default="p1", help="Panel of markers")
arg_parser=add(arg_parser, "--markers", default="subtype_markers", help="Marker lists.R")
arg_parser=add(arg_parser, "--regFile", help="Marker lists")
arg_parser=add(arg_parser, "--ref_markers", default="major_markers", help="Marker lists.R")

args = argparser::parse_args(arg_parser, argv = commandArgs(trailingOnly=TRUE))

expFile = with(args,list.files( f('{inDir}/features/MeanIntensity/'),
			pattern = f('{panel}_MeanIntensity_.*.fst'), full.names = T))
cat('Loading', expFile, '\n')
expDf = fst::read.fst(expFile)
					
# FORMAT NAMES
metals = colnames(expDf) %>% 
			gsub(paste0(".*", pars$features, "_?_(.*)"), "\\1", .) %>%
			gsub('^([0-9]+)([A-Za-z]{1,2}).*', '\\2\\1', .)

metalPattern =paste0(metals, collapse = ".*|\\.?") %>%
				paste0(., ".*")
colnames(expDf) = colnames(expDf) %>%
				gsub(paste0(".*", pars$features, "_?_(.*)"), "\\1", .) %>%
				gsub('^[0-9]+[A-Za-z]+_(.*)', '\\1', .)
if(! any(metals %in% colnames(expDf))) {
	print(metals[metals %in% colnames(expDf)])
	colnames(expDf) = colnames(expDf) %>%
						gsub(metalPattern, '', .)
}
expDf = subset(expDf, select = ! colnames(expDf) %in% pars$channels_exclude)

summaryDir = with(args, f('{inDir}/summary/{subset}_{markers}_{method}/tables'))
outDir = with(args, f("{inDir}/summary/{subset}_{markers}_{method}/intensity_plots"))

if(! dir.exists(outDir))
	dir.create(outDir, recursive = T)

print(summaryDir)

cellObjectsFile = list.files(summaryDir, 
				pattern = f("cell_objects_.*{args$panel}.fst"),
				full.names = T)
if(! length(cellObjectsFile)) {
	stop("ERROR: the cell objects table has not been created by the process subtype_exporter.")
} 
cellObjectsDf = fst::read.fst(cellObjectsFile)
expDf = subset(expDf, imagename %in% unique(cellObjectsDf))

markers = setdiff(colnames(expDf), c('ObjectNumber', 'imagename'))

densMatch = match(with(expDf, paste(imagename, ObjectNumber)),
                  with(cellObjectsDf, paste(imagename, object)))
cellTypes = cellObjectsDf$cellType[densMatch]
positive = cellObjectsDf$positive[densMatch]


cellTypeList = get_celltypes(marker_gene_list[[args$markers]])
cellTypeList = c(cellTypeList, cellTypes) %>% unique %>% 
	sapply(., toString)
cellTypeColors = rep('grey', length(cellTypeList))
names(cellTypeColors) = cellTypeList
paletteMatch = match(names(cellTypeColors), names(palette$cellTypeColors))
cellTypeColors[! is.na(paletteMatch)] = palette$cellTypeColors[paletteMatch[! is.na(paletteMatch)]]

markersList = sapply(cellTypeList, function(celltype) {
	get_celltype_markers(marker_gene_list[[args$markers]], celltype)
}) %>% unlist %>%  unique
markersList = intersect(markersList, colnames(expDf))

# Cluster summary for typing heatmaps
clusterSummary <- sapply(markers, function(x) {
  tapply(t(expDf[[x]]), as.factor(cellTypes), median)
})

clusterSize = table(cellTypes) %>% as.data.frame
clusterNorm = apply(clusterSummary, 2, to_zscore)

clusterPositive = sapply(rownames(clusterNorm), function(x) {
  print(x)
  densityFlt = subset(cellObjectsDf, cellType == x)
  sapply(markers, function(m) {
    pos = get_marker_frequency(data=densityFlt, marker = m, column = "positive")
    (grepl('\\+', pos ) %>% sum)/length(pos) * 100
  })
})

clusterSize = ddply(cellObjectsDf, .(cellType), summarise,
                   Freq = length(object))

if(max(abs(clusterNorm), na.rm = T) < 3) {
  density_scale = pretty(c(-2, 2))
  colRange=rev(brewer.pal(name = "RdBu", 11))[c(2, 4, 6, 8, 10)]
} else {
  density_scale = pretty(c(-3, 3))
  colRange=rev(brewer.pal(name = "RdBu", 11))[c(2, 3, 4, 6, 8, 9, 10)]
}
print("Plot the typing heatmap")
pdfOut = f("{outDir}/intensity_heatmap.pdf")
pdf(pdfOut, useDingbats = F, height = 10, width = 10)
  
clusterNormSub = clusterNorm[rownames(clusterNorm) != 'Excluded', ]
if(is.null(nrow(clusterNormSub))) next

rowOrder = order(match(rownames(clusterNormSub), cellTypeList))
colOrder = order(match(colnames(clusterNormSub), markersList))
clusterNormSub = clusterNormSub[rowOrder, colOrder]
clusterSizeSub = clusterSize[rowOrder, ]

heat = ComplexHeatmap::Heatmap(
	clusterNormSub[, colnames(clusterNormSub) %in% markersList],
	na_col = "grey90",
	cluster_columns = F,
	cluster_rows = F,
	row_title_gp = gpar(fontsize = 10),
	rect_gp=gpar(col = 'grey85'),
	col = circlize::colorRamp2(density_scale, colRange),
	row_title_rot = 90,
	row_title = 'Raw pixel intensity',
	border = "darkgrey",
	column_title_gp = gpar(fontsize = 10),
	row_names_gp = gpar(fontsize = 7),
	column_title = paste(args$panel),
	width = unit(length(markersList)/10, 'in'),
	height = unit(nrow(clusterNormSub)/10, "in"),
	clustering_method_rows = "ward.D2",
	# left_annotation = row_ha,
	heatmap_legend_param = list(title = "Median intensity\n[z-score]",
	                            ncol=1, nrow=4, by_row=T,
	                            direction='horizontal', fontsize=8),
	show_heatmap_legend = T,
	column_names_gp = gpar(fontsize = 8, angle=90))
ht = draw(heat)

markerOrder = colnames(clusterNormSub)[column_order(ht)]
cellTypeOrder = rownames(clusterNormSub)[row_order(ht)]

clusterPositive = clusterPositive[markerOrder, cellTypeOrder]
print(table(clusterPositive))
  if(! all(clusterPositive == 0) & ! all(clusterPositive == 1)) {
    
    print('Plotting bin mat')
    bin = ComplexHeatmap::Heatmap(t(clusterPositive), 
                                cluster_rows = F,
                                cluster_columns = F,
                                row_title_gp = gpar(fontsize = 10),
                                column_title = NULL,
								row_title = 'Positivity',
								row_title_rot=90,
                                rect_gp = gpar(col = 'grey85'),
                                column_title_gp = gpar(fontsize = 10),
                                row_names_gp = gpar(fontsize = 7),
                                column_names_gp = gpar(fontsize = 8),
                                border = "grey85",
                                na_col = 'white',
                                width = unit(length(markersList)/10, 'in'), 
                                height = unit(nrow(clusterNormSub)/10, "in"),
                                col = colorRampPalette(rev(brewer.pal(name = "RdBu", 5)[-c(2, 4, 5)]))(100),
								heatmap_legend_param = list(title = "Expressed marker", direction = 'horizontal', fontsize = 8))
    ht_list=heat %v% bin
	print("Drawing ht_list")
	ht_list=draw(ht_list)
    print('Binary')
  }
	print(cellTypeColors)
  
  row_ha = HeatmapAnnotation(
    "# Cells" = anno_barplot(border = F, 
                             x = clusterSizeSub$Freq,
                             gp = gpar(fill = '#0571B0', col='transparent',
							 fontsize = 8, title = expression("# Cells"))),
    cellType = rownames(clusterNormSub),
    which = 'row',
    gp = gpar(col = "black"),
	col = list(cellType = unlist(cellTypeColors)))
  full = ComplexHeatmap::Heatmap(
    clusterNormSub,
	na_col = "grey90",
    cluster_columns = F,
    cluster_rows = F,
    row_title_gp = gpar(fontsize = 10),
    rect_gp=gpar(col = 'grey85'),
    col = circlize::colorRamp2(density_scale, colRange),
    row_title_rot = 90,
	row_title = 'Raw pixel intensity',
    border = "darkgrey",
    column_title_gp = gpar(fontsize = 10),
    row_names_gp = gpar(fontsize = 7),
    column_title = paste(args$panel),
    width = unit(ncol(clusterNormSub)/10, 'in'), 
    height = unit(nrow(clusterNormSub)/10, "in"), 
    clustering_method_rows = "ward.D2",
    left_annotation = row_ha,
    heatmap_legend_param = list(title = "Median intensity\n[z-score]",
                                ncol=1, nrow=4,by_row=T, 
                                direction='horizontal', fontsize=8),
    show_heatmap_legend=T,
    column_names_gp = gpar(fontsize = 8, angle = 90))
  full = draw(full)
  print("Full map")
dev.off()

# Intensity boxplots
print("Plotting intensity heatmaps")
cellTypeStats = ddply(cellObjectsDf, .(cellType, positive), summarise, 
					count = length(cellType))
cellTypeStats$majorType = cellTypeStats$cellType
cellTypeStats$cluster = paste(cellTypeStats$positive, cellTypeStats$cellType)
clusters = paste(cellObjectsDf$positive[densMatch], cellObjectsDf$cellType[densMatch])

plot_heatmap(expDf[, setdiff(colnames(expDf), c('imagename', 'ObjectNumber'))], 
	clusters, runID = f("{outDir}/heatmap"), cellTypeStats, 
	plotDir = f("{outDir}"), plotPos = T)

cellTypeStats = ddply(cellObjectsDf, .(cellType), summarise, 
					count = length(cellType))
cellTypeStats$majorType = cellTypeStats$majorType
cellTypeStats$cluster = cellTypeStats$cellType
cellTypeStats$positive = cellTypeStats$cellType
cellTypeStats = cellTypeStats[order(match(cellTypeStats$cellType, cellTypeList)), ]
clusters = cellObjectsDf$cellType[densMatch]

plot_expression(expDf[, unique(c(markersList, colnames(expDf)))], clusters, 
				cellTypeStats, 
				pars[[args$method]], 
				plotDir = f("{outDir}"), magnitude = max(1, pars$magnitude/10))
					
# Plotting maps
posDir = with(args, f("{inDir}/summary/{subset}_{markers}_{method}/maps/scatter"))
if(! dir.exists(posDir))
	dir.create(posDir, recursive = T) 
pngOut = f("{posDir}/legend.png")
png(pngOut, width = 300, height = length(cellTypeColors)*30, units = 'px')
par(omi = c(0,0,0,0), mgp=c(0,0,0), mar=c(0,0,0,0), family = 'D')
plot(1:length(cellTypeColors), rep(2, length(cellTypeColors)), axes = F, type = "n")
legend("left", pt.bg = unlist(cellTypeColors), legend = names(cellTypeColors),
			pch = 22, box.lty=0, title = "Cell subtypes", cex = 1.5)
dev.off()	

if(cellTypeColors["Ambiguous"] == 'black') 
	bgColor = 'white'
imagenames = unique(cellObjectsDf$imagename)
for(image in imagenames) {
    cat("Plotting scatter plot of cell types for ", image, '\n')

	dataFlt = subset(cellObjectsDf, imagename == image)
    pngOut = f("{posDir}/{image}.png")
	height = max(dataFlt$centerY)
	width = max(dataFlt$centerX)
	if(width > maxDim | height > maxDim) {
		if(width > height) {
			height = height / reducetime 
			width = width / reducetime
		} else {
			width = width / reducetime
			height = height / reducetime
		}
	}

	cat(width, height, '\n')
  	png(pngOut, width = width, height = height, units = 'px')
  	par(omi = c(0,0,0,0), mgp=c(0,0,0),mar=c(0,0,0,0), family = 'D')
	g <- ggplot(dataFlt, aes(centerX, centerY, color = cellType))
	plot = g + geom_point(size = 1) +
		coord_equal() +
		scale_y_reverse() +
		xlab("") + ylab("") +
		theme_minimal() +
		guides(color = "none") +
		scale_color_manual(values = cellTypeColors) +
		theme(legend.position = 'right',
		      legend.background = element_rect(fill=bgColor),
		      axis.line = element_blank(), axis.ticks = element_blank(), 
		      axis.text = element_blank(), axis.title = element_blank(),
		      panel.background = element_rect(fill = bgColor, colour = NA),
		      plot.background = element_rect(fill = bgColor, colour = NA),
		      plot.title = element_text(color = "white"),
		      panel.grid = element_blank(),
		      panel.spacing = unit(0, "in"),
		      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"),
		      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit="in"),
		      # legend.justification = c(0, 0),
		      # legend.key.size = unit(x = 2, units =  "cm"),
		      legend.text = element_text(color = "white", size = 20),
		      legend.spacing = unit(0.1, "in")) +
		labs(x=NULL, y=NULL) 
	print(plot)
	dev.off()
}




stats = ddply(cellObjectsDf, .(cellType), summarise, TotalFreq = length(cellType))    
stats$cellType = factor(stats$cellType, levels = cellTypeList)
pdfOut = f("{outDir}/cell_types_pie_chart_v2.pdf")
pdf(pdfOut, useDingbats = F, height = 5, width = 5)
g <- ggplot(stats, aes(x="", y = TotalFreq, fill = cellType))
plot = g + geom_bar(stat="identity", position =  'fill', color = 'black', alpha = 1) + 
  cowplot::theme_cowplot() +
  coord_polar(theta = "y", start = 0, direction = -1) + 
  theme(axis.line = element_blank(), 
	    strip.background = element_blank(),
    	axis.text = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") + ylab("") +
  scale_y_continuous(labels=function(x) {
    paste0(x * 100, "%")
  }) +
  scale_fill_manual(values = cellTypeColors[match( names(cellTypeColors), cellTypeList)])
print(plot)
dev.off()




exp = setDT(expDf)
anaMatch = match(with(exp, paste(imagename, ObjectNumber)),
                  with(cellObjectsDf, paste(imagename, object)))
exp$positive = cellObjectsDf$positive[anaMatch]
exp$cellType = cellObjectsDf$cellType[anaMatch]
exp$confidence = grepl("Excluded ", cellObjectsDf$cluster[anaMatch])
exp$confidence = exp$confidence %>% gsub("TRUE", 'low confidence', .) %>%
  gsub('FALSE', 'high confidence', .)
 
print("median per celltype")
majorMarkers = intersect(markers, unlist(marker_gene_list[[args$ref_markers]]))
width = length(unique(cellTypes))
pdfOut = f("{outDir}/median_intensities_per_celltype.log10.pdf")
pdf(pdfOut, height = 7, width = width)
for(marker in majorMarkers) {
	
  cat("Intensity violin plots for ", marker, '\n')
  markerCols = intersect(colnames(exp), markers)
  avg = exp[, lapply(.SD, median), 
  				.SDcols = markerCols,
                by = .(imagename, cellType, confidence)]
  avg$panel = "panel"
  
  count = exp[, length(unique(ObjectNumber)),
                 by = .(imagename, cellType, confidence)]
  countMatch = match(with(avg, paste(imagename, cellType, confidence)), 
  					 with(count, paste(imagename, cellType, confidence)))
  avg$count = count$V1[countMatch]
  	
  cmbDF = reshape2::melt(avg, id.vars = c("imagename", "cellType", "confidence", "count", "panel"))
  cmbDF = reshape2::dcast(cmbDF, imagename + cellType + count + confidence + variable ~ panel)
 
  totalCount = ddply(cmbDF, .(panel, cellType, confidence, variable),
                       summarise,
                       total = sum(count))
  totalMatch = match(with(cmbDF, paste(panel, cellType, confidence, variable)),
                       with(totalCount, paste(panel, cellType, confidence, variable)))
  cmbDF$total = totalCount$total[totalMatch]
  cmbDF$pct = cmbDF$count / cmbDF$total * 100
	
  cmbDF$panel = cmbDF$panel * pars$magnitude / 10
  cmbDF$panel[cmbDF$panel < 1e-2] = 1e-2

  nrCelltypes = unique(cmbDF$variable)

  # Plotting an image if there are more than 5 cells to represent that image
  sub = subset(cmbDF, variable == marker & count > 5)
  sub = droplevels(sub)
  sub$cellType = factor(sub$cellType, levels = c(cellTypeList, setdiff(unique(sub$cellType), cellTypeList)))
    if(! marker %in% sub$variable) next
    if(all(is.na(sub[sub$variable == marker, "panel"])))
      next
    g <- ggplot(sub,  aes(cellType, panel))
    plot = g + 
     geom_violin(aes(fill = cellType), # varwidth = T, 
                 position = position_dodge(width = 1), 
                 color='black') +
      geom_quasirandom(aes(fill = cellType, size = count), dodge.width = 1, 
	  			color='black', alpha = .6, pch = 21) +
      stat_summary(geom = 'pointrange', aes(fill = cellType), # varwidth = T,
                   position = position_dodge(width = 1),
                   color='black', fun = median) +
	  scale_size_continuous(breaks = c(10, 100, 400, 1000), name = '# Cells') +
      theme_classic(base_size = 16) + 
      scale_fill_manual(values = cellTypeColors) +	  
      theme(legend.position="top",
            axis.text.x = element_text(color="black", angle = 90, hjust = 1),
            axis.text.y = element_text(color="black"),
            strip.background = element_blank(),
            legend.title = element_blank()) +
      scale_y_log10() +
	  guides(fill = 'none') +
      facet_wrap(confidence ~ .,  scale = 'free_y') + 
	  xlab("") + 
	  ylab("Median cell intensity per image") +
	  scale_x_discrete(labels = split_by_other) +
      ggtitle(marker) +
      theme(strip.text=element_text(angle=0))
    print(plot)

}
dev.off()

majorMarkers = intersect(markers, unlist(marker_gene_list[[args$ref_markers]]))
width = length(unique(cellTypes))

print("median positive")
pdfOut = f("{outDir}/median_intensities_per_positive_type.log10.pdf")
pdf(pdfOut, height = 7, width = width)
for(marker in majorMarkers) {

  exp$markerPos = get_marker_frequency(data = exp, marker, 'positive')
  markerCols = intersect(colnames(exp), markers)
  avg = exp[, lapply(.SD, median), 
  				.SDcols = markerCols,
                by = .(imagename, cellType, markerPos, confidence)]
  avg$panel = "panel"
  
  count = exp[, length(unique(ObjectNumber)),
                 by = .(imagename, cellType, markerPos, confidence)]
  countMatch = match(with(avg, paste(imagename, cellType, markerPos, confidence)), 
  					 with(count, paste(imagename, cellType, markerPos, confidence)))
  avg$count = count$V1[countMatch]
  	
  cmbDF = reshape2::melt(avg, id.vars = c("imagename", "cellType", "markerPos", "confidence", "count", "panel"))
  cmbDF = reshape2::dcast(cmbDF, imagename + cellType + markerPos + count + confidence + variable ~ panel)
 
  totalCount = ddply(cmbDF, .(panel, cellType, markerPos, confidence, variable),
                       summarise,
                       total = sum(count))
  totalMatch = match(with(cmbDF, paste(panel, cellType, markerPos, confidence, variable)),
                       with(totalCount, paste(panel, cellType, markerPos, confidence, variable)))
  cmbDF$total = totalCount$total[totalMatch]
  cmbDF$pct = cmbDF$count / cmbDF$total * 100
	
  cmbDF$panel = cmbDF$panel * pars$magnitude / 10
  cmbDF$panel[cmbDF$panel < 1e-2] = 1e-2

  nrCelltypes = unique(cmbDF$variable)
 
  sub = subset(cmbDF, variable == marker & count > 5)
  sub = droplevels(sub)
    if(! marker %in% sub$variable) next
    if(all(is.na(sub[sub$variable == marker, "panel"])))
      next
    g <- ggplot(sub,  aes(cellType, panel))
    plot = g + 
     geom_violin(aes(fill = markerPos), # varwidth = T, 
                 position = position_dodge(width = 1), 
                  color='black') +
      geom_quasirandom(aes(fill = markerPos, size = count), dodge.width = 1, 
	  			color='black', alpha = .6, pch = 21) +
      stat_summary(geom = 'pointrange', aes(fill = markerPos), # varwidth = T,
                   position = position_dodge(width = 1),
                   color='black', fun = median) +
	  scale_size_continuous(breaks = c(10, 100, 400, 1000), name = '# Cells') +
      theme_classic(base_size = 16) + 
      scale_fill_manual(values = c('blue', 'red')) +
	  
      theme(legend.position="top",
            axis.text.x=element_text(color="black", angle = 90, hjust = 1),
            axis.text.y=element_text(color="black"),
            strip.background=element_blank(),
            legend.title = element_blank()) +
      scale_y_log10() +
      facet_wrap(confidence ~ .,  scale = 'free_y') + 
	  xlab("") + ylab("Median cell intensity per image") +
	  scale_x_discrete(labels = split_by_other) +
      ggtitle(marker) +
      theme(strip.text=element_text(angle=0))
    print(plot)
}
dev.off()
