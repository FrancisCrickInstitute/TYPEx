#!/usr/bin/env Rscript
Sys.setenv("DISPLAY"=":0.0")
options(bitmapType="cairo-png")
library(tidyr)
library(ggplot2)
library(WGCNA)
library(gridExtra)

run='test_pub'
options(stringsAsFactors = F)
Sys.setenv('BASE_DIR' = glue::glue('/Users/angelom/Documents/projects/IMC/analyses/pipeline/{run}'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/conf/settings.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/utilities.R'))

arg_parser=argparser::arg_parser("Summarize typing results")
add=argparser::add_argument

arg_parser=add(arg_parser, arg="--wDir", default=f('/Volumes/angelom/analyses/imc/{run}/output/sampled/robustness'),
               help='Where typing results are')
arg_parser=add(arg_parser, arg="--outDir", default='analyses/robustness', help='Where typing results are')
args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))

inDir=args$wDir

outDir=f("{args$outDir}/{run}")
if(!dir.exists(outDir)) dir.create(outDir, recursive = T)

inFiles=list.files(inDir, pattern='matchedLabels.*.RData')

for(fileIn in inFiles) {
  
	print(fileIn)
	load(f("{inDir}/{fileIn}"))
	print('Loaded')
	analysisID=gsub('matchedLabels.(.*).RData', '\\1', fileIn)
	pdfOut = f("{outDir}/matchedClustersComparison.{analysisID}.pdf")
	pngOut = f("{outDir}/matchedClustersComparison.{analysisID}.tiff")
	clusterSizes = sapply(clusters, function(x) length(unique(x)))

	clusterStats = table(clusters[[1]])
	cellOrder=order(clusterStats, decreasing = T)
	clusterStats=clusterStats[cellOrder]
	clustersCmp=lapply(clusters, function(x) sum(x == clusters[[1]], na.rm = T))
	clusters = clusters[order(unlist(clustersCmp), decreasing = T)]

	colorMat = WGCNA::labels2colors(as.data.frame(clusters), 
	                                colorSeq=clusterColors[c(cellOrder, (length(cellOrder) + 1):length(clusterColors))])
	colorLabels=sapply(1:length(clusters), function(x) tapply(colorMat[,x], clusters[[x]], unique))
	colorLabels=unlist(colorLabels)
	colorLabels=colorLabels[!duplicated(names(colorLabels))]
	cellOrder = order(match(names(colorLabels), names(clusterStats)))
	colorLabels = colorLabels[cellOrder]
  colorMatVals = colorLabels
  names(colorMatVals) = colorMatVals
  
	clusterOrder = do.call(order, lapply(1:length(clusters), function(x) {
		matchInd=match(clusters[[x]], names(colorLabels))
		matchInd[is.na(clusters[[x]])] = NA
		matchInd
	}))
	
	nrCells=length(clusters[[1]])
	print(nrCells)
	clusterNames=sapply(1:length(clusters), function(x) {
		print(x)
		label=paste("Subampling", x - 1)
		if(x == 1) label = 'Reference'
		nrCells=sum(!is.na(clusters[[x]]))
		overlap = round(sum(clusters[[x]] == clusters[[1]], na.rm = T)/nrCells * 100, 2)
		paste0(label,'\n', clusterSizes[[x]], " clusters")
	})
	colnames(colorMat)=c('Ref', paste0('Run', 1:(length(clusters) - 1)))
	print(clusterNames)

	# Add cluster sizes and iteration # 
#	if(!file.exists(pngOut)) {
		print(pngOut)
		png(filename = pngOut, res = 150, width = 5, height=2, units = 'in', type = 'cairo')
		par(mar = c(0.5, 0.5, 0.5, 0.5))
		clP=plotOrderedColors(colors = colorMat, order = clusterOrder, rowLabels = names(colorMat))
		# g<-ggplot(data.frame(colorMat), aes(x = 1, Ref))
		# plot = g + geom_bar(aes(fill=Ref), position = 'fill')+ 
		#     scale_fill_manual(values = colorMatVals) + 
		#   guides(fill = 'none')
		# print(plot)
		dev.off()
    browseURL(pngOut)
#	}
 
	overlap=data.frame(t(sapply(1:length(clusters), function(x) {
		nrCells=sum(!is.na(clusters[[x]]))
		c(pct = sum(clusters[[x]] == clusters[[1]], na.rm = T) / nrCells * 100, total=nrCells, round=x, name =clusterNames[x])
	})))
	overlap$pct = as.numeric(as.character(overlap$pct))
	overlap$name = factor(overlap$name, levels = rev(overlap$name))
	pdf(pdfOut, width = 4, height = 3)
	g <- ggplot(overlap, aes(name, pct))
	plot <- g + geom_col(fill='#cedcefff', col='grey20') + #(stat='identity') + 
		cowplot::theme_cowplot() + 
		ylab(paste0("% cells in concordant clusters\n(", clusterSizes[[2]], 'cells)')) +
		xlab('') +
		coord_flip() +
		geom_text(aes(label=paste(round(pct, 1), '%' )), hjust=1)
	print(plot)
	dev.off()
	browseURL(pdfOut)
	# barplot(rev(overlap), main = '% Cells matching the reference', horiz = T, axisnames = F, col = 'aliceblue')
	# text(x = 50, y =  1:length(overlap), labels = names(rev(overlap)))
	# barplot(rev(clusterSizes), ylab = '% Rphenograph runs', horiz = T, axisnames = F, col = 'aliceblue')
	# text(x = 20, y =  1:length(overlap), labels = names(rev(overlap)))
	# plot(1:length(colorLabels), 1:length(colorLabels), type = "n", axes = F, xlab = "", ylab = "")
	# legend("center", legend = names(colorLabels), cex = 0.3, bty = "n",
	#	ncol = length(colorLabels)/4, horiz = F, pt.cex = 1.5, x.intersp = 1.2, y.intersp =1.2,
	#	pch = 20, col = colorLabels, xpd = NA)
	# print(labelNames)
	# legend("bottom", legend = labelNames, cex = 0.3, bty = "n",
#		ncol = length(colorLabels)/2, horiz = F, pt.cex = 1.5, x.intersp = 1.2, y.intersp = 1.2,
#		pch = 20, col = colorLabels, xpd = NA)
#	dev.off()

 # compare by cellType
	cellTypes=lapply(1:length(results), function(x) gsub('Smooth muscle cells', 'Myofibroblasts', as.character(results[[x]]$cellType)))
	clusterSizes = sapply(cellTypes, function(x) length(unique(x)))
	# Remove cells not analised in all subsampled runs
	cellCols=palette$cellTypeColors[unique(unlist(cellTypes))]
	cellCols=cellCols[order(names(cellCols))]	
	colorMat  = WGCNA::labels2colors(as.data.frame(cellTypes), colorSeq = cellCols)
	colorLabels=sapply(1:length(cellTypes), function(x) tapply(colorMat[,x], cellTypes[[x]], unique))
	colorLabels=unlist(colorLabels)
	colorLabels=colorLabels[!duplicated(names(colorLabels))]
	cellTypeOrder=names(palette$cellTypeColors)[names(palette$cellTypeColors) %in% names(colorLabels)]
	colorLabels=colorLabels[cellTypeOrder]

	select=apply(colorMat[, -1], 1, function(x) ! all(x == 'grey'))
	clusterOrder = do.call(order, lapply(1:length(cellTypes), function(x) {
		matchInd=match(cellTypes[[x]], names(colorLabels))
		matchInd[is.na(cellTypes[[x]])] = NA
		matchInd
	}))
	pngOut = f("{outDir}/matchedCellTypes.{analysisID}.png")
	pdfOut = f("{outDir}/matchedCellTypes.{analysisID}.pdf")

	nrCells=length(clusters[[1]])
	clusterNames=sapply(1:length(clusters), function(x) {
		print(x)
	  label=paste("Subampling", x - 1)
	  if(x == 1) label = 'Reference'
		nrCells=sum(!is.na(clusters[[x]]))
		overlap = round(sum(cellTypes[[x]] == cellTypes[[1]], na.rm = T)/nrCells * 100, 2)
        paste0(label,'\n', clusterSizes[[x]], " clusters")
	})
	colnames(colorMat)=c('Ref', paste0('Run', 1:(length(clusters) - 1)))
	#if(!file.exists(pngOut)) {
		tiff(filename = tiffOut, res = 300, width = 8, height=3, units = 'in', type = 'cairo')
		par(mar = c(0.5, 0.5, 0.5, 0.5))
		print(head(clusterOrder))
		plotOrderedColors(colors = colorMat, order = clusterOrder,  rowLabels = names(colorMat))
		dev.off()
 	#}
 browseURL(pngOut)
	overlap=sapply(1:length(clusters), function(x) {
		if(x==1) return(100)
		nrCells=sum(!is.na(clusters[[x]]))
	})
	names(overlap)=gsub("; ", ' ', clusterNames)
	
	overlap=data.frame(t(sapply(1:length(clusters), function(x) {
        nrCells=sum(!is.na(clusters[[x]]))
        c(pct = sum(cellTypes[[x]] == cellTypes[[1]], na.rm = T) / nrCells * 100, total=nrCells, round=x, name =clusterNames[x])
    })))
  overlap$pct = as.numeric(as.character(overlap$pct))
	txtOut = f("{outDir}/matchedCellTypes.{analysisID}.overlap.txt")
	write.table(overlap, file =txtOut, sep = '\t', row.names = F, quote = F)

    pdf(pdfOut, width = 4, height = 3)
    g <- ggplot(overlap, aes(name, pct))
    plot <- g + geom_col(fill='#cedcefff', col='grey20') + #(stat='identity') + 
        cowplot::theme_cowplot() + 
        ylab(paste0("% cells in concordant clusters\n(", clusterSizes[[2]], 'cells)')) +
        xlab('') +
        coord_flip() +
        geom_text(aes(label=paste(round(pct, 1), '%' )), hjust=1) 
    print(plot)

	plot(1:length(colorLabels), 1:length(colorLabels), type = "n", axes = F, xlab = "", ylab = "")
	legend("center", legend = names(colorLabels), cex = 0.3, bty = "n",
       ncol = length(colorLabels)/4, horiz = F, pt.cex = 1.5, x.intersp = 1.2, y.intersp =1.2,
       pch = 20, col = colorLabels, xpd = NA)

    dev.off()
 # browseURL(pdfOut)
    # # Plot number of cells per cluster (using the same colors)
    # pdf(paste0("Cluster_stats.", analysis_id, ".pdf"))
    # colnames(clusters) = names(results)
    # clusterStats = melt(clusters)
    # # resultMrg = data.frame(do.call(rbind, lapply(results, function(x) cbind(cluster = x$cluster, uniqueID = x$uniqueID))))
    # g <- ggplot(clusterStats, aes(as.factor(value), fill = as.factor(Var2)))
    # # g <- ggplot(resultsMrg, aes(as.factor(cluster), fill = as.factor(uniqueID)))
    # plot = g + geom_bar(position = "dodge", color = "darkgrey") + 
    #   theme_classic() +
    #   xlab("Cluster") + ylab("# Cells") +
    #   theme(legend.position = "top", legend.title = element_blank()) +
    #   scale_fill_manual(values = brewer.pal(length(results), "Set1"))
    # print(plot)
    # # using matched cluster names
    # # calculate F1 score
    # precision = data.frame(t(sapply(colnames(clusters)[-1], function(col) {
    #   tp = sum(clusters[,col] == clusters[, 1], na.rm = T)
    #   p = tp / sum(!is.na(clusters[,1]))
    #   return(c(p = p, id = col))
    # })), check.names = F)
    # precision$p = as.numeric(as.character(precision$p))
    # g <- ggplot(precision, aes( as.factor(id), y = p))
    # plot = g + geom_bar(stat = "identity", 
    #                     fill = "transparent", col = "black") +
    #   theme_classic(base_size = 18) +
    #   dththeme(axis.text = element_text(color = "black")) +:
    #   xlab("") + ylab("Precision score")
    # print(plot)
    # dev.off() 
}

# Plot the content of the unassigned cells
# Plot the content of the Immune clusters


