library(plyr)
 
plot_overlaps <- function(upDfConf, clusterLabels, runID, fdr = 0.05, effect_field='overlap') {
	
	plotDir=f("{runID}_aux/{effect_field}/plots")
	if(! dir.exists(plotDir))
		dir.create(plotDir, recursive=T)
	
	
	fileOut=f('{plotDir}/kstest.RData')
  
	if(file.exists(fileOut))  {
		load(fileOut)
		return(pos)
	}
	upDfConf = upDfConf[grep('summary', rownames(upDfConf), invert = T), ]
	print(colnames(upDfConf))


	library(doParallel)
	# Number of cores: the smaller number from nr markers and nr available cores
	ncores=min(8, ncol(upDfConf))
	cat('Using ', ncores, ' number of cores\n')
	cluster <- makeCluster(ncores)
	registerDoParallel(cluster)

	pos = foreach(marker = colnames(upDfConf), .combine='c', 
					.export = c("effect_field", "upDfConf", "runID", "plotDir", ls(globalenv())), 
					.packages=c('tidyr', 'ggplot2')) %dopar% {

		pdfOut=f("{plotDir}/{marker}_prob.pdf")
		rdOut=f("{runID}_aux/{effect_field}/{marker}.RData")
	
		if(file.exists(rdOut))	{
			load(rdOut)
			return(positive)
		}	
		positive=vector(mode='list')
		flt = data.frame(cbind(marker = upDfConf[[marker]],
				cluster = gsub(f("^(.*).{effect_field}.*$"), "\\1", rownames(upDfConf)),
				cluster2 = gsub(f("^.*.{effect_field}.(.*)$"), "\\1", rownames(upDfConf))),
				stringsAsFactors = F,
				check.names = F)
		flt$marker = flt$marker %>% as.character %>% as.numeric
		flt$cluster = flt$cluster %>% gsub('^x', '', ., ignore.case=T)
		flt$cluster2 = flt$cluster2 %>% gsub('\\.', ' ', .)
		flt$confidence = 'High'
		flt$confidence[grepl('^Excluded', flt$cluster2)] = 'Low'
		flt$cellType = gsub("\\ [0-9]+", "", flt$cluster)
		nrLowConfClusters = flt$cluster[grepl('^Excluded', flt$cluster)] %>%
			unique %>% length
		
		pdf(pdfOut, width=6, height = 6)
		for(clusterName in unique(flt$cluster)) {
			# if spceical characters in the name, which will be used as RegEx
			clusterNameRegEx=clusterName %>% 
				gsub('\\+', '\\\\+', .) %>% 
				gsub('\\.', '\\\\.', .)
			subsetIndices=gsub('^X', '', rownames(upDfConf)) %>%
				grep(f('^{clusterNameRegEx}\\.{effect_field}.'), .)
			background=setdiff(1:nrow(upDfConf), subsetIndices)
			
			if(nrLowConfClusters > 1) {
				# if stratification is performed
				if(grepl('^Excluded', clusterName))	{
					subsetIndices=intersect(subsetIndices,
						grep(f('\\.{effect_field}.Excluded.*'), rownames(upDfConf)))
						# relative to low conf
					background=intersect(background,
						grep(f('.{effect_field}.Excluded.*'), rownames(upDfConf)))
				}
			}
			
			x = flt$marker[subsetIndices]
			y = flt$marker[background]
			cat(clusterName, length(x), length(y), '\n')
			if(length(x) > 0 & length(y) > 0) {
				vartest=ks.test(x = x, y = y, alternative = 'l', exact = F)
			} else {
				vartest=list('p.value' = 1, statistic = NA)
			}
			pos=vartest[['p.value']] <= fdr
			positive[[clusterName]][[marker]] = c(vartest[['p.value']], vartest[['statistic']])
			
			flt$curr_cluster = 'Removed'
			flt$curr_cluster[background] = 'Bg'
			flt$curr_cluster[subsetIndices] = T

			g <- ggplot(flt, aes(marker, color=curr_cluster))
			plot = g + geom_density(aes(linetype=confidence)) +
				theme_classic() +
				ggtitle(paste(clusterName, marker, pos, 
					vartest[['p.value']], '\n',
					round(vartest[['statistic']], 3)))
			print(plot)
		}
		dev.off()
		save(positive, file = rdOut)
		return(positive)
	}
	stopCluster(cluster)
	pos=unlist(pos, recursive = F)
	save(pos, file = fileOut)
	return(pos)
}

plot_heatmap <- function(dfExp, clusters,  runID, labels, plotDir = "plots", plotPos = F) {

	if(! dir.exists(plotDir))
		dir.create(plotDir)
	print('Plotting heatmap')
	heatmapData = f('{plotDir}/heatmap.RData')
	save(dfExp, clusters, runID, labels, file=heatmapData)
	clusterSize = table(clusters) %>% as.data.frame

    clusterSummary <- sapply(colnames(dfExp), function(x) {
      tapply(t(dfExp[[x]]), as.factor(clusters), mean, na.rm = T)
    })
	#print(clusterSummary)
	minValue = min(clusterSummary[clusterSummary > 0], na.rm = T)
	clusterSummary[clusterSummary == 0] = minValue
	clusterSummary[is.nan(clusterSummary) | is.na(clusterSummary)] = minValue 

    clusterNorm=apply(clusterSummary, 2, to_zscore)
    density_scale=pretty(c(-2, 2))
	
	labelsMatch=match(rownames(clusterNorm), labels$cluster)
	labels$majorType_rev=gsub(' - .*', '', labels$majorType)
	celltypes = labels$majorType_rev[labelsMatch]
	# cases when no positivity but found most likely celltype
	celltypes[is.na(celltypes)] = rownames(clusterNorm)[is.na(celltypes)]

	pdfOut=f("{plotDir}/positivity_maps.pdf")
	pdf(pdfOut, useDingbats = F, height = 12, width = 20)
	if(grepl(majorDir, runID) | grepl(sampledDir, runID)) {
		subsets='all'
	} else {
		subsets = c(unique(celltypes), 'all')
	}
	
	for(subset in subsets) {
			
			cat('Plotting subset of clusters: ', subset, runID, '\n')
			clusterNormSub = clusterNorm
			clusterSizeSub = clusterSize
			
			if(subset != 'all') {
				clusterNormSub = clusterNorm[celltypes == subset, ]
			}
			if(is.null(nrow(clusterNormSub)))
				next
			if(! nrow(clusterNormSub))
				next
			selectedCols = apply(clusterNormSub, 2, function(x) length(unique(x)) > 1)
			clusterNormSub = subset(clusterNormSub, select=selectedCols)
			clusterSizeSub = clusterSize[clusterSize$clusters %in% rownames(clusterNormSub), ]
			if(is.null(nrow(clusterNormSub)))
				next
			
			# apply(clusterNormSub, 2, function(x) length(unique(x))) %>% print
			rowLabels = gsub(".*:", "", colnames(clusterNormSub))
			heat = ComplexHeatmap::Heatmap(t(clusterNormSub),
				na_col = "grey90",
				col = circlize::colorRamp2(density_scale, 
						rev(brewer.pal(name = "RdBu", length(density_scale)))),
				row_title_gp = gpar(fontsize=12),
				row_title_rot=90,
				border = "grey85",
				row_labels = rowLabels,
				column_title_gp = gpar(fontsize=8),
				name = subset,
				column_title = subset,
				row_title = "Mean raw pixel intensity",
				width = unit(nrow(clusterNormSub)/10, 'in'),
				height = unit(ncol(clusterNormSub)/4, "cm"),
				clustering_method_rows = "ward.D2",
				clustering_method_columns = 'ward.D2',
				clustering_distance_rows='spearman',
				clustering_distance_columns='spearman',
				heatmap_legend_param = list(title = "Intensity\n[z-score]", 
				ncol=1, nrow=4,by_row=T,
				direction='horizontal', fontsize=8),
				row_names_gp = gpar(fontsize = 8),
				show_heatmap_legend=T,
				column_names_gp = gpar(fontsize = 8, angle = 90))
		print(heat)

		print(unit(nrow(clusterSummary)/10, 'in'))
		print(unit(ncol(clusterSummary)/5, "cm"))
    	labels$positive=gsub('pos:(.*) neg:', '\\1', labels$positive) %>% 
			gsub('up:(.*) down:', '\\1', .)
		clusterLabels=rownames(clusterNormSub)
		markers = colnames(clusterSummary)[row_order(heat)]
	    binMat = sapply(clusterLabels, function(cluster) {
    	  positive=labels$positive[labels$cluster == cluster] %>%
		  						strsplit(split = '_') %>% unlist
	      sapply(markers %in% positive, sum)
    	})
		binMat = as.data.frame(binMat)
    	colnames(binMat) = clusterLabels
			
		labelMatch = match(clusterSizeSub$clusters, labels$cluster)
    	clusterSizeSub$label = labels$positive[labelMatch]
	    clusterSizeSub$clusters = factor(clusterSizeSub$clusters, 
			levels = rev(rownames(clusterNormSub)[column_order(heat)]))
    	rowLabels = sapply(clusterLabels, 
			function(x) 
				paste0(labels$cellType[labels$cluster == x][1], 
					"(", ifelse(plotPos, "n=", c(x, ",n=")), sum(clusters == x), ")")
			)
    	rowCellTypes=sapply(clusterLabels, 
					function(x) labels$cellType[labels$cluster == x][1])
			rowCellTypes[is.na(rowCellTypes)] = 'Excluded'
			cellTypeColors=rep(unlist(palette$cellTypeColors), 2)
			names(cellTypeColors) = c(names(palette$cellTypeColors), 
				paste('Excluded', names(palette$cellTypeColors)))
			cellTypeColors = cellTypeColors[match(rowCellTypes, names(cellTypeColors))]
			if(all(is.na(cellTypeColors)))	{
				cellTypeColors=rainbow(length(rowCellTypes))
			}
			cellTypeColors[is.na(cellTypeColors)] = 'grey'
			names(cellTypeColors) = rowCellTypes

	    	row_ha = HeatmapAnnotation("# Cells"=anno_barplot(border = F,
										x = clusterSizeSub$Freq[match(clusterSizeSub$clusters, clusterLabels)],
										gp=gpar(fill='#0571B0', col='transparent',
										fontsize=8, title=expression("# Cells"))), 
										which = 'column', cellType = sapply(clusterLabels, 
											function(x) labels$cellType[labels$cluster == x][1]),
												col=list(cellType = cellTypeColors),
												gp = gpar(col = "grey50"))
			
			if(! all(binMat == 0) & ! all(binMat == 1)) {

			print('Plotting bin mat')
			print(rowLabels)
		    bin=ComplexHeatmap::Heatmap(binMat, cluster_rows = F,
									cluster_columns = F,
        	                        bottom_annotation = row_ha,
    	                            row_title_gp = gpar(fontsize = 10),
	                                column_title = NULL,
									row_title = 'Positivity',
                    	            rect_gp=gpar(col='grey85'),
                	                column_title_gp = gpar(fontsize = 10),
            	                    row_names_gp = gpar(fontsize = 7),
        	                        column_names_gp = gpar(fontsize = 8),
    	                            border = "grey85",
	                                na_col='white',
									width=unit(nrow(clusterNormSub)/10, 'in'),
					                height=unit(ncol(clusterNormSub)/5, "cm"),
									col=colorRampPalette(rev(brewer.pal(name = "RdBu", 5)[-c(2, 4, 5)]))(100),
                	                row_labels = markers, #[rowPosSel],
            	                    column_labels = rowLabels,
        	                        heatmap_legend_param =
										list(title="Expressed marker",
											 direction='horizontal', fontsize=8))
			ht_list=heat %v% bin
			draw(ht_list, annotation_legend_side = "bottom", 
				heatmap_legend_side = "right", merge_legend=T)
			print('Binary')
		} else {
				cat('Skipping binary plot for ', subset, '\n')
			}
		}
		dev.off()

		clusterSize$cellType=labels$cellType[match(clusterSize$clusters, labels$cluster)]
		stats=ddply(clusterSize, .(cellType), summarise, TotalFreq = sum(Freq))    
		write.tab(stats, f('{runID}.major_stats.txt'))
		
		cellTypeColors = rep(unlist(palette$cellTypeColors), 2)
		names(cellTypeColors) = c(names(palette$cellTypeColors), 
			paste('Excluded', names(palette$cellTypeColors)))
		cellTypeColors=cellTypeColors[names(cellTypeColors) %in% stats$cellType]
		if(! length(cellTypeColors))	{
			cellTypeColors=rainbow(length(stats$cellType))
			names(cellTypeColors) = stats$cellType
		}
	
		cellTypeColors[is.na(cellTypeColors)] = 'grey'
		missing=unique(stats$cellType[! stats$cellType %in% names(cellTypeColors)])
		cat("WARNING: No color annotation for the following celltypes: ", missing, 
			". Amend in conf/celltype_colors.json.\n",
			file = f("{runID}.log"), append = T)
			
		cellOrder = order(match(stats$cellType, names(cellTypeColors)))
		stats$cellType = factor(stats$cellType, levels = unique(stats$cellType[cellOrder]))
		pdfOut = f("{plotDir}/cell_types_pie_chart.pdf")
		pdf(pdfOut, useDingbats = F, height = 5, width = 5)
	    g <- ggplot(stats, aes(x="", y = TotalFreq, fill = cellType))
    	plot = g + geom_bar(stat="identity", position =  'fill', color = 'black', alpha = 0.8) + 
	      cowplot::theme_cowplot() +
		  coord_polar(theta = "y", start = 0, direction = -1) + 
    	  theme(axis.line = element_blank(), 
        	    strip.background = element_blank(),
            	axis.text = element_blank(),
	            axis.ticks = element_blank()) +
    	  xlab("") + ylab("") +
	      scale_y_continuous(labels=function(x) {
    	    paste0(x * 100, "%")
	      }) +  # facet_wrap(. ~ confidence + cellassign_cluster, nrow= 2) +
		scale_fill_manual(values = cellTypeColors)
	    print(plot)
    	dev.off()
	
	print('Plotted heatmap')
    clusters_order = rownames(clusterSummary)[row_order(heat)]
    markers_order = colnames(clusterSummary)[column_order(heat)]
    return(list(cluster = clusters_order, marker = markers_order))
}

plot_binary <- function(labels, runID, clusters_order=NULL, markers_order=NULL) {

  labels$positive=gsub('pos:(.*) neg:', '\\1', labels$positive)
  if(!is.null(clusters_order)) {
    labels=labels[order(match(labels$cluster, clusters_order)), ]
  }
  markers=sapply(labels$positive,  function(x) strsplit(x, split='_')[[1]]) %>% unlist %>% unique
  if(!is.null(markers_order)) {
    markers=markers[order(match(markers, markers_order))]
  }
  binMat=sapply(labels$cluster, function(cluster) {
      positive=labels$positive[labels$cluster == cluster] %>% 
	  	strsplit(split = '_') %>% unlist
      sapply(markers %in% positive, sum)
  })

  pdfOut=paste0(runID, ".binmat.pdf")
  pdf(pdfOut, height=8)
  pheatmap(t(binMat), cluster_cols = F, cluster_rows = F,
           fontsize_row = 6, fontsize_col = 6,
           border_color = "grey85",
           # gaps_row = 1:length(labels_row),
           cellheight =  6,
           cellwidth = 6, na_col = 'white',
           labels_row = labels$cellType,
           labels_col = markers,
           angle_col = 90,
           color = colorRampPalette(rev(c('#1356bbff', "#0066FFFF", "aliceblue")))(100))
  dev.off()
}

plot_expression <- function(dfExp, clusters, clusterNames, p, magnitude=NULL, plotDir = 'plots') {
	
	if(! dir.exists(plotDir))
		dir.create(plotDir)
	
	dfExp = subset(dfExp, select = which(! colnames(dfExp) %in% c('imagename', 'ObjectNumber')))
	if(! is.null(magnitude))
		dfExp = to_magnitude(dfExp, magnitude)
	
	if(! all(clusters %in% clusterNames$cluster)) {
		print('WARNING: Not all clusters have been annotated based on positivity')
		dfExp = dfExp[clusters %in% clusterNames$cluster, ]
		clusters = clusters[clusters %in% clusterNames$cluster]
	}
	positive = clusterNames$positive[match(clusters, clusterNames$cluster)]
	clusterCount = length(unique(paste(clusters, positive)))
	cat('Number of clusters', clusterCount, '\n')

	dfExpMag = dfExp + 1
	
	medianVal = apply(dfExpMag, 2, median, na.rm = T)
	cat("Median intensity", medianVal, '\n')
	
	for(cellMarker in colnames(dfExpMag)) {
		cat(".")
		values = dfExpMag[[cellMarker]]
		clNamesOrder = names(table(clusters))

		cols = get_marker_frequency(data = clusterNames, marker=cellMarker, column = 'positive')
		cols = cols[match(clNamesOrder, clusterNames$cluster)]

		cols = sapply(cols, function(col) ifelse(grepl('\\+', col), 'red', 'transparent'))
		boxplot(values ~ clusters, col=cols, las=.5, cex.lab.x = 0.5,
			cex.main=.5, varwidth=F, outline=F,  log = "y", xaxt='n',
			xlab="", ylab = "", ylim=c(min(dfExpMag, na.rm = T), max(dfExpMag, na.rm = T)),
			main = with(clusterNames, paste("Cell marker:", cellMarker, "\n"),
									  paste0("n=", length(indices), " cells")),)
		axis(side=1, 1:length(clNamesOrder), clNamesOrder, las = 2, srt=45)
		abline(h=medianVal[[cellMarker]], lty=2, col='red')
  }
  cat("\n")
  dev.off()
  print("Plotted raw intensity per marker")
  
  
	
	pdf(file=f("{plotDir}/raw_intensities.per_cluster.pdf"), 
		width =  7 + round(clusterCount / 60),
		height=5 + round(clusterCount / 60))
	par(mfrow=c(2, 2))
	for(.row in 1:nrow(clusterNames)) {
		cat(".")
		indices=which(clusters == clusterNames$cluster[.row] & positive == clusterNames$positive[.row])
		inFlt = dfExpMag[indices, ]
		boxplot(inFlt, outline=FALSE,
			main = with(clusterNames, f("cluster:{cluster[.row]}\n",
							gsub("_", "+", gsub('pos:(.*)neg:', "\\1", positive[.row])), "\n")),
			ylab=paste0("n=", length(indices), " cells"),  cex.main=.5, las=1,
			ylim=c(min(dfExpMag, na.rm = T), max(dfExpMag, na.rm = T)), cex.lab.x=.5, xaxt='n', log='y')
		points(medianVal, pch = 23, col='red')
		axis(side=1, 1:ncol(inFlt), colnames(inFlt), las = 2, srt=45)
	}
	cat("\n")
	dev.off()
    print("Plotted raw intensity per cluster")
	
	pdf(file=f("{plotDir}/raw_intensities.per_marker.pdf"), width=7 + round(clusterCount / 60),
		height=4 + round(clusterCount / 60))
	# par(mfrow=c(2, 2))

}
