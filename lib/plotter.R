library(plyr)
 
plot_overlaps <- function(upDfConf, clusterLabels, runID, fdr = 0.05, effect_field='overlap') {
  
  fileOut=f('{runID}.kstest.RData')
  
	if(file.exists(fileOut))  {
    	load(fileOut)
    	return(positive)
	}

   positive=vector(mode='list')
   for(marker in colnames(upDfConf))  {

     flt=data.frame(cbind(marker=upDfConf[[marker]],
                         cluster=gsub(f("^(.*).{effect_field}.*$"), "\\1", rownames(upDfConf)),
                         cluster2=gsub(f("^.*.{effect_field}.(.*)$"), "\\1", rownames(upDfConf))),
                   stringsAsFactors = F, check.names = F)
    flt$marker=as.numeric(as.character(flt$marker))
	flt$cluster=gsub('^x', '', flt$cluster, ignore.case=T)
    flt$cluster2=gsub('\\.', ' ', flt$cluster2)
    flt$confidence='High'
    flt$confidence[grepl('^Excluded', flt$cluster2)] = 'Low'
    flt$cellType=gsub("\\ [0-9]+", "", flt$cluster)
    flt$cellType2=gsub(" [0-9]+", "", flt$cluster2)

	plotDir=f("{runID}_{effect_field}")
	if(! dir.exists(plotDir)) 
		dir.create(plotDir)
    pdf(f("{plotDir}/{marker}_prob.pdf"), width=6, height = 6)
    for(clusterName in unique(flt$cluster)) {
        subsetIndices=gsub('^X', '', rownames(upDfConf)) %>% 
			grep(f('^{clusterName}\\.{effect_field}.'), .)
        # Exclude comparisons to cells from the cluster 'Excluded'
        subsetIndices=setdiff(subsetIndices, grep(f('\\.{effect_field}.Excluded.*'), rownames(upDfConf)))
        if(clusterName != 'Excluded')
          subsetIndices=setdiff(subsetIndices, grep(f('^Excluded\\.{effect_field}'), rownames(upDfConf)))
        # do not compare low confidence to high confidence clusters from the same cell type (p=0 skew)
        if(length(grep('^Excluded', clusterName))) {
          cellType=gsub('^Excluded.(.*)\\.[0-9]+$', '\\1', clusterName)
          subsetIndices=setdiff(subsetIndices, grep(f('.{effect_field}.{cellType}\\.[0-9]+$'), rownames(upDfConf)))
        }
      x=flt$marker[subsetIndices]; y = flt$marker[-subsetIndices]
      if(length(x) && length(y)) {
		vartest=ks.test(x = x, y = y, alternative = 'l', exact=F)
      } else {
      	vartest=list('p.value' = 1, statistic=NA)
      }
      pos=vartest[['p.value']] <= fdr
      positive[[clusterName]][[marker]] = c(vartest[['p.value']], vartest[['statistic']])
      
      flt$curr_cluster=F
      flt$curr_cluster[subsetIndices]=T
	  
      g <- ggplot(flt, aes(marker, color=curr_cluster))
      plot = g + geom_density(aes(linetype=confidence)) +
        #stat_ecdf(aes(linetype=confidence)) +
        theme_classic() +
        ggtitle(paste(clusterName, marker, pos, vartest[['p.value']], '\n',
                      # criticalVal,
                      round(vartest[['statistic']], 3))) +
        guides(color=F) 
        # scale_color_manual(values=cellTypeColors)
      print(plot)
    }
    dev.off()
  }
  save(positive, file = fileOut)
  # assigned=do.call(rbind, positive)
  return(positive)
}

## FINISH
plot_overlap_distribution <- function(upDfConf, runID, positive, criticalValue, clusterLabels) {

  for(cluster in clusterLabels) {
    pdf(f("{runID}.{cluster}.AUC.pdf"), width=15)

    cluster=gsub(' ', '.', cluster)
    subsetIndices=grep(f('^{cluster}\\.AUC.'), gsub('^X', '', rownames(upDfConf)))
    
    # Exclude comparisons to cells from the cluster 'Excluded'
    subsetIndices=setdiff(subsetIndices, grep('\\.AUC.Excluded.*', rownames(upDfConf)))
    if(cluster != 'Excluded')
      subsetIndices=setdiff(subsetIndices, grep('^Excluded\\.AUC', rownames(upDfConf)))
    
    # do not compare low confidence to high confidence clusters from the same cell type (p=0 skew)
    if(length(grep('^Excluded', cluster)))  {
      cellType=gsub('^Excluded.(.*)\\.[0-9]+$', '\\1', cluster)
      subsetIndices=setdiff(subsetIndices, grep(f('.AUC.{cellType}\\.[0-9]+$'), rownames(upDfConf)))
    }
    
    selection=sapply(colnames(upDfConf), function(marker) {

      overlap=upDfConf[[marker]]
      vartest=ks.test(x = overlap[subsetIndices], y = overlap, alternative = 'g')
      positive = vartest$p.value > fdr & vartest$statistic <= criticalValue
      title=paste(marker, cluster, positive, vartest$p.value, '\n', vartest$statistic)
      
      plot(density(overlap), main = title)
      if(mean(overlap[subsetIndices]) > mean(overlap) + sd(overlap)) {
        lines(density(overlap[subsetIndices]), breaks = 500, main = title, 
			xlim = c(0, 1), col='yellow', add = T)
      } else {
        lines(density(overlap[subsetIndices]), breaks = 500, main = title, col ='orange', add = T)
      }
      threshold=qnorm(0.05, mean = 0.5, sd=sd(overlap), lower.tail = F)
      abline(v=threshold, col='red', lty=2, lwd=2)
      abline(v=mean(overlap) + sd(overlap), col='blue', lty=2, lwd=2)
      abline(v=mean(overlap[subsetIndices]), col='green', lwd=4)
      
      positive
    })
    dev.off()
   # cat(cluster, paste1(colnames(upDfConf)[selection]), '\n')
    paste0("pos:", paste1(colnames(upDfConf)[selection]), " neg:")
  }
}

plot_heatmap <- function(dfExp, clusters,  runID, labels) {

    print('Plotting heatmap')	
	heatmapData=f('{runID}.heatmap.RData')
	save(dfExp, clusters, runID, labels, file=heatmapData)

    clusterSize=table(clusters) %>% as.data.frame

    clusterSummary<-sapply(colnames(dfExp), function(x) {
      tapply(t(dfExp[[x]]), as.factor(clusters), median)
    })
    
   # clusterSummary[clusterSummary == 0]=min(clusterSummary[clusterSummary>0])
	#clusterSummary=log2(clusterSummary)

    # save(clusterSummary, labels, runID, clusterSize, file=f('{runID}.clusterSummary.RData'))
    # load("/Users/angelom/Documents/presentations/20_12_09_TME_meeting/major_mcsa/30_0.67_none_FALSE.clusterSummary.RData")
    
    clusterNorm=apply(clusterSummary, 2, to_zscore)
    density_scale=pretty(c(-4, 4))
	labelsMatch=match(rownames(clusterNorm), labels$cluster)

	celltypes = labels$majorType[labelsMatch]
	# cases when no positivity but found most likely celltype
	celltypes[is.na(celltypes)] = rownames(clusterNorm)[is.na(celltypes)]

    pdfOut=paste0(runID, ".maps.pdf")
	pdf(pdfOut, useDingbats = F, height = 12, width = 10)
	if(grepl(majorDir, runID) | grepl(sampledDir, runID)) {
		subsets='all'
	} else {
		subsets = unique(celltypes)
	}

	for(subset in subsets) {
		
		cat('Plotting subset of clusters: ', subset, '\n')
		clusterNormSub = clusterNorm
		clusterSizeSub=clusterSize

		if(subset != 'all') {
			clusterNormSub = clusterNorm[celltypes == subset, ]
			clusterSizeSub = clusterSize[clusterSize$clusters %in% rownames(clusterNormSub), ]
		}
		if(is.null(nrow(clusterNormSub)))
			next

    	heat=ComplexHeatmap::Heatmap(t(clusterNormSub),
                                 na_col = "grey90",
								 col = circlize::colorRamp2(density_scale, rev(brewer.pal(name = "RdBu", length(density_scale)))),
                                 #col=colorRampPalette(colors = rev(brewer.pal(name = "RdBu", 3)))(255),
                                 # col=colorRampPalette(colors = rev(brewer.pal(name = "Accent", 3)))(255),
                                 row_title_gp = gpar(fontsize=12),
                                 row_title_rot=0,
                                 border = "grey85",
                                 row_labels = gsub(".*:", "", colnames(clusterSummary)),
                                 column_title_gp = gpar(fontsize=8),
								 name = subset,	
								 column_title = subset,
                                 width=unit(nrow(clusterNormSub)/10, 'in'),
                                 height=unit(ncol(clusterNormSub)/5, "cm"),
                                 clustering_method_rows = "ward.D2",
                                 clustering_method_columns = 'ward.D2',
								 clustering_distance_rows='maximum',
								 clustering_distance_columns='maximum',
                                 heatmap_legend_param = list(title = "Intensity\n[z-score]", 
                                                             ncol=1, nrow=4,by_row=T,
                                                             direction='horizontal', fontsize=8),
                                 row_names_gp = gpar(fontsize = 8),
                                 show_heatmap_legend=T,
                                 # clustering_distance_rows='pearson',
                                 column_names_gp = gpar(fontsize = 8, angle=90))
	   # print(heat)
		print(unit(nrow(clusterSummary)/10, 'in'))
		print(unit(ncol(clusterSummary)/5, "cm"))
    	labels$positive=gsub('pos:(.*) neg:', '\\1', labels$positive) %>% 
			gsub('up:(.*) down:', '\\1', .)
		markers=sapply(labels$positive,  function(x) strsplit(x, split='_')[[1]]) %>% 
																unlist %>% unique
		clusterLabels=rownames(clusterNormSub)
		markers=colnames(clusterSummary)[row_order(heat)]
	    binMat=sapply(clusterLabels, function(cluster) {
    	  positive=labels$positive[labels$cluster == cluster] %>% 	
		  						strsplit(split = '_') %>% unlist
	      sapply(markers %in% positive, sum)
    	})
		binMat=as.data.frame(binMat)
    	colnames(binMat) = clusterLabels
	    rownames(binMat) = markers
	 
    	clusterSizeSub$label=labels$positive[match(clusterSizeSub$clusters, labels$cluster)]
	    clusterSizeSub$clusters=factor(clusterSizeSub$clusters, 
			levels=rev(rownames(clusterNormSub)[column_order(heat)]))
    	rowLabels=sapply(clusterLabels, function(x) paste0(labels$cellType[labels$cluster == x][1], 
															" (", x, 
															',n=', sum(clusters == x), ")"))
	
	    row_ha = HeatmapAnnotation("# Cells"=anno_barplot(border = F,
			x = clusterSizeSub$Freq[match(clusterSizeSub$clusters, clusterLabels)],
					gp=gpar(fill='#0571B0', col='transparent',
					fontsize=8, title=expression("# Cells"))), 
						which = 'column', cellType = sapply(clusterLabels, 
								function(x) labels$cellType[labels$cluster == x][1]),
							col=list(cellType = unlist(palette$cellTypeColors)),
						 gp = gpar(col = "grey50"))
	
    	if(!all(binMat == 0) & !all(binMat == 1)) {

			print('Plotting bin mat')
			print(binMat)
		    bin=ComplexHeatmap::Heatmap(binMat, cluster_rows = F,
									cluster_columns = F,
        	                        bottom_annotation = row_ha,
    	                            row_title_gp = gpar(fontsize = 10),
	                                column_title = NULL,
                    	            rect_gp=gpar(col='grey85'),
                	                column_title_gp = gpar(fontsize = 10),
            	                    row_names_gp = gpar(fontsize = 7),
        	                        column_names_gp = gpar(fontsize = 8),
    	                            border = "grey85",
	                                na_col='white',
                                	width=unit(nrow(clusterNormSub)/10, 'in'),
                            	    height=unit(ncol(clusterNormSub)/5, "cm"),
									col=colorRampPalette(rev(brewer.pal(name = "RdBu", 5)[-c(2, 4, 5)]))(100),
                    	            # col=colorRampPalette(rev(c('#1356bbff', "#0066FFFF", "aliceblue")))(100),
                	                row_labels = markers,
            	                    column_labels = rowLabels,
        	                        heatmap_legend_param = list(title="Expressed marker", direction='horizontal', fontsize=8))
	 			if(! pars$stratify & ! pars$subset == majorDir) {
	    			ht_list=heat %v%  bin
			    	draw(ht_list, annotation_legend_side = "bottom", heatmap_legend_side = "bottom", merge_legend=T)
				}
				print('Binary')
			} else {
				cat('Skipping binary plot for ', subset, '\n')
			}
		}

		clusterSize$cellType=labels$cellType[match(clusterSize$clusters, labels$cluster)]
		stats=ddply(clusterSize, .(cellType), summarise, TotalFreq = sum(Freq))    
		write.tab(stats, f('{runID}.major_stats.txt'))
		
		cellTypeColors=rep(unlist(palette$cellTypeColors), 2)
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
			". Amend in conf/celltype_colors.json\n", 
			file = f("{runID}.log"), append = T)
			print('CellType colors')
		cellOrder = order(match(stats$cellType, names(cellTypeColors)))
		stats$cellType = factor(stats$cellType, levels = unique(stats$cellType[cellOrder]))

	    g <- ggplot(stats, aes(x="", y = TotalFreq, fill = cellType))
    	plot = g + geom_bar(stat="identity", position =  'fill', color = 'grey40', alpha = 0.8) + 
	      cowplot::theme_cowplot() +coord_polar(theta="y", start=0, direction=-1) + 
    	  theme(axis.line=element_blank(), 
        	    strip.background = element_blank(),
            	axis.text=element_blank(),
	            axis.ticks=element_blank()) +
    	  xlab("") + ylab("") +
	      scale_y_continuous(labels=function(x) {
    	    paste0(x * 100, "%")
	      }) +  # facet_wrap(. ~ confidence + cellassign_cluster, nrow= 2) +
		scale_fill_manual(values=cellTypeColors)
	    print(plot)
    	dev.off()
	
	print('Plotted heatmap')
    clusters_order=rownames(clusterSummary)[row_order(heat)]
    markers_order=colnames(clusterSummary)[column_order(heat)]
    return(list(cluster=clusters_order, marker=markers_order))
}

plot_cluster_size <- function(clusters, labels, runID, clusters_order=NULL) {
 stats=table(clusters) %>% as.data.frame
 stats$label=labels$positive[match(stats$clusters, labels$cluster)]
 stats$label=gsub('pos:(.*) neg:', '\\1', stats$label)
 if(is.null(clusters_order)) {
   clusters_order=unique(stats$clusters[order(stats$Freq)])
 }
 stats$clusters=factor(stats$clusters, levels=rev(clusters_order))
 
 pdfOut=paste0(runID, ".cluster_size.pdf")
 pdf(pdfOut, height = 8, width=2.5)
 g <- ggplot(stats, aes(clusters, Freq))
 plot=g + geom_bar(stat='identity', fill='#1356bbff') + cowplot::theme_cowplot() +
   xlab('Clusters') + ylab('Frequency') + coord_flip()
 print(plot)
 dev.off()
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
      positive=labels$positive[labels$cluster == cluster] %>% strsplit(split = '_') %>% unlist
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

plot_expression <- function(dfExp, clusters, clusterNames, p, magnitude=NULL) {

  if(! is.null(magnitude))
    dfExp=to_magnitude(dfExp, magnitude)
  # TODO: e.g. no positivity but with majorType doesn't end up in clusterNames
  if(! all(clusters %in% clusterNames$cluster)) {
	  dfExp = dfExp[clusters %in% clusterNames$cluster, ]
  	  clusters = clusters[clusters %in% clusterNames$cluster]
 }
  clusterCount=length(unique(clusters))
  clusterNames$up = gsub("pos:(.*) neg:.*", "\\1", clusterNames$cellType_rev)
  positive=clusterNames$up[match(clusters, clusterNames$positive)]
  up=paste(clusters, positive)
  dfExpMag=dfExp + .01
  medianVal=apply(dfExpMag, 2, median)
  pdf(file=paste0(p$run_id, ".inMat.per_cluster.pdf"), width =  7 + round(clusterCount / 60),
		height=5 + round(clusterCount / 60))
  par(mfrow=c(2, 2))
  for(.row in 1:nrow(clusterNames)) {
    cat(".")
    indices=which(clusters == clusterNames$cluster[.row])
    inFlt=dfExpMag[indices, ]
    boxplot(inFlt, outline=FALSE, 
            main = with(clusterNames, f("cluster:{cluster[.row]}\n",
									gsub("_", "+", gsub('pos:(.*)neg:', "\\1", positive[.row])), "\n")),
              ylab=paste0("n=", length(indices), " cells"),  cex.main=1, las=1, 
            ylim=c(min(dfExpMag), max(dfExpMag)), cex.lab.x=1, xaxt='n', log='y')
    points(medianVal, pch = 23, col='red')
    axis(side=1, 1:ncol(inFlt), colnames(inFlt), las = 2, srt=90)
    # text(x=1:ncol(inFlt), y = par("usr")[3] - 0.45, labels=colnames(inFlt),
         # xpd=NA, srt=90, adj = 0.965, cex=1)
    # inPct=apply(dfExpMag, 2, to_percentile)[indices, ]
    # boxplot(inPct, main = with(clusterNames,
    #                            paste("cluster:", cluster[.row], "\n", 
    #                                  gsub(" down", "\ndown", markers[.row]))),
    #         xlab=paste0("n=", length(indices), " cells"), cex.main=1, las=2,
    #         ylim=c(0, 1),  outline=F, cex.lab=0.1 )
  }
  cat("\n")
  dev.off()
  
  pdf(file=paste0(p$run_id, ".inMat.per_cell.pdf"), width=7 + round(clusterCount / 60),
		height=4 + round(clusterCount / 60))
  # par(mfrow=c(2, 2))
  for(cellMarker in colnames(dfExpMag)) {
    cat(".")
    values=dfExpMag[[cellMarker]]
    clNamesOrder=names(table(clusters))

    cols=get_marker_frequency(data = clusterNames, marker=cellMarker, column = 'positive')
		cols=cols[match(clusterNames$cluster, clNamesOrder)]

    cols=sapply(cols, function(col) ifelse(grepl('\\+', col), 'red', 'transparent'))
    boxplot(values ~ clusters, col=cols, las=1, cex.lab = 0.4, 
            cex.main=2, varwidth=F, outline=F,  log = "y", xaxt='n',
						 xlab="", ylab = "", ylim=c(min(dfExpMag), max(dfExpMag)), 
            main = with(clusterNames, paste("Cell marker:", cellMarker, "\n"), 
                        paste0("n=", length(indices), " cells")),)
    axis(side=1, 1:length(clNamesOrder), clNamesOrder, las = 2, srt=90)
    abline(h=medianVal[[cellMarker]], lty=2, col='red')
    # boxplot(to_percentile(dfExp[[cellMarker]]) ~ up,
    #         main = with(clusterNames, 
    #                     paste("Cell marker:", cellMarker, "\n"),
    #                     paste0("n=", length(indices), " cells")), las=2,
    #         xlab="", ylim=c(0, 1), cex.lab = 0.1, 
    #         cex.main=.5, log = "y", varwidth=T, outline=F)
  }
  cat("\n")
  dev.off()
}
