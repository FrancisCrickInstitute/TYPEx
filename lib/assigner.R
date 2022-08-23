# Generic function

summarise_output <- function(inData, method, pars, runID, runOutput, columnNames, nfDir, colors=NULL) {
  
		ids=with(inData, paste(ObjectNumber, basename(file)))

		pars[[method]]$run_id=runID
		pars[[method]]$markers=pars$markers
		pars[[method]]$magnitude=pars$magnitude

		if(pars$markers %in% c("all", "phenotypic", "functional")) {
			markersList=marker_gene_list[[pars$defaultMarkersList]]
		} else {
			markersList=marker_gene_list[[pars$markers]]
		}

		print("Cell type assignment")
		clusterNameFile=f("{runID}.clusterNames.txt")
		clusterNames=assign_clusters(
			dfExp = inData[, ..columnNames], 
			clusters=runOutput,
			run=pars$run,
			runID=runID, 
			panel=pars$panel,
			magnitude=pars$magnitude,
			fromList=ifelse(method=="cellassign", pars$markers, NA),
			ilastik=pars$ilastik)
		print(clusterNames)
		print("Clusters assigned")
		clusterNames$cluster = gsub('\\.', ' ', clusterNames$cluster)
		majorTypes=NULL
		if(method != 'cellassign' & !all(is.null(colors))) {
			if(!all(colors=='black')) majorTypes=gsub(" [0-9]+", "", colors)
		}
		clusterNames$positive=sapply(clusterNames$positive, toString)
		positive=clusterNames$positive[match(runOutput, clusterNames$cluster)]
		positive[is.na(positive)] = 'pos: neg:'
		print("Getting tissue category")
		print(pars$tissAreaDir)
		tissue_categs=get_tissue_category(
			cellIDs = gsub(".txt", "", ids),
			run = pars$run, 
			panel = pars$panel,
			study=pars$study,
			cohort=pars$cohort,
			tissAreaDir=pars$tissAreaDir,
			class="regional")
		print(head(tissue_categs))
		print('Major types')
		majorTypes=assign_celltype(
			names = positive, 
			majorTypes=majorTypes,
			markers = marker_gene_list[[pars$ref_markers]],
			major = T, 
			region=tissue_categs$Tumour)
		print(table(paste(positive, majorTypes)))
		print("Getting marker expression")
		meanFeature=get_markers_expression(
			inData[, ..columnNames],
			clusters=runOutput,
            clusterNames = clusterNames, 
			fun="mean", 
            magnitude = pars$magnitude)

		print("Getting spatial coordinates")
		spatData=load_files(file.path(nfDir, "LocationCenter"), run=pars$run)
		spatMatch=match(with(inData, paste(ObjectNumber, basename(file))),
				  with(spatData, paste(ObjectNumber, basename(file))))

		print('Getting area')
		areaData=load_files(file.path(nfDir, "AreaShape"), run=pars$run)
		areaMatch=match(ids, with(areaData, paste(ObjectNumber, basename(file))))

		print("Getting probabilities")
		probabilities=get_probabilities(cellIDs = ids, clusters=runOutput,
								  probFile = f("{runID}.probs.txt"))

		if(method == "cellassign")  {
			celltypes=runOutput
			celltypes=gsub('CD([48])$', 'CD\\1 T cells', celltypes)
			clusterNames$cluster=gsub('CD8$', 'CD8 T cells', clusterNames$cluster)
			clusterNames$major=clusterNames$cluster
		} else {
			print('Cell type assignment')
			celltypes=assign_celltype(names = positive, markers = markersList, 
							  majorTypes = majorTypes,
							  region=tissue_categs$Tumour)
			if(pars$subset == subtypesDir) {
				clusterNames$major=gsub("[0-9]+$|Excluded$", "", clusterNames$cluster)
				clusterNames$major=gsub(' $', '', clusterNames$major)
			} else {
				clusterNames$major=NULL
			}
		}
		clusterNames$majorType=assign_celltype(
			names = clusterNames$positive,
			majorTypes=clusterNames$major,
			markers = marker_gene_list[[pars$ref_markers]],
			major = T)
		clusterNames$cellType=assign_celltype(
			names = clusterNames$cellType, 
			markers = markersList,
			majorTypes = clusterNames$majorType)
		
	  write.tab(clusterNames, file=clusterNameFile)
	  print('Marker expression heatmap')
	  nrClusters=length(unique(runOutput))
	  if(nrClusters > 300) {
		  sel=1:nrow(inData)/2
		  heatmapOrder=plot_heatmap(dfExp = inData[sel, ..columnNames],
    	                        clusters = runOutput[sel], runID = runID,
								labels=clusterNames)
		  sel=(nrow(inData)/2 + 1):nrow(inData)
		  heatmapOrder=plot_heatmap(dfExp = inData[sel, ..columnNames],
                                clusters = runOutput[sel], runID = runID,
                                labels=clusterNames)

	  } else {
		heatmapOrder=plot_heatmap(dfExp = inData[, ..columnNames],
                                clusters = runOutput, runID = runID,
                                labels=clusterNames)
		
	 }

		print('Reviewing cell assignments')
		print(pars$cellAssignFile)
		clusterNames$cellType_rev=review_cellType_by_major(cellTypes = clusterNames$cellType,
												 majorTypes = clusterNames$majorType,
												 positivity = clusterNames$positive, 
												 panel=pars$panel, study=pars$study, cohort=pars$cohort,
												 cellAssignFile=pars$cellAssignFile)
		clusterNames$majorType_rev=review_major_by_cellType(cellTypes = clusterNames$cellType,
												  majorTypes = clusterNames$majorType,
												  positivity = clusterNames$positive,
												  panel=pars$panel, study=pars$study,
												  cohort=pars$cohort, cellAssignFile=pars$cellAssignFile)


	if(pars$panel == 'p2' & pars$subset == 'subtypes' & pars$method == 'FastPG') {
	  celltypes[runOutput == 'Epithelial cells 26'] = 'Epithelial cells'
	  majorTypes[runOutput == 'Epithelial cells 26'] = 'Epithelial cells'
		
	  sel = runOutput %in% c('Macrophages 21', 'Macrophages 20', 
			'Excluded Macrophages 12', 'Excluded Macrophages 24')
	  celltypes[sel] = 'Macrophages'
	  majorTypes[sel] = 'Macrophages'


	  sel = runOutput == 'Epithelial cells 40'
	  celltypes[sel] = 'Epithelial cells'
      majorTypes[sel] = 'Epithelial cells'
}
	  #print("Plotting expression")
	  #plot_expression(dfExp = inData[, ..columnNames],
      #            clusters = runOutput,
      #            clusterNames, pars[[method]],
      #            magnitude = pars$magnitude)

		print("Wrapping up in data frame")
		outDF=data.frame(imagename       =gsub(".txt", "", basename(inData$file)),
					   object            =sapply(inData$ObjectNumber, toString),
					   panel             =pars$panel,
					   cluster           =sapply(runOutput, toString),
					   cellType          =celltypes,
					   runID             =gsub(paste0(pars$inDir, "/?"), "", runID),
					   majorType         =majorTypes,
					   positive          =positive,
					   centerX           =spatData$LocationCenter_X[spatMatch],
					   centerY           =spatData$LocationCenter_Y[spatMatch],
					   meanIntensity     =meanFeature,

					   area              =areaData$AreaShape_Area[areaMatch],
					   probability       =probabilities,
					   tissue_categs
		)
		write.tab(outDF, file=f("{runID}.clusters.txt"))
		fst::write_fst(outDF, path=f("{runID}.clusters.fst"), compress = 75)
		return(outDF)
}
  
assign_clusters <- function(dfExp, clusters, panel, run, runID,
			fdr=0.05, log2=F, magnitude=10**5,
			fromList=NA, ilastik=T)  {
  # panel="P1"; run="B"; fdr=0.05; log2=F; magnitude=10**6; runID="~/labwd/analyses/typing/test"
  effect_field='AUC'
  scranPkg=packageVersion("scran")
  if(scranPkg == "1.12.1") effect_field='overlap'
  if(length(unique(clusters)) == 1)
    return(data.frame(cluster=rep(NA, length(clusters)), name=rep(NA, length(clusters))))
  if(!is.null(magnitude)) dfExp=to_magnitude(dfExp, magnitude)
  if(log2) dfExp=log2(dfExp + .1)
  clusterSize=table(clusters)
  clusterLabels=names(clusterSize)
  
  # Not required for the latest version of the scran package
  if(!file.exists(f('{runID}.pairwise.RData'))) {

	print('Runing pairwise')
    # Get probabilities of overlap
    wilcoxTestUp=scran::pairwiseWilcox(t(dfExp), as.factor(clusters), direction="up")
	#save(wilcoxTestUp, file = f("${runID}.wtu.RData"))
	print(packageVersion("scran"))
    markersUp=scran::combineMarkers(wilcoxTestUp$statistics, 
		wilcoxTestUp$pairs, effect.field = effect_field) # previously overlap
	print('rownames')
    rowNames=rownames(markersUp[[1]])
    # Data frame with probabilities for each pairwise comparison
    upDf=do.call(cbind, lapply(markersUp, function(x) x[rowNames, ]))
    upDf=upDf[, grep(effect_field, colnames(upDf))]
    # Summarise per marker
    upStats=apply(upDf, 1, mean) + apply(upDf, 1, sd)
    upDf=as.data.frame(t(as.data.frame(upDf)))
    save(upDf, upStats, file = f('{runID}.pairwise.RData'))
  } else  {
    load(f('{runID}.pairwise.RData'))
  }
  
  if(! is.na(fromList)) {
	print('Assigning')
    positive=sapply(clusterLabels, function(cluster) {
	  cluster=gsub('CD([48])$', 'CD\\1 T cells', cluster)
      upRows=marker_gene_list[[fromList]][[cluster]]
      if(ilastik) upRows=intersect(upRows, marker_gene_list$mcsa)
      paste0("up:", paste1(upRows), " ", "down:")
    })
	return(data.frame(cluster=names(positive), positive=positive))
  }
	print('Assigning')
    upDfConf=upDf[grep(f('.{effect_field}.Excluded.*'), rownames(upDf), invert = T), ]
	upDfConf=upDfConf[grep('summary', rownames(upDfConf), invert = T),]
    # New way
    assigned=plot_overlaps(upDfConf, clusterLabels, runID, effect_field=effect_field)
	if(pars$subset %in% c(majorDir, sampledDir) | !pars$stratify ) {

		positive=determine_threshold(assigned, clusterSize, runID, confidence='all', greater=F)
		return(data.frame(rbind(cbind(cluster=names(positive), positive=positive))))
	} else {
			positiveLow=determine_threshold(assigned, clusterSize, runID, confidence='low', greater=F)
		    positiveHigh=determine_threshold(assigned, clusterSize, runID, confidence='high', greater=F)
  			return(data.frame(rbind(cbind(cluster=names(positiveLow), positive=positiveLow),
				cbind(cluster=names(positiveHigh), positive=positiveHigh))))
	}
    
  # print("Assigned markers to each cluster.\n")
  # print(positive[positive != xpressd])
  # dVal[positive != xpressd]

}

  
determine_threshold<-function(assigned, clusterSize, runID, confidence='low',
		
		breakStep=0.001, greater = F, markers=c('CD3', 'CD8a', 'CD4') ) {
		breaks=seq(0, 1, breakStep)

		# Determine separation by T cell markers
		combos=gsub(':', '_', get_combos(markers))
		dfD=t(data.frame(assigned))
		colnames(dfD)=c("pval", 'D')
		dfD=as.data.frame(dfD)
		dfD$marker=gsub('.*\\.([^.]+)$', '\\1', rownames(dfD))
		dfD$cluster=gsub('\\.([^.]+)$', '', rownames(dfD))
		dfD$cluster=gsub('\\.', ' ', dfD$cluster)
		dfD$cluster=gsub('^X([0-9]+)$', '\\1', dfD$cluster)
		dfD=subset(dfD, !is.na(D))
		mark=reshape2::dcast(formula = cluster ~ marker, data= dfD,
					   value.var = "D", fill = 0)
		# mark[, -1]=apply(mark[, -1], 2, to_min_max)
		if(confidence=='low') {
			mark=mark[ grep('Excluded', mark$cluster),]
		} else if(confidence == 'high') {
			mark=mark[ grep('Excluded', mark$cluster, invert=T),]
		}
		mark$sizes=clusterSize[match(mark$cluster, names(clusterSize))]
		mark=subset(mark, !is.na(sizes))

		stats=lapply(breaks, function(x) {
			pos=sapply(markers, function(m) {
				vals=rep(NA, nrow(mark))
				if(!greater) {
					ind=mark[, m] >= x  
				} else {
					ind=mark[, m] <= x
				}
				vals[ind]=m
				vals
			})
			if(nrow(mark) == 1) {
                past = paste1(pos[!is.na(pos)])
            } else {
				if(!all(is.na(pos)))
					past=apply(pos, 1, function(s) paste1(s[!is.na(s)]))
			}
			sapply(combos, function(x) {
				if(all(is.na(pos))) return(0)
				sum(mark$sizes[past == x])
			})
		})

		mat=as.data.frame(do.call(rbind, stats))
		if(all(is.na(mat) | mat == 0)) {
			return(sapply(setdiff(colnames(mark), c('cluster', 'sizes')), function(x) ''))
		}
		mat$D=breaks
		signRange=range(dfD$D[dfD$pval <= .1])
		print(signRange)

		pdfOut=f('{runID}.positivity_calling_kstest.{confidence}.pdf')
		pdf(pdfOut, height=4, width = 4.5)

		g <- ggplot(mat[mat$D >= signRange[1] & mat$D <= signRange[2],], aes(D))
		 plot= g +
    		geom_path(aes(y=`CD3_CD8a_CD4` + 1, color="CD3_CD8a_CD4", linetype='CD3_CD8a_CD4'), size=.7) +
	    	geom_path(aes(y=`CD3` + 1, color="CD3"), size=.7) +
	    	geom_path(aes(y=`CD8a_CD4` + 1, linetype='CD8a_CD4', color='CD8a_CD4'), size=.7) +
		    geom_path(aes(y=`CD3_CD8a` + 1 , color='CD3_CD8a'), size=.7) +
	    	geom_path(aes(y=`CD3_CD4` + 1, color='CD3_CD4'), size=.7) +
	    	geom_path(aes(y=`CD4` + 1, color='CD4', linetype='CD4'), size=.7) +
    		geom_path(aes(y=`CD8a` + 1, linetype='CD8a', color='CD8a'), size=.7) +
			scale_color_manual(values= unlist(palette$tcellLines)) +
			scale_linetype_manual(values=unlist(palette$tcellLinetype)) +
	    	cowplot::theme_cowplot() +
    		theme(legend.position = 'top', legend.title = element_blank()) +
		    guides(color=guide_legend(nrow=2), linetype=guide_legend(nrow=2)) +
    		# scale_y_log10() +
	    	ylab('Cell count')
		print(plot)
	
		#tcell_separation=mat$CD3_CD4 + mat$CD3_CD8a - mat$CD3_CD8a_CD4 - mat$`CD8a_CD4`
		tcell_separation=mat$CD3_CD4 + mat$CD3_CD8a + mat$CD4 - mat$CD3_CD8a_CD4 - mat$`CD8a_CD4` - mat$CD3 - mat$CD8a

		tcell_count=sapply(1:nrow(mat), function(x) 
		    mat$CD3_CD4[x] > 0 & mat$CD3_CD8a[x] > 0 & mat$CD4[x] & 
		    mat$D[x] >= signRange[1] & mat$D[x] <= signRange[2])
		if(all(!tcell_count)) tcell_count=sapply(1:nrow(mat), 
			 function(x) mat$CD3_CD4[x] + mat$CD3_CD8a[x]  > 0)

		intervals=diff(tcell_count) == 1
		if(sum(intervals) > 0) {
			# calculate the sizes of each tcell_count == true bin and select the largest bin
			sizes=sapply(1:sum(intervals), function(x){
			if(x == sum(intervals)) {
				end=max(which(tcell_count))
			} else {	
				end = min(which(!tcell_count)[which(!tcell_count) > which(intervals)[x]])
			}
			start=which(intervals)[x]
			sum(tcell_count[start:end])
			})
    		start=which(intervals)[which.max(sizes)]
    		if(which.max(sizes) == sum(intervals)) {
      			end = max(which(tcell_count))
   		 	} else {
				sel=which.max(sizes)
				end = min(which(!tcell_count)[which(!tcell_count) > which(intervals)[sel]])
			}
			tcell_count=rep(F, length(tcell_count))
			tcell_count[start:end] = T
 		 }

		if(!any(tcell_count)) return(NA)

		tcell_pct=sapply(1:nrow(mat), function(x)
			sum(mat$CD3_CD8a_CD4[x] + mat$`CD8a_CD4`[x] + mat$CD8a[x] + mat$CD3[x])/rowSums(mat[x, !colnames(mat) %in% c("V1", 'D')]))
			#sum(mat$CD3_CD8a_CD4[x] + mat$`CD8a_CD4`[x])/sum(mat$CD3_CD4[x] + mat$CD3_CD8a[x] + mat$CD3_CD8a_CD4[x] + mat$`CD8a_CD4`[x]))
		plot(tcell_pct, tcell_separation)
		# for the low confidence group
		# Dsel=min(which.min(tcell_pct), which(tcell_pct<0.05))
		
		#outsiderange = which(!tcell_count)
		#range_up = min(outsiderange[outsiderange > max(which(tcell_count))])
		#range_down = max(outsiderange[outsiderange < min(which(tcell_count))])
		#cat(range_up, range_down, '\n')
		#window=seq(0, 1, 0.01)
		#tcell_fun=sapply(window, function(x) {
		#	if(is.na(range_up)) return(NA)
	#		if(x > mat$D[range_up]) return(NA)
		#	if(is.na(range_down)) return(NA)
		#	if(x < mat$D[range_down]) return(NA)
		##	sum(tcell_pct[which(mat$D %in% x:(x+0.1) & tcell_count)])
		#})
		#d1 = diff(tcell_fun)
		#localInd=max(which(d1 < 0), na.rm = T)
		#localPeak=mat$D > window[localInd]	

		#plot(window, tcell_fun, type = 'l')
		#abline(v=window[localInd])

		# threshold=round(mat$D[Dsel], 2)
		#if(confidence == 'low' & !greater) {
		#	tcell_separation[mat$D < .5]=0  
		#} else if(!greater) {
		#	tcell_separation[mat$D < .3]=0
		#}
	    #	Dsel=which.max(tcell_separation)
		
		# Dsel=intersect(order(diff(tcell_separation)), which(!tcell_count))[1]
		
		# Dsel=intersect(order(diff(tcell_pct)), which(!tcell_count))[1]
		# Dsel=which(tcell_pct == min(tcell_pct[!tcell_count]))
		# Dsel=which(tcell_pct == min(tcell_pct[1:min(which(tcell_count)[-1])]))[1]
			#Dsel=intersect(which(tcell_pct == min(tcell_pct[which(tcell_count & localPeak)], na.rm = T)), which(tcell_count))
		Dsel=intersect(which(tcell_pct == min(tcell_pct[which(tcell_count)], na.rm = T)), which(tcell_count))
		#if(confidence == 'high') {
		#	Dsel = Dsel[1] 
		#} else {
			Dsel = max(Dsel)
		#}

		threshold=round(mat$D[Dsel], 2)
		plot(mat$D, tcell_pct, pch=19, cex=0.1)
		abline(v=threshold)
		abline(v=mat$D[which(!tcell_count)], col=rgb(red = 0.8, blue = 0.8, green = 0.8, alpha = 0.5))
		points(mat$D, tcell_pct, pch=19, cex=0.1, type='l')
		plot(mat$D, tcell_separation, pch=19, cex=0.1)
		abline(v=mat$D[Dsel], lty=2)
		cat('Selected D cutoff ', mat$D[Dsel], ' for confidence ', confidence, '\n')
		print(plot + geom_vline(xintercept=mat$D[Dsel], linetype='dotted', color = 'grey80'))
		dev.off()

		# browseURL(pdfOut)
		posAll=sapply(setdiff(colnames(mark), c('cluster', 'sizes')), function(m) {
			vals=rep(NA, nrow(mark))
			if(!greater) {
				ind=mark[, m] >= threshold
			} else {
				ind=mark[, m] <= threshold
			}
			vals[ind]=m
			vals
		})
		if(nrow(mark) == 1) {
            positive=paste1(posAll[!is.na(posAll)])
        } else if(any(!is.na(posAll))) {
            positive=apply(posAll, 1, function(x) paste1(x[!is.na(x)]))
        } else {
            positive=apply(posAll, 1, function(x) '')
        }
		cat('Output saved in ', f('{runID}.positivity_calling_kstest.{confidence}.pdf'), '\n')
		names(positive)=mark$cluster
		return(positive)
}

get_probabilities <- function(cellIDs, clusters, probFile) {
  
		if(!file.exists(probFile))
		return(rep(NA, nrow(inData)))
		probs=data.table::fread(probFile, sep = "\t")
		probMatch=match(cellIDs, probs$imagename)
		probs=probs[probMatch, ]
	#	print(colnames(probs))
		export=tapply(1:length(cellIDs), clusters, function(subset) {
			cluster=clusters[subset[1]]
			if(! cluster %in% colnames(probs)) {
				cat('No probability values for cluster', cluster, '\n')
				values=as.matrix(rep(NA, length(subset)), ncol=1)  # Excluded cells
			} else {
					values=as.matrix(probs[subset, ..cluster, with=F])
			}
			rownames(values)=cellIDs[subset]
			return(values)
		})
		export=do.call(rbind, export)
		export[match(cellIDs, rownames(export)), 1]
}

get_markers_expression <- function(dfExp, clusters, clusterNames, magnitude=NULL, fun="max") {

		if(!is.null(magnitude)) dfExp=to_magnitude(dfExp, magnitude)
		clusterNames$cluster=gsub('\\.', ' ', clusterNames$cluster)
		if('CD4' %in% clusters)
			clusterNames$cluster = gsub('CD([48]) T cells', 'CD\\1', clusterNames$cluster)
	#	print(clusterNames)
	#	print(table(clusters))
		names=clusterNames$positive[match(clusters, clusterNames$cluster)]
		names=gsub("up:(.*) down:.*", "\\1", names)
		names=gsub("pos:(.*) neg:.*", "\\1", names)
	#	print(table(clusters[is.na(names)]))
		names[is.na(names)] = ''
		getExpSummary=match.fun(fun)
		summary=tapply(1:nrow(dfExp), names, function(subset) {
				cluster=names[subset[1]]
				cat('Expression values for', cluster, '\n')
				if(cluster == "") {
					values=rep(NA, length(subset))
				} else {
					cols=strsplit(cluster, split = "_")[[1]]
					cols=cols[cols %in% colnames(dfExp)]
					values=apply(dfExp[subset, ..cols], 1, getExpSummary)
				}
				names(values)=rownames(dfExp)[subset]
				return(cbind(values))
		})
		summary=do.call(rbind, summary)
		summary[match(rownames(dfExp), rownames(summary)), 1]
}
