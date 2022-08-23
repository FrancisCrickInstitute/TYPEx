# Generic function
  
assign_cluster_positivity <- function(dfExp, clusters, panel, run, runID,
			fdr=0.05, log2=F, magnitude=10**5,
			fromList=NA)  {
				
	# # panel="P1"; run="B"; fdr=0.05; log2=F; magnitude=10**6; runID="~/labwd/analyses/typing/test"
	effect_field='AUC'
	
	scranPkg=packageVersion("scran")
	if(scranPkg == "1.12.1") effect_field='overlap'
	if(length(unique(clusters)) == 1)
		return(
			data.frame(cluster=rep(NA, length(clusters)), 
					   name=rep(NA, length(clusters)))
		)
	if(! is.null(magnitude)) {
		
		cat('Magnitude transformation ', magnitude, '\n')
		
		dfExp=to_magnitude(dfExp, magnitude)
	}
	
	if(log2) 
	  dfExp=log2(dfExp + .1)
	clusterSize=table(clusters)
	clusterLabels=names(clusterSize)


	if(! is.na(fromList)) {
		cat('NOTE: Assigning based on the list of markers for the identified cell type\n')
		positive=sapply(clusterLabels, function(cluster) {
		  expressedMarkers=get_celltype_markers(marker_gene_list[[fromList]], cluster)
		  paste0("up:", paste1(expressedMarkers), " ", "down:")
		})
		return(data.frame(
					cluster=names(positive),
					positive=positive))
	}
	
	# Not required for the latest version of the scran package
	if(! file.exists(f('{runID}.pairwise.RData'))) {

		print('Runing pairwise')
		# Get probabilities of overlap
		wilcoxTestUp=scran::pairwiseWilcox(
			t(dfExp), 
			as.factor(clusters), 
			direction="up")
		# save(wilcoxTestUp, file = f("${runID}.wtu.RData"))
		# print(packageVersion("scran"))
		markersUp=scran::combineMarkers(
			wilcoxTestUp$statistics, 
			wilcoxTestUp$pairs, 
			effect.field = effect_field) # previously overlap
		rowNames=rownames(markersUp[[1]])
		# Data frame with probabilities for each pairwise comparison
		upDf=do.call(cbind, lapply(markersUp, function(x) x[rowNames, ]))
		upDf=upDf[, grep(effect_field, colnames(upDf))]
		# Summarise per marker
		upStats=apply(upDf, 1, mean) + apply(upDf, 1, sd)
		upDf=data.frame(upDf, check.names = F)
		upDf = upDf %>% t %>%  as.data.frame
		save(upDf, upStats, file = f('{runID}.pairwise.RData'))
	
	} else  {
		load(f('{runID}.pairwise.RData'))
	}
  
	
	print('Assigning')
    upDfConf=upDf[grep(f('.{effect_field}.Excluded.*'), rownames(upDf), invert = T), ]
	upDfConf=upDfConf[grep('summary', rownames(upDfConf), invert = T),]
    # New way
    assigned=plot_overlaps(upDfConf, clusterLabels, runID, effect_field=effect_field)
	if(pars$subset %in% c(majorDir, sampledDir) | ! pars$stratify ) {
	    print('Determining threshold without stratification')
		pars$threshold %>% print
		positive=determine_threshold(
			assigned, 
			clusterSize, 
			runID, 
			confidence='all',
			greater=F,
			markers=pars$threshold[["markers"]])
		return(
			data.frame(
				rbind(
					cbind(
						cluster=names(positive), positive=positive))
		))
	} else {
		print('Determining threshold with stratification')
		positiveHigh=determine_threshold(
				assigned,
				clusterSize,
				runID,
				confidence='high',
				greater=F, 
				markers=pars$threshold[["markers"]])
				
		positiveLow=determine_threshold(
				assigned, 
				clusterSize, 
				runID, 
				confidence='low',
				greater=F,
				markers=pars$threshold[["markers"]])

		return(
			data.frame(rbind(cbind(cluster=names(positiveLow), 
									positive=positiveLow),
							cbind(cluster=names(positiveHigh),
									positive=positiveHigh))
				)
			)
	}
}

determine_threshold <- function(assigned, clusterSize, runID, confidence='low',
		breakStep=0.001, greater = F, markers=c('CD3', 'CD8a', 'CD4')) {
		breaks=seq(0, 1, breakStep)
		cat("Determining threshold with markers: ", markers, '\n')
		# Determine separation by T cell markers
		if(! is.null(markers)) {
			combos=gsub(':', '_', get_combos(markers))
		} else {
			stop("ERROR: Cell type-specific markers not provided for positivity calling")
		}
		


		dfD=data.frame(assigned, check.names = F) %>% t %>% 
			as.data.frame
		colnames(dfD)=c("pval", 'D')
		dfD$marker=rownames(dfD) %>%
			gsub('.*\\.([^.]+)$', '\\1', .)
		dfD$cluster=rownames(dfD) %>% 
			gsub('\\.([^.]+)$', '', .) %>% 
			gsub('\\.', ' ', .) %>% 
			gsub('^X([0-9]+)$', '\\1', .)
		dfD=subset(dfD, !is.na(D))
		mark=reshape2::dcast(formula = cluster ~ marker, data= dfD,
					   		value.var = "D", fill = 0)

		# mark[, -1]=apply(mark[, -1], 2, to_min_max)
		if(confidence == 'low') {
			mark=mark[ grep('Excluded', mark$cluster),]
		} else if(confidence == 'high') {
			mark=mark[ grep('Excluded', mark$cluster, invert=T),]
		}
		if(nrow(mark) == 0) {
			cat('No clusters with ', confidence, ' confidence\n')
			return(NULL)
		}
		mark$sizes=clusterSize[match(mark$cluster, names(clusterSize))]
		mark=subset(mark, ! is.na(sizes))

		stats=lapply(breaks, function(x) {
			pos=sapply(markers, function(m) {
				vals=rep(NA, nrow(mark))
				if(! greater) {
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
				if(! all(is.na(pos)))
					past=apply(pos, 1, function(s) paste1(s[! is.na(s)]))
			}
			sapply(combos, function(x) {
				if(all(is.na(pos))) return(0)
				sum(mark$sizes[past == x])
			})
		})

		mat=do.call(rbind, stats) %>% as.data.frame
		if(all(is.na(mat) | mat == 0)) {
			sapply(setdiff(colnames(mark), c('cluster', 'sizes')), '') %>%
				return
		}
		mat$D=breaks
		signRange=range(dfD$D[dfD$pval <= .1])
		print(signRange)

		pdfOut=f('{runID}.positivity_kstest.{confidence}.pdf')
		pdf(pdfOut, height=4, width = 4.5)

		frequences=sapply(c("high_frequency", 'rare', 'low_frequency'), 
			function(x) pars$threshold[[x]]) %>% unlist
		if(! all(frequences %in% colnames(mat))) {
			cat(frequences[! frequences %in% colnames(mat)], ' is missing in ', 
				colnames(mat), '\n')
			stop("ERROR: Verify that the marker combinations in stratification.json as", 
				pars$stratify_label, "are valid")
		}
		tcell_separation=apply(mat[, pars$threshold[['high_frequency']]], 1, sum) - 
						 apply(mat[, pars$threshold[['rare']]], 1, sum)
		print(colnames(mat))
		
		# only using high frequency as low may not be found in a small dataset
		tcell_count=sapply(1:nrow(mat), function(x) {
			all(
				sapply(pars$threshold[['high_frequency']], 
					function(combination) mat[x, combination] > 0)
			) &  mat$D[x] >= signRange[1] & 
			mat$D[x] <= signRange[2]
		})
		
		# Previously excluded CD4 requirement here
		if(all(! tcell_count)) 
			tcell_count=sapply(1:nrow(mat), 
				function(x) sum(mat[x, pars$threshold[['high_frequency']] ])  > 0)
		
		intervals=diff(tcell_count) == 1
		if(sum(intervals) > 0) {

			# calculate the sizes of each tcell_count == true bin and select the largest bin
			sizes=sapply(1:sum(intervals), function(x)	{
				if(x == sum(intervals)) {
					end=max(which(tcell_count))
				} else {
					end = min(which(!tcell_count)[which(! tcell_count) > which(intervals)[x]])
				}
				start=which(intervals)[x]
				sum(tcell_count[start:end])
			})
    		start=which(intervals)[which.max(sizes)]
    		if(which.max(sizes) == sum(intervals))	{
      			end = max(which(tcell_count))
   		 	} else {
				sel=which.max(sizes)
				end = min(which(! tcell_count)[which(! tcell_count) > which(intervals)[sel]])
			}
			tcell_count=rep(F, length(tcell_count))
			tcell_count[start:end] = T
 		 }
		 
 		tcell_pct=sapply(1:nrow(mat), function(x)
			sum(mat[x, pars$threshold[['high_frequency']]]) / 
						rowSums(mat[ x, ! colnames(mat) %in% c("V1", 'D')])
		)
 		plot(tcell_pct, tcell_separation)

		if(! any(tcell_count))	{
			cat('WARNING: No T cell count -> taking a threshold of 0.5,', 
				' which matches wilcox.test differential expression at p=0.05\n',
				append = T, file = f("{runID}.log"))
			threshold=0.5
		} else {
			cat('INFO: Minimum tcell_count', min(tcell_pct[which(tcell_count)], na.rm = T), '\n',
				append = T, file = f("{runID}.log"))
			cat('INFO: Range', range(tcell_pct[which(tcell_count)], na.rm = T), '\n',
				append = T, file = f("{runID}.log"))
			cat('INFO: Minimum tcell_count', min(tcell_pct[which(tcell_count)], na.rm = T), '\n')
			cat('INFO: Range', range(tcell_pct[which(tcell_count)], na.rm = T), '\n')
			Dsel=intersect(which(tcell_pct == min(tcell_pct[which(tcell_count)], na.rm = T)), which(tcell_count))
			print(tcell_pct[Dsel])
			cat('INFO: Dsel options', Dsel, '\n')
			Dsel = max(Dsel)
			threshold=round(mat$D[Dsel], 2)
			plot(mat$D, tcell_pct, pch=19, cex=0.1)
			abline(v=threshold, color = 'red', lty =2)
			abline(v=mat$D[which(!tcell_count)], 
				col=rgb(red = 0.8, blue = 0.8, green = 0.8, alpha = 0.5))
			points(mat$D, tcell_pct, pch=19, cex=0.1, type='l')
			plot(mat$D, tcell_separation, pch=19, cex=0.1)
			abline(v=mat$D[Dsel], lty=2)
			cat('Selected D cutoff ', mat$D[Dsel], 
				' for confidence ', confidence, '\n')
			cat('INFO: Selected D cutoff ', mat$D[Dsel], 
				' for confidence ', confidence, '\n', append = T,
				file = f("{runID}.log"))
				
			g <- ggplot(mat, aes(D))
				         plot= g +
				            geom_path(aes(y=`CD3_CD4_CD8` + 1, color="CD3_CD8a_CD4", linetype='CD3_CD8a_CD4'), size=.7) +
				            geom_path(aes(y=`CD3` + 1, color="CD3"), size=.7) +
				            geom_path(aes(y=`CD4_CD8` + 1, linetype='CD8a_CD4', color='CD8a_CD4'), size=.7) +
				            geom_path(aes(y=`CD3_CD8` + 1 , color='CD3_CD8a'), size=.7) +
				            geom_path(aes(y=`CD3_CD4` + 1, color='CD3_CD4'), size=.7) +
				            geom_path(aes(y=`CD4` + 1, color='CD4', linetype='CD4'), size=.7) +
				            geom_path(aes(y=`CD8` + 1, linetype='CD8a', color='CD8a_CD4a'), size=.7) +
				            scale_color_manual(values= unlist(palette$tcellLines)) +
				            scale_linetype_manual(values=unlist(palette$tcellLinetype)) +
				            cowplot::theme_cowplot() +
				            theme(legend.position = 'top', legend.title = element_blank()) +
				            guides(color=guide_legend(nrow=3), linetype=guide_legend(nrow=3)) +
				            # scale_y_log10() +
				            ylab('Cell count')
			print(plot)
			print(plot + geom_vline(xintercept=threshold, linetype='dotted', color = 'grey80'))

			dev.off()
		}
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
		
		if(nrow(mark) == 1)	{
            positive=paste1(posAll[!is.na(posAll)])
        } else if(any(!is.na(posAll))) {
            positive=apply(posAll, 1, function(x) paste1(x[!is.na(x)]))
        } else {
            positive=apply(posAll, 1, function(x) '')
        }
		cat('Output saved in ', f('{runID}.positivity_kstest.{confidence}.pdf'), '\n')
		names(positive)=mark$cluster
		return(positive)
}

assign_celltype <- function(names, markers, majorTypes=NULL, region=NULL, major=F) {

	markerNames=unlist(markers) %>% unique
	cellTypeNames=get_celltypes(markers)
	if(all(majorTypes %in% c('', 'Excluded')))
		majorTypes=NULL
	
	if(!is.null(majorTypes)) 
		names=paste(majorTypes, names, sep =  "|")
	if(!is.null(region)) 
		names=paste(region, names, sep =  "=")
	
	clusters=sapply(unique(names), toString)
	markerFrequency=unlist(markers) %>% table
	
	
	markersList=get_node_list(markers)
	
	
	celltypes=sapply(clusters, function(name) {
		if(is.na(name))
			name=""
		region=strsplit(name, split="=")[[1]][1]
		if(length(grep("\\|", region)))
			region="none"
		if(is.na(region) | region=='NA' | region=='')
			region="none"
	    if(region == name)
			region='none'

	    majorType=gsub("^[^=]+=", "", name) %>% 
			gsub("^([^|]+)\\|.*", "\\1", .)
	    if(is.null(majorTypes))
			majorType='none'
	    if(length(grep("Excluded", majorType)))
			majorType="none"
	    if(length(grep("\\|", majorType)))
			majorType="none"

	    name=gsub("^[^|]+\\|", "", name) %>%
			gsub("^[^=]+=", "", .) %>%
			gsub("^[^:]+:([^ ]+).*", "\\1", .)

	    if(name %in% c("", "|", "NA"))
	      if(majorType == "none" & region == "none") {
	        return("Unassigned")
	      } else if(majorType != "none" & 
		  	majorType != 'Epithelial cells') {
	        return(majorType)
	      } else if(region != 'none' & ! is.na(region) & region > 0) {
	        # added if region == 1
	        return("Epithelial cells - Tissue Segmentation")
	      } else if(majorType != 'none') {
	        return(majorType)
	    }

	    if(name == "Excluded")
			return(name)
		
		# cat(name, region, majorType, '\n')
	
		# Split, format and select the celltype-specifc marker
	    cellMarkers=strsplit(name, split="_")[[1]]
	    cellMarkers=markers_format(cellMarkers)
	    cellMarkers=intersect(cellMarkers, markerNames)
	
		if(! length(cellMarkers)) {

	      if(grepl('Epithelial cells', majorType) & 
			 region != 'none' & region != name & 
		  	 region > 0) 
			  return("Epithelial cells - Tissue Segmentation")
	      if(majorType != "none") 
			  return(majorType)
	      if(major & majorType != 'none') 
			  return('Unassigned')
	      if(region != 'none' & region != name & region > 0) 
			  return("Epithelial cells - Tissue Segmentation")
	      if(majorType == 'Mesenchymal cells') 
			  return(majorType)
	      return("Unassigned")
	    }
		
	    types=sapply(cellTypeNames, function(celltype) {
			cellspecific_markers=markersList[[celltype]]
			all(cellMarkers %in% cellspecific_markers) &
				length(cellMarkers) == length(cellspecific_markers)
		})
	    if(sum(types) == 1)  return(names(types)[types])

	    intersect=sapply(cellMarkers, function(.mark1) {
	      # if(.mark1 == "alphaSMA") return(1)
	      sel1=sapply(cellTypeNames, function(celltype) 
			  	.mark1 %in% markersList[[celltype]]
		)
	      sel1=cellTypeNames[sel1]
	      sum(sapply(cellMarkers, function(.mark2) {
	        if(.mark1==.mark2) return(0)
				sel2=sapply(cellTypeNames, function(celltype) 
					.mark2 %in% markersList[[celltype]]
				)
	        sel2=cellTypeNames[sel2]
	        length(intersect(sel1, sel2))
	      }))
	    })
		
		if(all(intersect == length(cellMarkers) - 1) & length(intersect) > 1) {
			overlap=sapply(names(markers), function(celltype) {
				all(cellMarkers %in% markersList[[celltype]])
			})
			if(any(overlap)) return(names(overlap)[overlap])
		}
		
		specificity=sapply(cellMarkers, function(.mark) {
			# Vimentin is the only marker that distinguishes uniquely 
			# a cell type but is found in others, EMT
			if(.mark == "Vimentin" & ! major) 
				return(F)
			specific=unlist(markers) %in% .mark

			return(sum(specific) == 1)
		})
		
		# Consider if the cell belondgs to the tumour region
		if(sum(specificity) > 1 & majorType == 'none' & 
		  region != 'none' & region != name & region > 0) 
			return("Epithelial cells - Tissue Segmentation")
		if(sum(specificity) > 1 & majorType == 'none')
			return("Ambiguous")
		
		#overlap=sapply(cellTypeNames, function(celltype) {
		#	sum(cellMarkers %in% get_celltype_markers(markers, celltype))
		#})
		
		assigned=sapply(cellTypeNames, function(celltype) {
			cellspecific_markers=markersList[[celltype]]
			if(length(cellspecific_markers) == 0) return(0)

			sum(sapply(cellspecific_markers, function(marker) {
				if(! marker %in% names(markerFrequency)) return(0)
					1 / markerFrequency[[marker]] * sum(marker %in% cellMarkers) *
				sum(cellMarkers %in% cellspecific_markers) / length(cellspecific_markers)
			}))
		})
		# sort(assigned)
		maxas=max(assigned, na.rm  = T)
		
		if(sum(assigned==maxas & maxas > 0, na.rm = T) > 1) {
			if(! majorType %in% c("none")) 
				return(majorType) 
			if(region != 'none' & region != name & region > 0) 
				return("Epithelial cells - Tissue Segmentation")
			return("Ambiguous")
		}
		if(all(assigned == 0) & majorType != 'none') {
			if(region != 'none' & region != name & region > 0) 
				return("Epithelial cells - Tissue Segmentation")
			if(majorType != 'Mesenchymal cells') 
				return("Unassigned")
			return(majorType)
		}
		
		if(length(intersect) > 1 & all(intersect == 0) & majorType != 'none') {
			if(! majorType %in% c("none") & !major)
				return(majorType)
			if(majorType != 'none' & region != 'none' & region != name & region > 0 & ! major)
				return("Epithelial cells - Tissue Segmentation")
			if(majorType != "none") 
				return(majorType)
			return("Ambiguous")
		}
		
		if(major & sum(assigned == maxas, na.rm = T) > 1 & majorType != "none") 
			return(majorType)
		if(major & sum(assigned == maxas, na.rm = T) > 1 & 
				majorType == "none" & 
				region != 'none' & region != name &
				region > 0) 
			return("Epithelial cells - Tissue Segmentation")
			
		if(major & sum(assigned==maxas, na.rm = T) > 1 & majorType == "none") 
			return('Unassigned')
		
		assignedCelltype=names(assigned)[which(assigned==maxas)]
		# See if on the subtree of the majorType
		if(! majorType %in% c('none', 'Unassigned'))
			if(! isChild(tree=markers, parent=majorType, child=assignedCelltype))
  			 	return('Ambiguous - For Review')
		  
		return(assignedCelltype)
	})
	celltypes[match(names, clusters)]
}
