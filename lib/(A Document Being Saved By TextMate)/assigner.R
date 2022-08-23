# Generic function
  
assign_cluster_positivity <- function(dfExp, clusters, panel, run, runID,
			fdr=0.05, log2=F, magnitude=10**5,
			fromList=NA, ilastik=T)  {
				
	# # panel="P1"; run="B"; fdr=0.05; log2=F; magnitude=10**6; runID="~/labwd/analyses/typing/test"
	effect_field='AUC'
	scranPkg=packageVersion("scran")
	if(scranPkg == "1.12.1") effect_field='overlap'
	if(length(unique(clusters)) == 1)
	return(
		data.frame(cluster=rep(NA, length(clusters)),
				   name=rep(NA, length(clusters))))
	if(!is.null(magnitude)) 
	  dfExp=to_magnitude(dfExp, magnitude)
	if(log2) 
	  dfExp=log2(dfExp + .1)
	clusterSize=table(clusters)
	clusterLabels=names(clusterSize)

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
  
	if(! is.na(fromList)) {
	positive=sapply(clusterLabels, function(cluster) {
	  upRows=marker_gene_list[[fromList]][[cluster]]
	  if(ilastik) upRows=intersect(upRows, marker_gene_list$major_markers)
	  paste0("up:", paste1(upRows), " ", "down:")
	})
	return(data.frame(
				cluster=names(positive),
				positive=positive))
	}
	print('Assigning')
    upDfConf=upDf[grep(f('.{effect_field}.Excluded.*'), rownames(upDf), invert = T), ]
	upDfConf=upDfConf[grep('summary', rownames(upDfConf), invert = T),]
    # New way
    assigned=plot_overlaps(upDfConf, clusterLabels, runID, effect_field=effect_field)
	if(pars$subset %in% c(majorDir, sampledDir) | ! pars$stratify ) {

		positive=determine_threshold(
			assigned, 
			clusterSize, 
			runID, 
			confidence='all',
			greater=F)
		return(
			data.frame(
				rbind(
					cbind(
						cluster=names(positive), positive=positive))
		))
	} else {
		print('Determining threshold with stratification')
		positiveLow=determine_threshold(
				assigned, 
				clusterSize, 
				runID, 
				confidence='low',
				greater=F)
				print('positive high')
		positiveHigh=determine_threshold(
				assigned,
				clusterSize,
				runID,
				confidence='high',
				greater=F)
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
		breakStep=0.001, greater = F, markers=c('CD3', 'CD8a', 'CD4') ) {
	
		breaks=seq(0, 1, breakStep)

		# Determine separation by T cell markers
		combos=gsub(':', '_', get_combos(markers))
		dfD=data.frame(assigned, check.names = F) %>% t %>% 
			as.data.frame
		colnames(dfD)=c("pval", 'D')
		dfD$marker=gsub('.*\\.([^.]+)$', '\\1', rownames(dfD))
		dfD$cluster=gsub('\\.([^.]+)$', '', rownames(dfD)) %>% 
			gsub('\\.', ' ', .) %>% gsub('^X([0-9]+)$', '\\1', .)
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
		mark=subset(mark, ! is.na(sizes))

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
				if(! all(is.na(pos)))
					past=apply(pos, 1, function(s) paste1(s[! is.na(s)]))
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

		pdfOut=f('{runID}.positivity_kstest.{confidence}.pdf')
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
		tcell_separation=
			mat$CD3_CD4 + 
			mat$CD3_CD8a + mat$CD4 - 
			mat$CD3_CD8a_CD4 - mat$`CD8a_CD4` -
			mat$CD3 - mat$CD8a

		tcell_count=sapply(1:nrow(mat), function(x) 
		    mat$CD3_CD4[x] > 0 & 
			mat$CD3_CD8a[x] > 0 & 
			mat$CD4[x] & 
		    mat$D[x] >= signRange[1] & 
			mat$D[x] <= signRange[2]
		)
		
		if(all(!tcell_count)) tcell_count=sapply(1:nrow(mat), 
			 function(x) mat$CD3_CD4[x] + mat$CD3_CD8a[x]  > 0)

		intervals=diff(tcell_count) == 1
		if(sum(intervals) > 0) {

			# calculate the sizes of each tcell_count == true bin and select the largest bin
			sizes=sapply(1:sum(intervals), function(x)	{
				if(x == sum(intervals)) {
					end=max(which(tcell_count))
				} else {
					end = min(which(!tcell_count)[which(!tcell_count) > which(intervals)[x]])
				}
				start=which(intervals)[x]
				sum(tcell_count[start:end])
			})
    		start=which(intervals)[which.max(sizes)]
    		if(which.max(sizes) == sum(intervals))	{
      			end = max(which(tcell_count))
   		 	} else {
				sel=which.max(sizes)
				end = min(which(!tcell_count)[which(! tcell_count) > which(intervals)[sel]])
			}
			tcell_count=rep(F, length(tcell_count))
			tcell_count[start:end] = T
 		 }
		 
 		tcell_pct=sapply(1:nrow(mat), function(x)
 			sum(mat$CD3_CD8a_CD4[x] + mat$`CD8a_CD4`[x] + mat$CD8a[x] +
 				mat$CD3[x])/rowSums(mat[ x, ! colnames(mat) %in% c("V1", 'D')])
		)
 		plot(tcell_pct, tcell_separation)

		if(! any(tcell_count))	{
			cat('No T cell count ->  taking a threshold of 0.5,', 
				' which matches wilcox.test differential expression at p=0.05\n')
			threshold=0.5
		} else {
			
			Dsel=intersect(which(tcell_pct == min(tcell_pct[which(tcell_count)], na.rm = T)), which(tcell_count))
			Dsel = max(Dsel)

			threshold=round(mat$D[Dsel], 2)
			plot(mat$D, tcell_pct, pch=19, cex=0.1)
			abline(v=threshold)
			abline(v=mat$D[which(!tcell_count)], 
				col=rgb(red = 0.8, blue = 0.8, green = 0.8, alpha = 0.5))
			points(mat$D, tcell_pct, pch=19, cex=0.1, type='l')
			plot(mat$D, tcell_separation, pch=19, cex=0.1)
			abline(v=mat$D[Dsel], lty=2)
			cat('Selected D cutoff ', mat$D[Dsel], 
				' for confidence ', confidence, '\n')
			print(plot + geom_vline(xintercept=mat$D[Dsel],
				linetype='dotted', color = 'grey80'))
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

	    if(name=="Excluded")
			return(name)
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
			cellspecific_markers=get_celltype_markers(markers, celltype)
			all(cellMarkers %in% cellspecific_markers) &
				length(cellMarkers) == length(cellspecific_markers)
		})
	    if(sum(types) == 1)  return(names(types)[types])

	    intersect=sapply(cellMarkers, function(.mark1) {
	      # if(.mark1 == "alphaSMA") return(1)
	      sel1=sapply(cellTypeNames, function(celltype) 
			  	.mark1 %in% get_celltype_markers(markers, celltype)
		  )
	      sel1=cellTypeNames[sel1]
	      sum(sapply(cellMarkers, function(.mark2) {
	        if(.mark1==.mark2) return(0)
				sel2=sapply(cellTypeNames, function(celltype) 
					.mark2 %in% get_celltype_markers(markers, celltype)
				)
	        sel2=cellTypeNames[sel2]
	        length(intersect(sel1, sel2))
	      }))
	    })
		
		
		if(all(intersect == length(cellMarkers) - 1) & length(intersect) > 1) {
			overlap=sapply(names(markers), function(celltype) {
				all(cellMarkers %in% get_celltype_markers(markers, celltype))
			})
			if(any(overlap)) return(names(overlap)[overlap])
		}
		
		specificity=sapply(cellMarkers, function(.mark) {
			# Vimentin is the only marker that distinguishes uniquely 
			# a cell type but is found in others, EMT
			if(.mark == "Vimentin" & !major) 
				return(F)
			specific=unlist(markers) %in% .mark
			
			# #specific=sapply(names(markers), function(celltype) {
				#  all(markers[[celltype]] == .mark)
				#})
			return(sum(specific) == 1)
		})
		
		# Consider if the cell belondgs to the tumour region
		if(sum(specificity) > 1 & majorType == 'none' & 
		  region != 'none' & region != name & region > 0) 
			return("Epithelial cells - Tissue Segmentation")
		if(sum(specificity) > 1 & majorType == 'none')
			return("Ambiguous")
		
		overlap=sapply(cellTypeNames, function(celltype) {
			sum(cellMarkers %in% get_celltype_markers(markers, celltype))
		})
		
		assigned=sapply(cellTypeNames, function(celltype) {
			cellspecific_markers=get_celltype_markers(markers, celltype)
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


assign_celltype_original <- function(names, markers, majorTypes=NULL, region=NULL, major=F) {

	markerNames=unlist(markers) %>% unique
	cellTypeNames=get_celltypes(markers)
	if(all(majorTypes %in% c('', 'Excluded')))
		majorTypes=NULL
	
	if(! is.null(majorTypes))
		names=paste(majorTypes, names, sep =  "|")
	if(! is.null(region))
		names=paste(region, names, sep =  "=")
	
	clusters=sapply(unique(names), toString)
	markerFrequency=unlist(markers) %>% table
	
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
	      } else if(majorType != "none" & majorType != 'Epithelial cells') {
	        return(majorType)
	      } else if(region != 'none' & ! is.na(region) & region > 0) {
	        # added if region == 1
	        return("Epithelial cells - Tissue Segmentation")
	      } else if(majorType != 'none') {
	        return(majorType)
	    }

	    if(name=="Excluded") return(name)
	    cellMarkers=strsplit(name, split="_")[[1]] %>%
			markers_format %>% 
			intersect(., markerNames)
	    # if(majorType != "none") cellMarkers=setdiff(cellMarkers, nonspecificAb[[majorType]]) 
	    if(!length(cellMarkers)) {

	      if(grepl('Epithelial cells', majorType) & region != 'none' & region != name & region > 0) 
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
			cellspecific_markers=get_celltype_markers(markers, celltype)
			all(cellMarkers %in% cellspecific_markers) &
				length(cellMarkers) == length(cellspecific_markers)
		})
	    if(sum(types) == 1)  return(names(types)[types])

	    intersect=sapply(cellMarkers, function(.mark1) {
	      # if(.mark1 == "alphaSMA") return(1)
	      sel1=sapply(cellTypeNames, function(celltype) 
			  	.mark1 %in% get_celltype_markers(markers, celltype)
		  )
	      sel1=cellTypeNames[sel1]
	      sum(sapply(cellMarkers, function(.mark2) {
	        if(.mark1==.mark2) return(0)
				sel2=sapply(cellTypeNames, function(celltype) 
					.mark2 %in% get_celltype_markers(markers, celltype)
				)
	        sel2=cellTypeNames[sel2]
	        length(intersect(sel1, sel2))
	      }))
	    })
		
		
		if(all(intersect == length(cellMarkers) - 1) & length(intersect) > 1) {
			overlap=sapply(names(markers), function(celltype) {
				all(cellMarkers %in% get_celltype_markers(markers, celltype))
			})
			if(any(overlap)) return(names(overlap)[overlap])
		}
		
		specificity=sapply(cellMarkers, function(.mark) {
			# Vimentin is the only marker that distinguishes uniquely 
			# a cell type but is found in others, EMT
			if(.mark == "Vimentin" & !major) 
				return(F)
			specific=unlist(markers) %in% .mark
			return(sum(specific) == 1)
		})
		
		# Consider if the cell belondgs to the tumour region
		if(sum(specificity) > 1 & majorType == 'none' & region != 'none' & region != name & region > 0) 
			return("Epithelial cells - Tissue Segmentation")
		if(sum(specificity) > 1 & majorType == 'none')
			return("Ambiguous")
		
		overlap=sapply(cellTypeNames, function(celltype) {
			sum(cellMarkers %in% get_celltype_markers(markers, celltype))
		})
		
		assigned=sapply(cellTypeNames, function(celltype) {
			cellspecific_markers=get_celltype_markers(markers, celltype)
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
		
		if(major & sum(assigned==maxas, na.rm = T) > 1 & majorType != "none") 
			return(majorType)
		if(major & sum(assigned==maxas, na.rm = T) > 1 & 
				majorType == "none" & 
				region != 'none' & region != name &
				region > 0) 
			return("Epithelial cells - Tissue Segmentation")
			
		if(major & sum(assigned==maxas, na.rm = T) > 1 & majorType == "none") 
			return('Unassigned')
		names(assigned)[which(assigned==maxas)]
	})
	celltypes[match(names, clusters)]
}

