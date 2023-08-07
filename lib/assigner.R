# Generic function
  

determine_threshold <- function(assigned, clusterSize, runID, confidence='low',
		
		breakStep=0.0001, greater = F, markers=c('CD3', 'CD8a', 'CD4')) {
		breaks=seq(0, 1, breakStep)
		smallImageSetCutoff = 100000
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
		if(! all(markers %in% colnames(mark))) {
			print(markers)
		   print(colnames(mark))
			stop('The marker names defined for thresholding in typing_params.json are not in the cell objects file')
		}

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
		mat$D = breaks
		signRange = range(dfD$D[dfD$pval <= .1])
		cat('Sign range:', signRange, '\n')

		plotDir = f('{runID}_plots')
		if(! dir.exists(plotDir))
			dir.create(plotDir)
		
		pdfOut = f('{plotDir}/threshold.{confidence}_confidence.pdf')
		pdf(pdfOut, height=4, width = 4.5)

		frequences=sapply(c("high_frequency", 'rare', 'low_frequency', 'variable'), 
			function(x) pars$threshold[[x]]) %>% unlist
		if(! all(frequences %in% colnames(mat))) {
			cat(frequences[! frequences %in% colnames(mat)], ' is missing in ', 
				colnames(mat), '\n')
			stop("ERROR: Verify that the marker combinations in ",
				"typing_params.json are valid")
		}
		
		# If there are few images, variable is likely to be low
		if(sum(clusterSize) > smallImageSetCutoff) {
			present=sapply(c("high_frequency", 'variable'),
	            function(x) pars$threshold[[x]]) %>% unlist
		} else {
			# If there are few images, variable is likely to be low
			present=sapply(c("high_frequency"),
	            function(x) pars$threshold[[x]]) %>% unlist
		}
		rare=sapply(c("low_frequency", 'rare'),
            function(x) pars$threshold[[x]]) %>% unlist
		lower_freq=sapply(c("low_frequency", 'variable'),
			function(x) pars$threshold[[x]]) %>% unlist
		# Convert to DF in case one column/combination specified
		tcell_separation=apply(as.data.frame(mat[, present]), 1, sum) -
					apply(as.data.frame(mat[, rare]), 1, sum)
					
		tcell_count=sapply(1:nrow(mat), function(x) {
			presentCombos = present
			all(
				sapply(
					presentCombos, function(combination) mat[x, combination] > 0
				)
			) & 
			mat$D[x] >= signRange[1] & 
			mat$D[x] <= signRange[2]
		})
		print(table(tcell_count))
		# only using high frequency low may not be found in a small dataset
		if(all(! tcell_count)) 
			tcell_count=sapply(1:nrow(mat), 
				function(x) all(mat[x, pars$threshold[['high_frequency']] ]) > 0 & 
							mat$D[x] >= signRange[1] & 
							mat$D[x] <= signRange[2])
		print(table(tcell_count))
		
		if(all(! tcell_count)) 
			tcell_count=sapply(1:nrow(mat), 
				function(x) sum(mat[x, pars$threshold[['high_frequency']] ]) > 0)
		
		signInd=sapply(1:nrow(mat), function(x) {
			mat$D[x] >= signRange[1] & mat$D[x] <= signRange[2]
		})
		
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
    		start=max(which(intervals)[which.max(sizes)], min(signInd))
    		if(which.max(sizes) == sum(intervals))	{
      			end = max(which(tcell_count))
   		 	} else {
				sel=which.max(sizes)
				end = min(which(! tcell_count)[which(! tcell_count) > which(intervals)[sel]])
			}
			tcell_count=rep(F, length(tcell_count))
			tcell_count[start:end] = T
 		 }
		
 		minimise_combos = sapply(c('rare'), #, 'low_frequency'
			function(x) pars$threshold[[x]]) %>% unlist
		maximise_combos = sapply(c('high_frequency'), #variable
            function(x) pars$threshold[[x]]) %>% unlist
		cat("CHECK: minimise these combos: ", minimise_combos, "\n",
			append = T, file = f("{runID}.log"))
 		tcell_pct=sapply(1:nrow(mat), function(x)
			sum(mat[x, minimise_combos]) /
				sum(mat[ x, c(minimise_combos, maximise_combos)])
		)
		print(colnames(mat))
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
			Dsel = intersect(which(tcell_pct == min(tcell_pct[which(tcell_count)], na.rm = T)), 
						   which(tcell_count))
						   
			cat('INFO: Dsel index options', Dsel, '\n')
			Dsel = max(Dsel)
			threshold=mat$D[Dsel] # %>% round(., 2)
			plot(mat$D, tcell_pct, pch=19, cex=0.1, ylim = c(0, 1), 
				main = f("{confidence} confidence D cutoff = {threshold}"), 
				xlab = 'D score', bty = "n",
				ylab = 'Ratio between proportions of\nrare to high cell populations')
			abline(v=threshold, color = 'red', lty =2)
			abline(v=mat$D[which(!tcell_count)], 
				col = "#e6e6e6ff")
			points(mat$D, tcell_pct, pch=19, cex=0.1, type='l')

			cat('Selected D cutoff ', mat$D[Dsel], 
				' for confidence ', confidence, '\n')
			cat('INFO: Selected D cutoff ', mat$D[Dsel], 
				' for confidence ', confidence, '\n', append = T,
				file = f("{runID}.log"))
			cat("CHECK: % rare cell types (e.g. DP): ",
				sum(mat[Dsel, grepv("_", pars$threshold[['rare']])]) / 
					sum(mat[Dsel, grep('V1|^D',names(mat), invert=T)]) * 100, '\n',
				file = f("{runID}.log"), append = T)
			cat("CHECK: % rare cell types (e.g. Double positive T cells): ",
				sum(mat[Dsel,  grepv("_", pars$threshold[['rare']])]), 
				sum(mat[Dsel, grep('V1|^D', names(mat), invert=T)]) , '\n',
					file = f("{runID}.log"), append = T)
			cat("CHECK: % rare cell types (e.g. Double positive T cells): ", 
				sum(mat[Dsel, grepv("_", pars$threshold[['rare']])]) / 
					sum(mat[Dsel, grep('V1|^D', names(mat), invert=T)]) * 100, '\n')

			# Plotting the counts for each possible threshold
			groups=list('rare'=rare, 'present'=present, 'combined'=c(rare, present))
			for(group in names(groups)) {
				g <- ggplot(mat, aes(x = D))
				for(combo in groups[[group]])	{
		            g = g + geom_path(aes_string(y= f("{combo} + 1"), 
						color=f('"{combo}"')), size=.7)
				}
				groups[[group]]= groups[[group]] %>% gsub('(.*)', '`\\1`', .)
				string_formula=f(paste(groups[[group]], collapse = "+", sep = '+'), "+ 1")
				g = g + geom_path(aes_string(y=string_formula), color="black", 
												linetype="dashed", size=.7)
				plot = g + 	scale_color_manual(values = unlist(palette$tcellLines)) +
							scale_linetype_manual(values = unlist(palette$tcellLinetype)) +
							cowplot::theme_cowplot() +
							theme(legend.position = 'top', legend.title = element_blank()) +
							guides(color=guide_legend(nrow = 3), linetype=guide_legend(nrow=3)) +
							scale_y_log10() +
							ylab(f('Cell count [ {group} ]'))
				print(plot + geom_vline(xintercept=threshold, linetype='dotted', color = 'black'))
				
			}
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


assign_cluster_positivity <- function(dfExp, clusters, panel, run, runID,
			fdr=0.05, log2=F, magnitude=10**5,
			fromList=NA)  {
				
	effect_field='AUC'
	
	scranPkg=packageVersion("scran")
	if(scranPkg == "1.12.1") effect_field='overlap'
	if(length(unique(clusters)) == 1)
		return(
			data.frame(cluster=rep(NA, length(clusters)), 
					   name=rep(NA, length(clusters)))
		)
	head(dfExp)	%>% print
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
	outDir=f('{runID}_aux/{effect_field}')
	if(! dir.exists(outDir))
		dir.create(outDir, recursive = T)
	compFile=f('{outDir}/pairwise.RData')
	# Not required for the latest version of the scran package
	if(! file.exists(compFile)) {

		print('Running pairwise comparisons')
		# Get probabilities of overlap
		wilcoxTestUp = scran::pairwiseWilcox(
			t(dfExp), 
			as.factor(clusters), 
			direction="up")
		# print(packageVersion("scran"))
		markersUp = scran::combineMarkers(
			wilcoxTestUp$statistics, 
			wilcoxTestUp$pairs, 
			effect.field = effect_field)
		rowNames = rownames(markersUp[[1]])
		# Data frame with probabilities for each pairwise comparison
		upDf = do.call(cbind, lapply(markersUp, function(x) x[rowNames, ]))
		upDf = upDf[, grep(effect_field, colnames(upDf))]
		# Summarise per marker
		upStats = apply(upDf, 1, mean) + apply(upDf, 1, sd)
		upDf = data.frame(upDf, check.names = F)
		upDf = upDf %>% t %>%  as.data.frame
		save(upDf, upStats, file = compFile)
	} else  {
		load(compFile)
	}
  
	print('Assigning')
	assigned=plot_overlaps(upDf, clusterLabels, runID, effect_field=effect_field)
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
		
		positiveLow=determine_threshold(
				assigned, 
				clusterSize, 
				runID, 
				confidence='low',
				greater=F,
				markers=pars$threshold[["markers"]])
				
		positiveHigh=determine_threshold(
				assigned,
				clusterSize,
				runID,
				confidence='high',
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


assign_celltype <- function(names, markers, 
	majorTypes=NULL, region=NULL, major=F, mostLikelyCellType=NULL) {
	
	markerNames=unlist(markers) %>% unique
	cellTypeNames=get_celltypes(markers)
	
	# Consider only high confidence majorType calls
	# unassigned is not. the same as NULL/none
	if(all(majorTypes %in% c('', 'Excluded')))
		majorTypes = NULL
	if(any(majorTypes %in% mostLikelyCellType)) {
			majorTypes[majorTypes %in% mostLikelyCellType] = 'Unassigned'
	}
	if(! is.null(majorTypes))
		majorTypes = majorTypes %>% 
			gsub('Ambiguous .*', 'Ambiguous', .)
	names = paste(majorTypes, names, sep =  "|")
	if(! is.null(region))
		names = paste(region, names, sep =  "=")
	
	clusters = sapply(unique(names), toString)
	markerFrequency = unlist(markers) %>% table
	
	markersList = get_node_list(markers)
	
	celltypes = sapply(clusters, function(name) {
		if(is.na(name))
			name=""
		region = strsplit(name, split="=")[[1]][1]
		if(is.na(region) | region=='NA' | 
			region=='' | region == name |
			grepl("\\|", region))
			region="none"

	    majorType = gsub("^[^=]+=", "", name) %>% 
			gsub("^([^|]+)\\|.*", "\\1", .)
	    if(is.null(majorTypes) | 
			grepl("Excluded", majorType) |
			grepl("\\|", majorType))
			majorType='none'
		
	    name=gsub("^[^|]*\\|", "", name) %>%
			gsub("^[^=]+=", "", .) %>%
			gsub("^[^:]+:([^ ]+).*", "\\1", .)
		
		considerTissueSeg=
			majorType == 'none' &
			region != 'none' & 
			region > 0

	    if(name %in% c("", "|", "NA"))
	      if(majorType == "none" & region == "none") {
	        return("Unassigned")
			} else if(majorType != "none") {
				if(majorType %in% cellTypeNames || grepl("Unassigned|Ambiguous", majorType))
	       	 		return(majorType)
				return(f("{majorType} - Other"))
	      } else if(considerTissueSeg) {
  	        # added if region == 1
	        return("Epithelial cells - Tissue Segmentation")
		}

	    if(name == "Excluded")
			return(name)
	
		# Split, format and select the celltype-specifc marker
	    cellMarkers = strsplit(name, split = "_")[[1]]		
	    cellMarkers = intersect(cellMarkers, markerNames)
		cat(name, majorType, region, cellMarkers, major, '\n')
		
		if(! length(cellMarkers)) {
	      if(majorType != "none") {
			if(majorType %in% cellTypeNames || grepl("Unassigned|Ambiguous", majorType))
       	 		return(majorType)
			return(f("{majorType} - Other"))
		  }
	      if(considerTissueSeg)
			  return("Epithelial cells - Tissue Segmentation")
	      return("Unassigned")
	    }
	    types = sapply(cellTypeNames, function(celltype) {
			cellspecific_markers = markersList[[celltype]]
			all(cellMarkers %in% cellspecific_markers) &
				length(cellMarkers) == length(cellspecific_markers)
		})
	    if(sum(types) == 1) {
			# complete overlap of all markers with the ones specified
			return(names(types)[types])
		}
		# Is the combination of markers specific for one cell type
	    intersect = sapply(cellMarkers, function(.mark1) {
			sel1 = sapply(cellTypeNames, function(celltype)
			  	.mark1 %in% markersList[[celltype]])
			sel1 = cellTypeNames[sel1] %>% unique
			matching = sapply(cellMarkers, function(.mark2) {
				if(.mark1 == .mark2)
					return(0)
				sel2 = sapply(cellTypeNames, function(celltype)
					.mark2 %in% markersList[[celltype]]
				)
				sel2=cellTypeNames[sel2] %>% unique
					if(!length(intersect(sel1, sel2))) return(0)
						# cat(.mark1, .mark2, intersect(sel1, sel2), '\n')
				intersect(sel1, sel2) %>% length
			})
			sum(matching)
	    })
			
		if(all(intersect == length(cellMarkers) - 1) & length(intersect) > 1) {
			overlap=sapply(names(markers), function(celltype) {
				all(cellMarkers %in% markersList[[celltype]])
			})
			if(any(overlap))
				return(names(overlap)[overlap])
		}
		
		# Check if there are specific cell markers for different cell types
		specificity=sapply(cellMarkers, function(.mark) {
			specific=unlist(markers) %in% .mark
			return(sum(specific) == 1)
		})
		if(sum(specificity) > 1) {
	    	specific_cellTypes = sapply(cellMarkers[specificity], function(.mark1) {
				sel1 = sapply(cellTypeNames, function(celltype)
					.mark1 %in% markersList[[celltype]]
				)
				sel1 = cellTypeNames[sel1] %>% unique
			})
			nrSpecficCellTypes = specific_cellTypes %>%
				unlist %>% unique %>% length
			# Consider if the cell belondgs to the tumour region
			if(nrSpecficCellTypes > 1  & considerTissueSeg)
				return("Epithelial cells - Tissue Segmentation")
			if(nrSpecficCellTypes > 1 & grepl('^Ambiguous|^none|Unassigned', majorType)) {
				print(specific_cellTypes)
				if(! major) 
					return('Ambiguous')
				return("Ambiguous")
			}
		}

		assigned=sapply(cellTypeNames, function(celltype) {
			cellspecific_markers=markersList[[celltype]]
			if(length(cellspecific_markers) == 0) 
				return(0)

			weight = sum(sapply(cellspecific_markers, function(marker) {
				if(! marker %in% names(markerFrequency))
						return(0)
				# The more frequently expressed the marker, the less weight it has
				1 / markerFrequency[[marker]] * sum(marker %in% cellMarkers)
			}))
			weight = weight  * sum(cellMarkers %in% cellspecific_markers) / 
						length(cellspecific_markers)
		})
		
		if(all(assigned == 0)) {
			if(considerTissueSeg)
				return("Epithelial cells - Tissue Segmentation")
			if(! majorType %in% c('none', 'mostLikelyCellType')) {
				if(majorType %in% cellTypeNames)
	       	 		return(majorType)
				return(f("{majorType} - Other"))
			}
			return("Unassigned")
		}
		print(sort(assigned))
		if(! major & ! grepl('none|Unassigned|Ambiguous', majorType)) {
			# reduce possible cell types to children of major
			selectChild = sapply(names(assigned), function(x) {
				if(assigned[x] == 0)
					return(F)
				 isChild(tree = markers, parent = gsub(" - Other", '', majorType), child=x)
			})
			if(any(selectChild)) {
				cat('Truncating assigned based on', majorType, 'to',
					 names(assigned)[selectChild], '\n')
				assigned=assigned[selectChild]
			}
			if(length(assigned) == 1) {
				return(names(assigned))
			}
		}
		
		maxas = max(assigned, na.rm  = T)
		assignedCelltype = names(assigned)[which(assigned == maxas)]
		# if  we rely on any specific markers in the absense of majorType, make sure the assigned type has them
		if(sum(specificity) == 1 & grepl('none|Unassigned', majorType))	{
			assignedMarkers = get_celltype_markers(tree=markers, assignedCelltype)
			cat("Checking if assgined markers ", assignedMarkers, " in assigned cellType'n")
			print(names(specificity)[specificity])
			if(! any(names(specificity)[specificity] %in% assignedMarkers))
				return('Ambiguous')
		}
		
		if(length(assignedCelltype) > 1)  {
			if(! majorType %in% c("none"))	{
				if(majorType %in% cellTypeNames || grepl("Unassigned|Ambiguous", majorType))
	       	 		return(majorType)
				return(f("{majorType} - Other"))
			}
			if(considerTissueSeg) 
				return("Epithelial cells - Tissue Segmentation")
			return("Ambiguous Equal score")
		}
		print('Intersect')
		print(intersect)
		if(length(intersect) > 1 & all(intersect == 0) & 
				! majorType %in% c('none', "Unassigned")) {
			if(! majorType %in% c("none") & ! major) {
				if(majorType %in% cellTypeNames || grepl("Unassigned|Ambiguous", majorType))
	       	 		return(majorType)
				return(f("{majorType} - Other"))
			}
			if(considerTissueSeg & ! major)
				return("Epithelial cells - Tissue Segmentation")
			return("Ambiguous Intersect")
		}
		# See if on the subtree of the majorType
		if(! grepl('none|Unassigned', majorType)) {
			cat("CHECK:", majorType, name, '\n')
			if(! isChild(tree=markers,  parent = gsub(" - Other", '', majorType), child=assignedCelltype)) {
				if(! major) 
					return(f('{assignedCelltype}'))
				if(majorType != 'none') {
					if(majorType %in% cellTypeNames || grepl("Unassigned|Ambiguous", majorType))
		       	 		return(majorType)
					return(f("{majorType} - Other"))
				}
			}
		}
		return(assignedCelltype)
	})
	celltypes[match(names, clusters)]
}

