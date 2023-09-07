
# FILE FORMATS
convert_mat2fcs <- function(inMat)  {
  inMatFcs <- new("flowFrame", exprs=as.matrix(inData))
  return(inMatFcs)
}

# INPUTS
get_region_info <- function(panel, cellIDs, featurePattern, regFile) {

	panel=gsub("_.*", "", panel)
	print(regFile)
	if(! file.exists(regFile))
		stop('The sample annotation file does not exists --sample_file')
	regData = read.delim(f(regFile), stringsAsFactors = F)
	
	# Specific for Tx where regData can have images from multiple panels
	panelSelection = grepl(f("^{toupper(panel)}_"), regData$imagename)
	if(any(panelSelection))
		regData=subset(regData, panelSelection)
	cellRegions = get_cellRegionID(cellIDs)
	ind = match(cellRegions, regData$imagename)
	colNames = unique(grep(f("^{featurePattern}$"), colnames(regData)))
	return(regData[ind, colNames])
}

get_tiff_path <- function(cellIDs, inDir, pattern=".tiff") {
 
  sapply(cellIDs, function(cellID) {
    tiff=f("{inDir}/{cellID}{pattern}")
    if(! file.exists(tiff))
		return(NA)
    return(tiff)
  })
}

get_cellRegionID <- function(cellIDs) {
  
  cellRegions=gsub(".txt", "", basename(cellIDs))
  return(cellRegions)
}

get_markers_underscore <- function(markers, positive) {
	
	hasUnderscore = grep("_", markers, value = T) 
	
	if(length(hasUnderscore) & grepl(paste0(hasUnderscore, sep = "|"), positive)) {
		cols = strsplit(positive, split = paste0(c("_", hasUnderscore), sep = "|"))[[1]]
		cols = c(cols, hasUnderscore)
	} else {
		cols = strsplit(positive, split = "_")[[1]]
	}
	cols = cols[cols %in% markers]
	return(cols)
}

# VISUALIZATION UTILS
plot_spatial <- function(X, Y, clusters, exclude.zero=F)  {
  
  # X=Location_Center_X
  # Y=Location_Center_Y
  ind=1:length(X)
  if(exclude.zero) 
    ind=which(clusters > 0)
  plot_scatter(X[ind], -Y[ind], clusters[ind], xlab="X", ylab="Y")
}

plot_scatter <- function(x, y, color, xlab, ylab,legend=T)  {
  
  tmp=data.frame(x=x, y=y, color=color)
  if(!is.numeric(color))
    tmp$color <- gsub("\\.", " ", color)
  plot <- ggplot(tmp,  aes(x=x, y=y, color=color)) +
    geom_point(size=0.00001, alpha=0.2) +
    theme_minimal() + 
	theme(
		axis.text=element_blank(),
		panel.grid=element_blank(),
        legend.position="top",
        legend.title=element_blank()) +
    guides(color=guide_legend(override.aes=list(size=2, alpha=1))) +
    xlab(xlab) + ylab(ylab) 
  if(!legend) {
    plot=plot + guides(color=F)
  }
  if(length(unique(color)) > 100) {
    plot <- plot + scale_color_gradientn( 
      colours=colorRampPalette(rev(brewer.pal(n=11, name="Spectral")))(50))
  } else {
    inputColors=c(color_clusters, tol21rainbow, rainbow(n=length(color)))
    # colValues=WGCNA::labels2colors(color, colorSeq=inputColors)
    plot <- plot + scale_color_manual(values=inputColors)
    # brewer.pal(name="Dark2", n=length(unique(tmp$color))))
  }
  print(plot)
}

# Get info
get_tissue_category <- function(cellIDs, panel, tissAreaDir, class="regional") {
	panel=gsub("_.*", "", panel)
	inFile=f("{tissAreaDir}/cell_info/{panel}_{class}_cell_info.csv")
	print(inFile)
	if(! file.exists(inFile)) {
		if(class == 'background')
			return(data.frame(
					Background=rep(NA, length(cellIDs)), 
					region=NA
					)
			)
		if(class == 'regional')
			return(
				data.frame(
					Background=rep(NA, length(cellIDs)),
					# Immune=rep(NA, length(cellIDs)),
					Stroma=rep(NA, length(cellIDs)),
					Tumour=rep(NA, length(cellIDs)),
					region=NA)
			)
			output=data.frame(rep(NA, length(cellIDs)), rep(NA, length(cellIDs)))
			colnames(output) = c(f('non-{class}'), class)
			return(output)
	}
	
	data=data.table::fread(inFile)
	if(class == "background")
		data$regionID = abs(data$regionID - 255)/255
	
	tmp = subset(data, regionID > 0 & region != "Background")
	categs = tmp[, lapply(.SD, paste1), .SDcols="region", by=c("imagename", "ObjectNumber")]
	
	wideDf = data.table::dcast(data, imagename + ObjectNumber + centerX + centerY ~ region, value.var="regionID")
	if(nrow(categs) > 0 & class == 'regional') {
		wideDf=categs[wideDf, on = .(imagename=imagename, ObjectNumber= ObjectNumber)]
	} else {
		if(class == 'regional') {
			wideDf[, region := lapply(.SD, function(x) ifelse(x > 0, "Tissue", "Background")), 
	           .SDcols=unique(data$region)]
		} else {
			wideDf[, region := lapply(.SD, function(x) ifelse(x > 0, class, f('non-{class}'))), 
	           .SDcols=unique(data$region)]
		}
	}
	ind = match(cellIDs, with(wideDf, paste(ObjectNumber, imagename)))
	categories = c(names(table(data$region)), "region")
	cbind(wideDf[ind, categories, with=F])
}

get_tissue_area <- function(samples, panel, tissAreaDir, pattern="tissue_area.csv",
                            category='Tissue') {
								
  panel=gsub("_.*", "", panel)
  file=file.path(tissAreaDir, paste(tolower(panel), pattern, sep="_"))
  if(! file.exists(file)) file=file.path(tissAreaDir, pattern)
  if(! file.exists(file)) return(NA)
  print(file)
  tarea = read.csv(file)
  print(head(tarea))
  sapply(samples, function(sample)  {
    if(category == 'Tissue') {
      idx=grep(paste0("(Tissue.*|Tumour.*|Stroma.*|Immune.*)", sample, ".tiff"), tarea$Label)  
    } else {
      idx=grep(f("{category}.*{sample}.tiff"), tarea$Label)  
    }
    
    if(length(idx) == 0) return(NA)
    sum(tarea$Total.Area[idx])/ 10 ** 6
  })
}

get_tissue_categ_area <- function(samples, panel, category) {

  file=f(categArea)
  if(! file.exists(file)) return(NA)
  if(category == "") return(NA)
  tarea = read.delim(file)
  sapply(samples, function(sample)  {
    idx=grep(f("{category}_{sample}.tiff"), tarea$Slice)
    if(length(idx) == 0) return(NA)
    tarea$Total.Area[idx]/ 10 ** 6
  })
}

get_filtered_objects <- function(inData, celltypeModelFile) {
	
  if(all(is.na(inData$probability))) {
	  print('WARNING: Skipping stratification by confidence.')
	  return(NULL)
  } else if(! file.exists(f(celltypeModelFile))) {
	  cat('No stratification -> take lower quartile\n')
	  thres=tapply(inData$probability, inData$cellType, function(subset) quantile(subset, .25, na.rm = T))
		labels=names(table(inData$cellType))
		names(thres) = labels
		print(thres)
		values=sapply(1:nrow(inData), 
			function(x) is.na(thres[[inData$cellType[x]]]) |
				inData$probability[x] > thres[[inData$cellType[x]]])
		print(table(values))
		return(values)
 	}
	inData$meanIntensity=log2(to_magnitude(inData$meanIntensity, pars$magnitude) + 1)
	load(f(celltypeModelFile))
	print(table(inData$cellType))
	predictions=predict(model, newdata=inData, type="response", allow.new.levels=T)
	predictions > modelCutoff
}

# Comparison with flow
get_marker_frequency <- function(data, marker, column="positivity") {
	markerPattern=paste(paste0(marker, c("_", "$", " ")), collapse="|")
	values=rep(f("{marker}-"), nrow(data))
	values[grep(markerPattern, data[[column]])] = f("{marker}+")
	return(values)
}

review_cellType_by_major <- function(cellTypes, majorTypes, positivity, cellAssignFile) {
  
  # A function to be developed to consider the hierarchy
  # I a major type is assigned but the cell type is from another brach, reassign to major
  # consider markers that are epxressed by different cell types e.g. CD16 in T cells, CD38 in T cells
  # if in the same branch, assign to subtype (lowest level)
  positivity=gsub('pos:(.*) neg:', '\\1', positivity)
  print(cellAssignFile)
  
  if(! file.exists(cellAssignFile)) {
    return(cellTypes)
  }
  
  guide = read.delim(cellAssignFile, stringsAsFactors = F)
  guide = subset(guide, !is.na(newCellType) & !is.na(newMajorType))
  ind = match(paste(cellTypes, majorTypes, positivity),
         with(guide, paste(OldCellType, OldMajorType, positive)))
  if(! all(is.na(ind)))
    cellTypes[! is.na(ind)] = guide$newCellType[ind[! is.na(ind)]]
  cellTypes=gsub(' - Other', '', cellTypes)
  cellTypes[grepl('^Epithelial cells - Tissue.*', cellTypes)] = 'Epithelial cells'
  return(cellTypes)
}

review_major_by_cellType <- function(cellTypes, majorTypes, positivity, 
		cellAssignFile, subtypeMarkersList, majorMarkersList) {
  
  # Stromal
  if(file.exists(cellAssignFile)) {
   
	  positivity=gsub('pos:(.*) neg:', '\\1', positivity)
	  # majorTypes[which(majorTypes == "Smooth muscle cells")] = "Myofibroblasts"
	  guide=read.csv(f(cellAssignFile), sep='\t', stringsAsFactors = F)
	  guide=subset(guide, !is.na(newCellType) & !is.na(newMajorType))
	  ind=match(paste(cellTypes, majorTypes, positivity),
	         with(guide, paste(OldCellType, OldMajorType, positive)))
	  if(!all(is.na(ind))) {
	    majorTypes[! is.na(ind)] = guide$newMajorType[ind[! is.na(ind)]] 
	  }
  }
  ambiguousIndices = grepl('Ambiguous|Unassigned', majorTypes) & 
  	! grepl('Ambiguous|Unassigned', cellTypes)
	
	if(any(ambiguousIndices)) {
		
	    majorReview = sapply(unique(cellTypes[ambiguousIndices]), function(cellType) {
	  		selectedMarkers = get_celltype_markers(subtypeMarkersList, cellType)
	  		assign_celltype(paste(selectedMarkers, sep = "_", collapse = "_"), majorMarkersList, major = T)
	    })
	    majorTypes[ambiguousIndices] = sapply(which(ambiguousIndices), function(x) {
	  	  cellType = gsub('For Review ', '', cellTypes[x])
	  	  if(cellType %in% c('Ambiguous', 'Unassigned'))
	  		  return(cellType)
		  if(majorReview[x] %in% c("Unassigned"))
			  return(cellType)
		  if(cellType == 'CD4+ Myeloid cells')
			  majorReview[x] = 'CD4+ Myeloid cells'
	  	  majorReview[x]
	    })
	}
	majorTypes[grepl('^Epithelial cells - Tissue.*', majorTypes)] = 'Epithelial cells'
	
	return(majorTypes)
}

summary_table <- function(inDf, pars, regFile) {
	
	print(head(inDf))
  markers = setdiff(unique(unlist(lapply(unique(inDf$positivity),
                                       function(x) strsplit(toString(x), split="_")[[1]]))), c('NA', ''))
									   
  dfMarker = sapply(markers, function(marker) {
    get_marker_frequency(inDf, marker)
  })
  dfMarker=data.frame(imagename=   inDf$imagename,
                      cellCount=   inDf$cellCount,
                      cellDensity= inDf$cellDensity,
                      cellType=    inDf$cellType,
                      dfMarker,
                      majorType=    inDf$majorType, stringsAsFactors = F)

   
  if(length(pars$experimental_condition) > 0 & file.exists(regFile)) {
	  # get info
	  info = sapply(pars$experimental_condition, function(factor) {
		 	 get_region_info(panel = pars$panel, 
		  					cellIDs = inDf$imagename,
							regFile = regFile,
							featurePattern = f("{factor}$"))
	  })
	  colnames(info) = pars$experimental_condition
	  if(any(is.na(info)))
		  cat("WARNING: check whether the imagenames match the sample annotation table. Missing values for some.\n")
	  dfMarker = cbind(dfMarker, info) %>% as.data.frame
	  inDf = cbind(inDf, info)  %>% as.data.frame

	  dfMarkMelt = reshape2::melt(dfMarker, id.vars=c("imagename", pars$experimental_condition, "cellCount", "cellDensity"))
	  ids = c("imagename", pars$experimental_condition)
  } else {
	  dfMarkMelt = reshape2::melt(dfMarker, id.vars=c("imagename",  "cellCount", "cellDensity"))
	  ids = c("imagename")
	  
  }
  dfMarkStats = ddply(dfMarkMelt, c(ids, "variable", "value"), summarise,
                    count = sum(cellCount, na.rm = T),
                    density = sum(cellDensity, na.rm = T))

  # Output for cell density table
  inDf$cellID=paste(inDf$cellType, sep = ":")
  print(head(inDf))
  typeStats = ddply(inDf, .(imagename, cellID), summarise,
                  cellDensity = ifelse(all(is.na(cellDensity)), NA, sum(cellDensity, na.rm = T)),
                  cellCount = sum(cellCount, na.rm = T),
                  cellPercentage = sum(cellPercentage, na.rm = T), .drop=F)

  if(length(pars$experimental_condition) > 0 & file.exists(regFile)) {
	  # get info
	  info = sapply(pars$experimental_condition, function(factor) {
		 	 get_region_info(panel = pars$panel, 
		  					cellIDs = typeStats$imagename,
							regFile = regFile,
							featurePattern = f("{factor}$"))
	  })
	  colnames(info) = pars$experimental_condition
	  typeStats = cbind(imagename = typeStats[, "imagename"], info, typeStats[, ! colnames(typeStats) %in%"imagename"])
  }
	
  posStatsList = lapply(markers, function(marker) {
    print(marker)
    inDf$selection = get_marker_frequency(inDf, marker, 'positivity')
    stats=ddply(inDf, .(imagename, cellID), summarise,
                cellDensity=sum(cellDensity[grep("\\+", selection)], na.rm = T))
    stats$marker=marker
    stats
  })

  posStatsList = do.call(rbind, posStatsList)
  posStats = reshape2::dcast(data = posStatsList, formula = imagename + cellID ~ marker,
                           value.var =  "cellDensity")
  markerColSel = colnames(posStats) %in% markers
  colnames(posStats)[markerColSel] = paste0('cellDensity_', colnames(posStats)[markerColSel])
  
  typeMatch = match(with(typeStats, paste(imagename, cellID)), with(posStats, paste(imagename, cellID)))
  typeStats$majorType = gsub(':.*', "", typeStats$cellID)
  typeStats$cellType = gsub('^[^:]+:', "", typeStats$cellID)
  typeStats = typeStats[, c(ids, 'cellType', 'cellDensity', 'cellCount', 'cellPercentage')]
  typeStats = cbind(typeStats, posStats[typeMatch, markerColSel])
  typeStats[is.na(typeStats)] = 0
  return(typeStats)
}
