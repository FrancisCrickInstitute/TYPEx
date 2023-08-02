
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
	regData=read.delim(f(regFile), stringsAsFactors = F)
	
	# Specific for Tx where regData can have images from multiple panels
	panelSelection = grepl(f("^{toupper(panel)}_"), regData$imagename)
	if(any(panelSelection))
		regData=subset(regData, panelSelection)
	cellRegions=get_cellRegionID(cellIDs)
	ind=match(cellRegions, regData$imagename)
	colNames=unique(grep(f("^{featurePattern}$"), colnames(regData)))
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
  
  cellRegions=gsub("^([^\\.]+).*", "\\1", basename(cellIDs))
  return(cellRegions)
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

overlay_cell_types <- function(imagename, imageDf, cellTypeCol, cellTypes = NULL,
                               title=NULL, legend=F, local=T,
                               pattern = "composite_with_outlines.png",
                               xstart=NULL, xend=NULL, ystart=NULL, yend=NULL, ptx = 1,
                               imageType="outlines") {
  # cellTypes = NULL; title=NULL; legend=F; local=T; run="B"; pattern = "composite_with_outlines.png"
  # xstart=NULL; xend=NULL; ystart=NULL; yend=NULL; ptx = 1
  inDir=file.path(Sys.getenv('BASE_DIR'), 'output/composites')
  rawImgFile=get_tiff_path(cellIDs = imagename, inDir=inDir)
  compositeImgFile=NA 
  if(is.na(rawImgFile)) {
    cat("file not found:", imagename, inDir, "\n")
    outlineImg=NULL
  } else {
	library(tiff)
    imageExtension=gsub(".*\\.([^.]+)$", "\\1", rawImgFile)
    readImg=match.fun(paste0("read", toupper(imageExtension)))
    outlineImg=readImg(rawImgFile)
  }
  if(is.na(compositeImgFile)) {
    cat("file not found:", imagename, compositeImgFile, "\n")
    compositeImg=NULL
  } else {
    imageExtension=gsub(".*\\.([^.]+)$", "\\1", compositeImgFile)
    readImg=match.fun(paste0("read", toupper(imageExtension)))
    compositeImg=readImg(compositeImgFile)
  }
  if(! is.null(cellTypes)) {
    if(cellTypes == "Background") {
      rowSel=imageDf$background.Background == 1
    } else if(cellTypes %in% "filtered") {
      rowSel=! imageDf$filtered
    } else if(cellTypes =="nonfiltered") {
      rowSel=imageDf$filtered
    } else if(cellTypes=="unassigned") {
      unassigned=union(intersect(which(imageDf[[cellTypeCol]] == ""), 
                           grep("Excluded", imageDf$cluster)),
                       which(imageDf[[cellTypeCol]] == 'Unassigned'))
      rowSel=rep(F, nrow(imageDf))
      rowSel[unassigned] = T
    } else if(cellTypes=="assigned") {
      rowSel=imageDf[[cellTypeCol]] != 'Unassigned'
    } else {
      rowSel=!is.na(imageDf[[cellTypeCol]])
    }
  } else {
    rowSel=!is.na(imageDf[, cellTypeCol])
  }
  
  if(!is.null(xstart))
    rowSel = rowSel & imageDf$centerX >= xstart
  if(!is.null(xend))
    rowSel = rowSel & imageDf$centerX <= xend
  if(!is.null(ystart))
    rowSel = rowSel & imageDf$centerY >= ystart
  if(!is.null(yend))
    rowSel = rowSel & imageDf$centerY <= yend
  if(!sum(rowSel, na.rm = T)) return(NULL)
  imageDf=imageDf[which(rowSel), ]
  cellStats=ddply(imageDf, .(get(cellTypeCol)), summarise, count=length(positive))
  cellStats=cellStats[order(cellStats$count, decreasing = T), ]
  imageDf=droplevels(imageDf[ imageDf[[cellTypeCol]] %in% cellStats[,1], ])
  imageDf[[cellTypeCol]]=factor(imageDf[[cellTypeCol]], levels=cellStats[, 1])
  if(imageType == "outlines")  {
    image=compositeImg
  } else if(imageType  == "composite") {
    image=outlineImg
  } else {
    "unknown"
  }
  # plot_image_cell_overlay(imageDf, image=compositeImg, clusters=cellTypeCol,
                          # legend=legend, title=title, ptx=ptx)
  plot_image_cell_overlay(imageDf, image=image, clusters=cellTypeCol,
                          legend=legend, title=title, ptx=ptx)

  # xstart=xstart, yend=yend, xend=xend, ystart=ystart)
  # xstart=xstart, yend=yend, xend=xend, ystart=ystart)
}

plot_image_cell_overlay <- function(cellData, image=NULL, clusters, title=NULL, legend = F,
                                    xstart=NULL, ystart=NULL, xend=NULL,  yend=NULL, ptx=1, 
                                    legendTextSize=7) {
  legend.position = "right"
  nrCategs=length(unique(cellData[[clusters]]))
  if(nrCategs>20) legendTextSize=4
  if(nrCategs ==1 ) legend.position='top'
  # cellData[, ..clusters := sapply(.SD, function(x) toString(x)), .SDcols=clusters]
  if(length(clusters))
  plot=ggplot(data=cellData, aes_string("centerX", "centerY", color=clusters)) +
    coord_equal() +
    scale_y_reverse() +
    xlab("") + ylab("") +
    theme_minimal() +
    theme(legend.position = legend.position,
          legend.background = element_rect(fill="black"),
          axis.line = element_blank(), axis.ticks = element_blank(), 
          axis.text = element_blank(), axis.title = element_blank(),
          panel.background = element_rect(fill="black", colour = NA),
          plot.background = element_rect(fill="black", colour = NA),
          plot.title=element_text(color="white"),
          panel.grid = element_blank(),
          panel.spacing = unit(0, "in"),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit="in"),
          # legend.justification = c(0, 0),
          # legend.key.size = unit(x = 2, units =  "cm"),
          legend.text=element_text(color="white", size = legendTextSize),
          legend.spacing = unit(0.1, "in")) +
    labs(x=NULL, y=NULL) 
  if(clusters == "cellType")
	colors=c(palette$cellTypeColors, palette$cellTypingStatusCols, palette$backgroundCols)
	colors=colors[names(colors) %in% unique(cellData[["cellType"]])]
    plot = plot + scale_color_manual(values=colors)
  if(!is.null(image)) {
    xmin=ifelse(is.null(xstart), 0, xstart)
    xmax=ifelse(is.null(xstart), dim(image)[2], xend)
    ymin=ifelse(is.null(ystart), 0, ystart)
    ymax=ifelse(is.null(ystart), -dim(image)[1], yend)
    plot=plot + annotation_raster(image, xmin = xmin, xmax = xmax, ymin = ymin, ymax=ymax)
   # plot = plot + annotation_raster(image, xmin = 0, xmax = max(imageDf$centerX), ymin =-max(imageDf$centerY), ymax=0)
  }
  if(!is.function(title))
    plot = plot + ggtitle(title)
  if(!legend) {
    plot=plot + guides(colour = legend)
  } else if(nrCategs>20) {
    plot = plot + guides(color = guide_legend(override.aes = list(shape = 15, size =.5),
                                              legend.spacing = unit(0, "in"), ncol = 2))
  } else {
    plot = plot + guides(color = guide_legend(override.aes = list(shape = 15, size = 2), ncol = 1))
  }
  plot = plot + geom_point(size=ptx, stroke=0.1) # shape = 1
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

get_positive_combination <- function(data, combo, column='positive')  {
  
  markers=sapply(combo, function(x)
    strsplit(x, split = '\\+|-')[[1]]) %>%
    unlist %>% unique
  pos = lapply(markers, function(x)
    get_marker_frequency(data, x, column=column))
  pos=do.call(paste, pos) %>%
    gsub(' ', '', .)
  pos==combo
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

summary_table <- function(inDf) {
  markers=setdiff(unique(unlist(lapply(unique(inDf$positivity),
                                       function(x) strsplit(toString(x), split="_")[[1]]))), c('NA', ''))
  dfMarker=sapply(markers, function(marker) {
    get_marker_frequency(inDf, marker)
  })
  dfMarker=data.frame(imagename=   inDf$imagename,
                      cellCount=   inDf$cellCount,
                      cellDensity= inDf$cellDensity,
                      cellType=    inDf$cellType,
                      dfMarker,
                      majorType=    inDf$majorType, stringsAsFactors = F)
  dfMarkMelt=reshape2::melt(dfMarker, id.vars=c("imagename", "majorType", "cellCount", "cellDensity"))
  dfMarkStats=ddply(dfMarkMelt, .(imagename, variable, value), summarise,
                    count=sum(cellCount, na.rm = T),
                    density = sum(cellDensity, na.rm = T))

  # Output for cell density table
  inDf$cellID=paste(inDf$majorType, inDf$cellType, sep = ":")
  typeStats=ddply(inDf, .(imagename, cellID), summarise,
                  cellDensity=sum(cellDensity, na.rm = T),
                  cellCount=sum(cellCount, na.rm = T),
                  cellPercentage=sum(cellPercentage, na.rm = T), .drop=F)
  posStatsList=lapply(markers, function(marker) {
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
  typeStats = typeStats[, c('imagename', 'majorType', 'cellType', 'cellDensity', 'cellCount', 'cellPercentage')]
  typeStats = cbind(typeStats, posStats[typeMatch, markerColSel])
  typeStats[is.na(typeStats)] = 0
  return(typeStats)
}
