get_critical_value <- function(n1, n2, file = NULL, d=10000, p=0.001) {
  set.seed(123)
  if(!is.null(file))  {
    if(file.exists(file)) {
    data=read.csv(file, sep = '\t')
    sel=which(data$n1 == n1 & data$n2 == n2)
    if(length(sel))
      return(data$threshold[sel])
    }
  }
  # d <- 10000
  tab=lapply(1:d, function(i) {
    if(!i %% 10000) print(i)
    x1 <- runif(n1, min=-1, max=1) # + runif(1, min=-1, max=1) 
    x2 <- runif(n2, min=-1, max=1) # + runif(1,  min=-1, max=1) 
    out <- ks.test(x1, x2) #, alternative = 'g'
    c(out$statistic, pval=out$p.value)
  })
  tab=data.frame(do.call(rbind, tab))
  threshold=quantile(tab$D, p)
  if(!is.null(file))
    write.tab(cbind(n1, n2, threshold), file = file)
  return(threshold)
}

get_critical_value_per_marker <- function(x, y, d=1000, p=0.05) {
  set.seed(123)
  # d <- 10000
    tab=lapply(1:d, function(i) {
      if(!i %% 10000) print(i)
        x1=sample(y, size=length(x), replace = F)
        # x1<-runif(length(x), min=0, max=1)
        # x2<-runif(length(y), min=0, max=1)
        # runif(1,  min=-1, max=1)
        out<-ks.test(x = x1, y = y, alternative='g')
        
        # if(out$statistic > 0.03) {
        #   lines(density(x1), lty=2, color='red')
        # } else {
        #   lines(density(x1), lty=2)
        # }
        c(out$statistic, pval=out$p.value)
        # c(med=median(x1))
    })
    # dev.off()
    tab=data.frame(do.call(rbind, tab))
    # threshold=sum(tab$p < p)/d
    # threshold=mean(tab$D)
    # threshold=quantile(tab$D, 1-p)
    # sum(tab$D < vartest$statistic)/nrow(tab)
    # threshold=sum(tab$D < vartest$statistic)/nrow(tab)
    threshold=quantile(tab$D, p)
    return(threshold)
}

get_conf <- function(n, y, d=1000, p=0.05) {
  set.seed(555)
  # d <- 10000
  # plot(ecdf(y))
  tab=lapply(1:d, function(i) {
        if(!i %% 10000) print(i)
        x1=sample(y, size=n, replace = F)
        # x2<-runif(n, min=0, max=1)
        # runif(1,  min=-1, max=1)
        # lines(ecdf(x1), lty=2, col='red')
        out<-ks.test(x = x1, y = setdiff(y, x1), alternative='l', exact = F)
        c(out$statistic, pval=out$p.value)
    })
    # lines(ecdf(y), lty=2, col='red')
    # dev.off()
    tab=data.frame(do.call(rbind, tab))
    # threshold=sum(tab$p < p)/d
    threshold=mean(tab$D)
    # threshold=quantile(tab$D, 1-p)
    # sum(tab$D < vartest$statistic)/nrow(tab)
    # threshold=sum(tab$D.. < vartest$statistic)/nrow(tab)
    # threshold=quantile(tab$D, p)
    return(threshold)
}

# FILE FORMATS
convert_mat2fcs <- function(inMat)  {
  inMatFcs <- new("flowFrame", exprs=as.matrix(inData))
  return(inMatFcs)
}

# FORMATTING
markers2symbols <- function(values) {
  # values=an array of protein names, the names of the ion Abs
  values <- gsub("Vimentin", "VIM", values)
  values <- gsub("Collagen", "COL1A1", values)
  values <- gsub("pancytokeratin", "panCK", values)
  values <- gsub("alphaSMA", "ACTA2", values)
}

markers_format <- function(values)  {
  # values = an array of protein names, the names of the ion Abs
  return(gsub("vimentin", "Vimentin",
              gsub("CASP3", "casp3",
                   gsub("^aSMA.*","alphaSMA",
                        gsub("pancytokeratin", "panCK", values)))))
}

step_backward_lme4_model=function(features, data) {
  
  none=lme4::glmer(as.formula(f("control ~ ", paste(features, collapse='+'),
                                " + (1|cellType)")), data=data, family=binomial)
  print(BIC(none))
  BICvalues=lapply(features, function(feature) {
    model=lme4::glmer(as.formula(f("control ~ ", paste(setdiff(features, feature), collapse='+'),
                                   " + (1|cellType)")), data=data, family=binomial)
    print(BIC(model))
    model
  })
  subset=sapply(BICvalues, BIC) > BIC(none)
  print(subset)
  if(sum(subset) == length(subset)) return(none)
  if(sum(subset) > 0) return(step_backward_lme4_model(features[subset], data))
}

step_forward_lme4_model=function(features, data, model=NULL)  {
  
  if(is.null(model)) {
    BICvalues=sapply(features, function(feature) {
      print(feature)
      model=lme4::glmer(as.formula(f("control ~ ", feature, " + (1|cellType)")),
                        data=data, family=binomial)
    })
  } else {
    selectedFeatures=rownames(summary(model)$vcov)[-1]
    BICvalues=sapply(features, function(feature) {
      print(feature)
      model=lme4::glmer(as.formula(f("control ~ ", paste(c(feature, selectedFeatures), collapse="+"),
                                     " + (1|cellType)")), data=data, family=binomial)
    })
  }
  bics=sapply(BICvalues, BIC)
  if(!is.null(model))
    if(!any(bics < BIC(model))) return(model)
  minBIC=names(bics)[bics==min(bics)]
  step_forward_lme4_model(setdiff(features, minBIC), data, BICvalues[[minBIC]])
}

get_bg_intensities <- function(markers, panel, run, cohort, study) {
  panel=gsub('_.*', "", panel)
  bgFile=f("tissue_segmentation/regional_intensity/range_{panel}_regional_median_MeanIntensity_{run}.txt")
  if(! file.exists(bgFile))
    return(rep(NA, length(markers)))
  bgDf=read.csv(bgFile, sep = "\t")
  bgDf$range[match(f("{markers}_Background"), rownames(bgDf))]
}

get_tissue_category <- function(cellIDs, run, panel, study, cohort, tissAreaDir, class="regional") {
  panel=gsub("_.*", "", panel)
  inFile=f("{tissAreaDir}/cell_info/{panel}_{class}_cell_info_{run}.csv")
  print(inFile)
  if(!file.exists(inFile)) {
    if(class == 'background')
      return(data.frame(Background=rep(NA, length(cellIDs)), region="Background"))
    if(class != 'background')
      return(data.frame(Background=rep(NA, length(cellIDs)),
                   # Immune=rep(NA, length(cellIDs)),
                   Stroma=rep(NA, length(cellIDs)),
                   Tumour=rep(NA, length(cellIDs)),
                   region="Background"))
  }
  data=data.table::fread(inFile)
  if(class=="background")
    data$regionID = abs(data$regionID - 255)/255
  tmp=subset(data, regionID > 0 & region != "Background")
  categs=tmp[, lapply(.SD, paste1), .SDcols="region", by=c("imagename", "ObjectNumber")]

  wideDf=data.table::dcast(data, imagename + ObjectNumber + centerX + centerY ~ region, value.var="regionID")
  if(nrow(categs) > 0) {
    wideDf=categs[wideDf, on = .(imagename=imagename, ObjectNumber= ObjectNumber)]
  } else {
    wideDf[, region := lapply(.SD, function(x) ifelse(x>0, "Background", "Tissue")), 
           .SDcols=unique(data$region)]
  }
	print(head(f("{wideDf$ObjectNumber} {toupper(panel)}_{wideDf$imagename}")))
  print(table(wideDf$imagename))
  ind=match(cellIDs, f("{wideDf$ObjectNumber} {toupper(panel)}_{wideDf$imagename}"))
   print(head(ind))
	print(head(cellIDs))

  categories=c(names(table(data$region)), "region")
  cbind(wideDf[ind, categories, with=F])
}

# INPUTS
get_column_names <- function(colnames, features, markers_select=NULL, channels_exclude=NULL) {
  
  columnPattern=paste(paste0("_", features, "_"), collapse="|")
  if(!is.null(markers_select)) {
    columnPattern=paste(sapply(features, function(f)
      paste0("_", f, "_", markers_select, "_")), collapse="|")
  }
  columnNames=grep(columnPattern, colnames(inData), value=T)
  if(!is.null(channels_exclude)) {
    excludePattern=paste0("_", channels_exclude, "_", collapse="|")
    setdiff(columnNames,
            grep(excludePattern, colnames(inData), value=T))
  }
  return(columnNames)
}

get_region_info <- function(panel, cellIDs, featurePattern, regFile) {
  panel=gsub("_.*", "", panel)
  regData=read.csv(f(regFile), sep="\t", stringsAsFactors = F)
  regData=regData[grep(paste0("^", toupper(panel), "_"), regData$imagename), ]
  # metaData=data.table::fread(metaFile, sep=",")
  # regMatch=match(regData$TMA_ID_actual_forlink, metaData$TMA_position)
  # cat("Unmatched regions:", unique(regData$TMA_ID_actual[is.na(regMatch)]), "\n")
  # regData=cbind(regData, metaData[regMatch, ])
  cellRegions=get_cellRegionID(cellIDs)
  ind=match(cellRegions, regData$imagename)
  colNames=unique(grep(paste0("^", featurePattern, "$"), colnames(regData)))
  return(regData[ind, colNames])
}

get_composite <- function(run, cellIDs, pattern="composite_with_outlines.png", local=T) {
  
    if(!local) {
      panels=unique(gsub("^(P.)_.*", "\\1", cellIDs))
      inDir=file.path(rubiconDir, "outputs/NextFlow") #, */*/mask_outlines/
      runs=unique(grep(paste0("_", run, "$"), list.files(f("{inDir}/{panels}")), value = T))
      maskDirs=lapply(panels, 
                      function(p) grep("key_markers", 
                                       list.files(f("{inDir}/{p}/{runs}/mask_outlines"), full.names = T), value = T))
      names(maskDirs)=panels
    } else {
      maskDirs=lapply(c("P1", "P2"), function(p) f("/Volumes/pdfci/data/rubicon/composites/{p}/{run}") )#f("analyses/typing/review/composites"))
      names(maskDirs)=c("P1", "P2")
    }
   
    sapply(cellIDs, function(cellID) {
      cellID=gsub("-", "_", cellID)
      panel=gsub("^(P.)_.*", "\\1", cellID)
      tiff=list.frec(maskDirs[[panel]], pattern = paste0(".*", cellID, "_.*", pattern))[1]
    })
}

get_tiff_path <- function(cellIDs, pattern=".tiff", local=T, inDir='output/tissue_seg/composites') {
 
  sapply(cellIDs, function(cellID) {
    tiff=file.path(inDir, paste0(cellID, pattern))
    if(!file.exists(tiff)) return(NA)
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
    theme_minimal() + theme(axis.text=element_blank(),
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

#get_duplicates <- function(imagename, regFile) {
  
#  meta=read.csv(f(regFile), sep = "\t")
#  tmaPosition=meta$TMA_ID_corrected[which(meta$imagename==imagename)]
#  tmaSelect=which(meta$TMA_ID_corrected == tmaPosition)
#  replicates=setdiff(meta$imagename[tmaSelect], imagename)
#  return(replicates)
#}

overlay_cell_types <- function(imagename, imageDf, cellTypeCol, cellTypes = NULL,
                               title=NULL, legend=F, local=T, run="B",
                               pattern = "composite_with_outlines.png",
                               xstart=NULL, xend=NULL, ystart=NULL, yend=NULL, ptx = 1,
                               imageType="outlines") {
  # cellTypes = NULL; title=NULL; legend=F; local=T; run="B"; pattern = "composite_with_outlines.png"
  # xstart=NULL; xend=NULL; ystart=NULL; yend=NULL; ptx = 1
  inDir=file.path(Sys.getenv('BASE_DIR'), 'output/composites')
  rawImgFile=get_tiff_path(cellIDs = imagename, local = local, inDir=inDir)
  compositeImgFile=get_composite(run=run, cellIDs = imagename)
  if(is.na(rawImgFile)) {
    cat("file not found:", imagename, run, inDir, "\n")
    outlineImg=NULL
  } else {
	library(tiff)
    imageExtension=gsub(".*\\.([^.]+)$", "\\1", rawImgFile)
    readImg=match.fun(paste0("read", toupper(imageExtension)))
    outlineImg=readImg(rawImgFile)
  }
  if(is.na(compositeImgFile)) {
    cat("file not found:", imagename, run, compositeImgFile, "\n")
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
get_tissue_area <- function(samples, panel, study, cohort, tissAreaDir, pattern="tissue_area.csv",
                            category='Tissue') {
  panel=gsub("_.*", "", panel)
  file=file.path(tissAreaDir, paste(study, cohort, tolower(panel), pattern, sep="_"))
  if(!file.exists(file)) file=file.path(tissAreaDir, pattern)
  if(!file.exists(file)) return(NA)
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
  if(!file.exists(file)) return(NA)
  if(category == "") return(NA)
  tarea = read.csv(file, sep = "\t")
  sapply(samples, function(sample)  {
    idx=grep(paste0(category, "_", sample, ".tiff"), tarea$Slice)
    if(length(idx) == 0) return(NA)
    tarea$Total.Area[idx]/ 10 ** 6
  })
}

get_filtered_objects <- function(inData, subset, run, markers, panel, cohort, celltypeModelFile) {
  if(! file.exists(f(celltypeModelFile))) return(NA)
  if(all(is.na(inData$probability))) return(NA)
  load(f(celltypeModelFile))
  print(table(inData$cellType))
  #return(inData$probability > 0.7)
  predictions=predict(model, newdata=inData, type="response", allow.new.levels=T)
  # predictions=predict(model, newdata=inData, type="response")
  predictions > modelCutoff
}

# get_TMA_areas <- function(inDir,
#                      feature="LocationCenter",
#                      x="LocationCenter_X",
#                      y="LocationCenter_Y")  {
#   
#   inData=load_files(inDir, feature)
#   samples=unique(gsub(".txt", "", basename(inData$file)))
#   test=sapply(samples, function(sample) {
#     xrange=range=c(min(inData[, x], na.rm=T), max(inData[, x], na.rm=T))
#     yrange=c(min(inData[,y], na.rm=T), max(inData[, y], na.rm=T))
#     # assuming circular shape, rough area estimation
#     diam=max(c(xrange[2]-xrange[1]), (yrange[2] - yrange[1]))
#     (diam/2000) ** 2 * pi
#   })
# }

# get_stack_metals <- function(panel, run, stack)  {
#   
#   panelDir=file.path(rubiconDir, "outputs", "NextFlow", panel)
#   runs=list.files(panelDir, pattern=paste0("_", run), full.names=T)
#   metaData=file.path(runs[1], "mcd/metadata.csv") %>% read.csv
#   sapply(metaData$metal[metaData[stack] == 1], toString)
# }

get_stack_markers <- function(panel, run, stack="ilastik_stack")  {
  
  panelDir=file.path(rubiconDir, "outputs", "NextFlow", panel)
  runs=list.files(panelDir, pattern=paste0("_", run), full.names=T)
  stacks=list.files(file.path(runs, "results/imctools"),
                      pattern=stack, include.dirs=T, 
                      recursive=T, full.names=T)
    markers=gsub("[^_]+_([^.]+).tiff", "\\1", list.files(path=stacks))
  counts=table(markers)
  markers=names(counts)[counts == max(counts)]
  markers=gsub("^aSMA", "alphaSMA",
               gsub("vimentin", "Vimentin",
                    gsub("CASP3", "casp3",
                         gsub("pancytokeratin", "panCK", markers))))
  markers %>% return
}

assign_celltype<-function(names, markers, majorTypes=NULL, region=NULL, major=F) {

  # prior to calling, names should have been processed by filter_positivity
  # Exclude also nonspecific markers here e.g. FAP1 if major types are immune?
  if(all(majorTypes %in% c('', 'Excluded'))) majorTypes=NULL
  # if(major) markers=marker_gene_list$mcsanov
  # major: celltype pos: neg:
  if(!is.null(majorTypes)) names=paste(majorTypes, names, sep =  "|")
  if(!is.null(region)) names=paste(region, names, sep =  "=")
  markerFrequency=table(unlist(markers))
  clusters=sapply(unique(names), toString)
  celltypes=sapply(clusters, function(name) {
    if(is.na(name)) name=""
    region=strsplit(name, split="=")[[1]][1]
    if(length(grep("\\|", region))) region="none"
    if(is.na(region) | region=='NA' | region=='') region="none"
    if(region == name) region='none'
    
    majorType=gsub("^([^|]+)\\|.*", "\\1", gsub("^[^=]+=", "", name))
    if(is.null(majorTypes)) majorType='none'
    if(length(grep("Excluded", majorType))) majorType="none"
    if(length(grep("\\|", majorType))) majorType="none"
    
    name=gsub("^[^|]+\\|", "", gsub("^[^=]+=", "", name))
    name=gsub("^[^:]+:([^ ]+).*", "\\1", name)
    
    if(name %in% c("", "|", "NA"))
      if(majorType == "none" & region == "none") {
        return("Unassigned") 
      } else if(majorType != "none" & majorType != 'Epithelial cells') {
        return(majorType)
      } else if(region != 'none' & !is.na(region) & region > 0) {
        # added if region == 1
        return("Epithelial cells Region")
      } else if(majorType != 'none') {
		return(majorType)
	}
    
    if(name=="Excluded") return(name)
    cellMarkers=strsplit(name, split="_")[[1]]
    cellMarkers=markers_format(cellMarkers)
    cellMarkers=intersect(cellMarkers, unlist(markers))
    # if(majorType != "none") cellMarkers=setdiff(cellMarkers, nonspecificAb[[majorType]]) 
    if(!length(cellMarkers)) {
	
	  if(majorType == 'Epithelial cells' & region != 'none' & region != name & region > 0) return("Epithelial cells Region")
      if(majorType != "none") return(majorType)
      if(major & majorType != 'none') return('Unassigned')
      if(region != 'none' & region != name & region > 0) return("Epithelial cells Region")
      if(majorType == 'Mesenchymal cells') return(majorType)
      return("Unassigned")
    }
    
    types=sapply(names(markers), function(celltype) {
      all(cellMarkers %in% markers[[celltype]]) &
        length(cellMarkers) == length(markers[[celltype]])
    })
    if(sum(types) == 1)  return(names(types)[types])
    
    intersect=sapply(cellMarkers, function(.mark1) {
      # if(.mark1 == "alphaSMA") return(1)
      sel1=sapply(names(markers), function(celltype) .mark1 %in% markers[[celltype]])
      sel1=names(sel1)[sel1]
      sum(sapply(cellMarkers, function(.mark2) {
        if(.mark1==.mark2) return(0)
        sel2=sapply(names(markers), function(celltype) .mark2 %in% markers[[celltype]])
        sel2=names(sel2)[sel2]
        length(intersect(sel1, sel2))
      }))
    })
    
    if(all(intersect==length(cellMarkers) - 1) & length(intersect) > 1) {
      overlap=sapply(names(markers), function(celltype) {
        all(cellMarkers %in% markers[[celltype]])
      })
      if(any(overlap)) return(names(overlap)[overlap])
    }
    specificity=sapply(cellMarkers, function(.mark) {
      # Vimentin is the only marker that distinguishes uniquely 
      # a cell type but is found in others, EMT
      if(.mark == "Vimentin" & !major) return(F)
	  specific=unlist(markers) %in% .mark

      #specific=sapply(names(markers), function(celltype) {
      #  all(markers[[celltype]] == .mark)
      #})
      return(sum(specific) == 1)
    })
    
	if(sum(specificity) > 1 & majorType == 'none' & region != 'none' & region != name & region > 0) return("Epithelial cells Region")
    if(sum(specificity) > 1 & majorType == 'none') return("Ambiguous")
    
    overlap=sapply(names(markers), function(celltype) {
      sum(cellMarkers %in% markers[[celltype]])
    })
    assigned=sapply(names(markers), function(celltype) {
      if(length(markers[[celltype]]) == 0) return(0)
      sum(sapply(markers[[celltype]], function(marker) {
        if(! marker %in% names(markerFrequency)) return(0)
        1 / markerFrequency[[marker]] * sum(marker %in% cellMarkers) *
          sum(cellMarkers %in% markers[[celltype]])/length(markers[[celltype]])
      }))
    })
    sort(assigned)
    maxas=max(assigned, na.rm  = T)
    # if(sum(assigned > 0) == 1) 
      # return(names(assigned)[assigned>0])
    if(sum(assigned==maxas & maxas > 0, na.rm = T) > 1) {
      if(majorType != "none") return(majorType)
      if(region != 'none' & region != name & region > 0) return("Epithelial cells Region")
      return("Ambiguous")
    }
    if(all(assigned == 0) & majorType != 'none') {
      if(region != 'none' & region != name & region > 0) return("Epithelial cells Region")
      if(majorType != 'Mesenchymal cells') return("Unassigned")
      return(majorType)
    }
    if(length(intersect) > 1 & all(intersect == 0) & majorType != 'none') {
	  if(majorType != 'none' & region != 'none' & region != name & region > 0)
		return("Epithelial cells Region")
      if(majorType != "none") return(majorType)
        return("Ambiguous")
    }
    if(major & sum(assigned==maxas, na.rm = T) > 1 & majorType != "none") return(majorType)
    if(major & sum(assigned==maxas, na.rm = T) > 1 & majorType == "none") return('Unassigned')
    names(assigned)[which(assigned==maxas)]
   })
  celltypes[match(names, clusters)]
}

get_pixel_prop <- function(cellIDs, run, panel, markers, 
                           feature="prop_over_threshold", fun="max") {
  
  pixFile=list.frec(f(pixelDir), pattern=f("Run_{run}_{panel}_allmarkers_overmean_spillovermax_hpr.fst"))
  if(!length(pixFile)) return(NA)
  if(!file.exists(pixFile)) return(NA)
  pixData=fst::read_fst(pixFile)
  pixCellIDs=with(pixData, paste(object, imagename))
  pixMatch=match(cellIDs, pixCellIDs)
  if(all(is.na(pixMatch))) stop('Cell ID format:  "<cellObjectID> <imageID>"')
  pixFlt=pixData[pixMatch, grep(feature, colnames(pixData))]
  colnames(pixFlt)=gsub(paste0("^([^_]+)_.*", feature), "\\1", colnames(pixFlt)) %>% markers_format
  if(nrow(pixFlt) == 0) return(rep(NA, length(cellIDs)))
  if(ncol(pixFlt) == 0) return(rep(NA, length(cellIDs)))
  # Format markers
  markers=gsub("pos:(.*) neg:.*", "\\1", markers)
  markers=gsub("up:(.*) down:.*", "\\1", markers)
  getPixSummary=match.fun(fun)
  cols=strsplit(markers[1], split = "_")[[1]]
  if(length(unique(markers)) & length(cols) == 1 & all(cols != "")) {
      return(pixFlt[, cols])
    } else {
      summary=tapply(1:nrow(pixFlt), markers, function(subset) {
        cluster=markers[subset[1]]
        if(cluster == "") {
          values=rep(NA, length(subset))
        } else {
          cols=strsplit(cluster, split = "_")[[1]]
          if(any(! cols %in% colnames(pixFlt))) {
            values=rep(NA, length(subset))
          } else if(length(cols) == 1) {
            values=pixFlt[subset, cols]
          } else if(length(cols) > 0) {
            if(fun=="paste") {
              values=apply(pixFlt[subset, cols], 1, getPixSummary, collapse="_")
            } else {
              values=apply(pixFlt[subset, cols], 1, getPixSummary)
            }
          }
        }
        names(values)=cellIDs[subset]
        return(cbind(values))
      })
      if(length(summary) == 1) return(summary)
      summary=do.call(rbind, summary)
      return(summary[match(cellIDs, rownames(summary)), 1])
  }
}


filter_positivity <- function(markers, celltypes, panel) {
  if(is.null(celltypes)) return(markers)
  panel=gsub("_.*", "", panel)
  if(!length(markers)) return(NULL)
  abReview<-read.csv(f(abReviewFile), sep = "\t", row.names = 1, check.names = F)
  # exclude the marker if it's considered nonspecific for a certain celltype
  sapply(1:length(markers), function(idx) {
    if(! celltypes[idx] %in% colnames(abReview)) return(markers[idx])
    values=strsplit(markers[idx], split = "_")[[1]]
    values = setdiff(values, rownames(abReview)[abReview[, celltypes[idx]] > 0])
    paste(values, collapse = "_")
  })
}

typing_merger <- function(inFiles, inDir, panel="P1", 
                          run="all", markers="all", method="all",
                          subset="full", decimal=6, assignCellType=T) {

  inDirPattern=analysisPath
  value.var=c("cluster", "names")
  if(markers=="all")  {
    markers_list=marker_gene_list[["majoril"]]
  } else {
    markers_list=marker_gene_list[[args$markers]]
  }
  if(method == "all")
    inDirPattern=gsub(".method.", ".*", inDirPattern)
  if(run == "all")
    inDirPattern=gsub(".run.", ".*", inDirPattern)
  inFiles=grep(f(inDirPattern), inFiles, value = T)  
  df_list=vector(mode="list")
  for(file in inFiles)  {
    export=data.table::fread(file)
    if(assignCellType) {
      if(method == "cellassign") {
        export$cellType=export$cluster
      } else {
        export$cellType=assign_celltype(export$names, markers_list)
      }
      value.var=c(value.var, "cellType")
    }
    df_list[[file]]= export
    cat("Processed file", file, dim(df_list[[file]]), "\n")
  }
  cat("Loaded", length(df_list), "\n")
  if(length(df_list) == 0) 
    return(NULL)
  # panel run markers method pars 
  df_merge=data.table::rbindlist(df_list, fill = T)
  file_base=sapply(df_merge$file, basename)
  df_merge$panel=panel
  df_merge$file =gsub("(.*).txt", "\\1", file_base)
  df_merge$names=gsub("up:(.*) down:.*", "\\1", df_merge$names)
  df_merge$centerX=round(df_merge$centerX, decimal)
  df_merge$centerY=round(df_merge$centerY, decimal)
  # save.image(rdataOut)
  df_per_cell<-data.table::dcast(df_merge,
                                 panel + file + centerX + centerY + object ~ runID,
                                 value.var=value.var)
  return(df_per_cell)
}


# Comparison with flow
get_marker_frequency <- function(data, marker, column="positivity") {
  markerPattern=paste(paste0(marker, c("_", "$", " ")), collapse="|")
  values=rep(f("{marker}-"), nrow(data))
  values[grep(markerPattern, data[[column]])] = f("{marker}+")
  return(values)
}

#get_image_annotation <- function(images) {
#  
#  subsetReview=read.csv(celltypeReviewFile, sep = "\t", stringsAsFactors = F)
#  subsetReview$control[match(images, subsetReview$imagename)]
#}

load_typing_info <- function(inDir, runs, methods, panel, 
                             markers="mcsa", info="Positivity",
                             resultPattern="resultPattern") {
  data=vector(mode="list")
  for(run in strsplit(runs, split=",")[[1]])  {
    for(method in strsplit(methods, split=",")[[1]]) {
        analysisID=paste(method, run, panel, sep="_")
        files=list.files(inDir, pattern=f("{info}.*_{analysisID}.txt"), full.names = T)
        if(!length(files)) next
        for(file in files) {
          params=gsub(f("{info}.{resultPattern}_(.*).{analysisID}.txt"), "\\1", basename(file))
          tmp=read.csv(file, sep = "\t", stringsAsFactors = F)
          data[[paste(run, method, params)]] = cbind(tmp, method=method, params=params, run=run)
        }
    }
  }
  data=do.call(rbind, data)
  return(data)
}


review_cellType_by_major <- function(cellTypes, majorTypes, positivity, panel, study, cohort, cellAssignFile) {
  
  # A function to be developed to consider the hierarchy
  # I a major type is assigned but the cell type is from another brach, reassign to major
  # consider markers that are epxressed by different cell types e.g. CD16 in T cells, CD38 in T cells
  # if in the same branch, assign to subtype (lowest level)
  # CD56_CD79a_panactin_VISTA is assigned to  Lymphocytes
  positivity=gsub('pos:(.*) neg:', '\\1', positivity)
	print(cellAssignFile)
  if(!file.exists(cellAssignFile)) {
    return(cellTypes)
  }
  
  guide=read.csv(cellAssignFile, sep='\t', stringsAsFactors = F)
  guide=subset(guide, !is.na(newCellType) & !is.na(newMajor))
  ind=match(paste(cellTypes, majorTypes, positivity),
            with(guide, paste(cellType, majorType, positivity)))
  if(!all(is.na(ind))) {
    cellTypes[!is.na(ind)] = guide$newCellType[ind[!is.na(ind)]]
  }
  cellTypes=gsub('^CD4$', 'CD4 T cells', cellTypes)
  cellTypes=gsub('^CD8$', 'CD8 T cells', cellTypes)
  cellTypes=gsub('^T cells$', 'T cells - Other', cellTypes)
  cellTypes=gsub('^Leukocytes$', 'Leukocytes - Other', cellTypes)
  cellTypes=gsub('^Smooth muscle cells$', 'Myofibroblasts', cellTypes)
  cellTypes[grepl('^Epithelial cells Region$', majorTypes)] = 'Epithelial cells'
  return(cellTypes)
}

review_major_by_cellType <- function(cellTypes, majorTypes, positivity, panel, study, cohort, cellAssignFile) {
  # Stromal
  if(!file.exists(f(cellAssignFile))) {
   return(majorTypes)
  }
  positivity=gsub('pos:(.*) neg:', '\\1', positivity)
  # majorTypes[which(majorTypes == "Smooth muscle cells")] = "Myofibroblasts"
  guide=read.csv(f(cellAssignFile), sep='\t', stringsAsFactors = F)
  guide=subset(guide, !is.na(newCellType) & !is.na(newMajor))
  ind=match(paste(cellTypes, majorTypes, positivity),
        with(guide, paste(cellType, majorType, positivity)))
  if(!all(is.na(ind))) {
    majorTypes[!is.na(ind)] = guide$newMajor[ind[!is.na(ind)]]
  }
  majorTypes=gsub('^CD4$', 'CD4 T cells', majorTypes)
  majorTypes=gsub('^CD8$', 'CD8 T cells', majorTypes)
  majorTypes=gsub('^T cells$', 'T cells - Other', majorTypes)
  majorTypes=gsub('^Leukocytes$', 'Leukocytes - Other', majorTypes)
  majorTypes=gsub('^Smooth muscle cells$', 'Myofibroblasts', majorTypes)
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
  posStatsList=do.call(rbind, posStatsList)
  posStats=reshape2::dcast(data = posStatsList, formula = imagename + cellID ~ marker,
                           value.var =  "cellDensity")
  markerColSel=colnames(posStats) %in% markers
  colnames(posStats)[markerColSel] = paste0('cellDensity_', colnames(posStats)[markerColSel])
  typeMatch=match(with(typeStats, paste(imagename, cellID)), with(posStats, paste(imagename, cellID)))
  typeStats$majorType=gsub(':.*', "", typeStats$cellID)
  typeStats$cellType=gsub('^[^:]+:', "", typeStats$cellID)
  typeStats=typeStats[, c('imagename', 'majorType', 'cellType', 'cellDensity', 'cellCount', 'cellPercentage')]
  typeStats=cbind(typeStats, posStats[typeMatch, markerColSel])
  typeStats[is.na(typeStats)] = 0
  return(typeStats)
}
