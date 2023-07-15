# Generic function
run_method <- function(inData, method, pars, runID, wDir, regFile, nfDir,
                     feature="MeanIntensity", colors=NULL)  {
						 
  # inData     - a data frame object with cells in rows and markers	 in columns
  # method     - accepts the value of Rphenograph, X-shift, CellAssign, kmeans FlowSOM
  # feature    - the value that summarizes the intensity signal (mean, integrated, median)
  # parameters - method-specific parameters + contains markers to be excluded/included
  
  ids=with(inData, paste(ObjectNumber, basename(imagename)))
  
  funForward=match.fun(paste0("run_", method))
  pars[[method]]$run_id=runID
  pars[[method]]$markers=pars$markers
  pars[[method]]$magnitude=pars$magnitude

  # Load info
  areaData=load_files(file.path(nfDir, "AreaShape"), run=pars$run)
  areaMatch=match(ids, with(areaData, paste(ObjectNumber, basename(imagename))))
  
  # if the analysis has been done
  inMatFile=f("{runID}_aux/inMat.RData")
  clFile=f("{runID}_aux/clusters.RData")
  
  # Select the columns with the selected feature, .e.g MeanIntensity
  columnNames=colnames(inData)
  if(grepl(feature, colnames(inData)))
	  columnNames=grepv(feature, columnNames)
  columnNames =setdiff(columnNames, c(pars$channels_exclude, "ObjectNumber", "imagename"))

  metals = columnNames %>%
    gsub(paste0(".*", feature, "_?_(.*)"), "\\1", .) %>%
	gsub('^([0-9]+)([A-Za-z]+)_.*', '\\2\\1', .)
  metalPattern =paste0(metals, collapse = ".*|") %>%
	paste0(., ".*")
  colnames(inData)=colnames(inData) %>%
	gsub(paste0(".*", feature, "_?_(.*)"), "\\1", .) %>%
	gsub('^[0-9]+[A-Za-z]+_(.*)', '\\1', .) %>%
	gsub(metalPattern, '', .)
  columnNames=columnNames %>%
  	gsub(paste0(".*", feature, "_?_(.*)"), "\\1", .) %>%
	gsub('^[0-9]+[A-Za-z]+_(.*)', '\\1', .) %>%
	gsub(metalPattern, '', .) # remove metals

  if(! file.exists(inMatFile))	{
    # Keep only numeric columns
    rowsKeep=rep(TRUE, nrow(inData))

    if(! is.null(pars$markers) & pars$subset != subtypesDir) {
		selectedMarkers=marker_gene_list[[pars$markers]] %>% 
			unlist %>% unique
		indices=match(selectedMarkers, colnames(inData))
		keep=colnames(inData)[indices[! is.na(indices)]]
		cat("INFO: Column names matching the markers list:", keep, "\n",
			file=f("{runID}.log"), append=T)
	  columnNames=intersect(columnNames, keep)
    }
	cat("Columns:", columnNames, '\n')
	
    # Exlude cells with an area smaller or equal to pars$area_exclude
    if(! is.null(pars$area_exclude) & pars$area_exclude > 0) {
      cat("INFO: Number of cells with area smaller than", pars$area_exclude, "is",
          sum(areaData$AreaShape_Area[areaMatch] <= pars$area_exclude), "\n",
          file=f("{runID}.log"), append=T)
      rowsKeep=rowsKeep & areaData$AreaShape_Area[areaMatch] > pars$area_exclude
    }
	cat('Area', sum(rowsKeep), '\n')

    # Exclude cells with intensity smaller or equal to pars$min_expression for all markers
    rowsKeep=apply(inData[, ..columnNames], 1, 
		function(x) any(x > pars$min_expresssion) & ! any(is.na(x)))
	cat('Expression', sum(rowsKeep), '\n')

	cat("INFO: Number of cells with intensity higher than", pars$min_expresssion, 
		" for all markers of interest",
        sum(! rowsKeep), "\n", file=f("{runID}.log"), append=T)
    
    # Include batches for cellassign; # rows must match cov matrix
    if(method == "cellassign") {
      # add batch effects
	  tma_values = sapply(pars$batch_effects, function(factor) {
		  get_region_info(panel=pars$panel, 
		  				cellIDs=inData$imagename,
						regFile=regFile,
						featurePattern=f("{factor}$"))
	  })
      if(all(is.na(tma_values))) {
		tma_values=rep("one", nrow(inData))
		cat("WARNING:  Batch effects will not be calculated for",
			"the probabilistic cellassign model",
	        sum(! rowsKeep), "\n", file=f("{runID}.log"), append=T)
        #stop("ERROR: Missing TMA values")
	  } else {
	  	tma_values=apply(tma_values, 1, paste, collapse = '_')
  	  }

	  uniqueTMAvalues = unique(tma_values[!is.na(tma_values) & rowsKeep])
	  print(uniqueTMAvalues)
      if(length(uniqueTMAvalues) > 1) {
		  print(head(tma_values))
	  	
		print("TMA values for batch effect correction")
      	rowsKeep=rowsKeep & ! is.na(tma_values)  # & !is.na(tissue_values)
        pars[[method]]$X = model.matrix(~ 0 + tma_values[rowsKeep])   
		cat('Batch', sum(rowsKeep), '\n')
      } else {
		  cat('WARNING: Batch effects will not be calculated for',
		  	"the probabilistic cellassign model\n",
			file=f("{runID}.log"), append=T)
		  pars[[method]]$X=NULL
      }
      # create_covariate_matrix(list(tma_values[rowsKeep], tissue_values[rowsKeep]))
    }
    
    if(! dir.exists(dirname(runID)))
      dir.create(dirname(runID), recursive=T)
	
    cat("Number of cells below intensity and area cutoffs", 
		pars$markers, " is ", sum(! rowsKeep), "\n", append=T)
	cat('Restricting to', sum(rowsKeep), 'rows and ', length(columnNames), 'columns\n')
	inMat=inData[rowsKeep, ..columnNames]
	
    colors=colors[rowsKeep]
	stopifnot(ncol(inMat) & nrow(inMat))
  	cat('INFO: range of values in the input matrix is ', range(inMat, na.rm = T), 
  	    	'using magnitude ', pars$magnitude, '\n')
	
	if(pars[[method]]$transformation != "none") {
		transformInput=match.fun(paste0("to_", pars[[method]]$transformation))
		cat('INFO: Transforming the matrix using -', pars[[method]]$transformation, 
			'- as specified in typing_params.json for the corresponding tool.\n')
			
		if(pars[[method]]$transformation == 'magnitude') {
	        inMat=transformInput(inMat, dynamic_range = pars$magnitude)
		} else {
	        inMat=transformInput(inMat)
		}
    }
	runOutput=rep("Excluded", nrow(inData))
    rownames(inMat) = ids[rowsKeep]
    cat("INFO: Markers for analysis:\n", colnames(inMat), "\n", 
        file=f("{runID}.log"), append=T)
	if(method == "MC")  {
      inMat=apply(inMat, 2, function(x) tapply(x, as.factor(colors), median))
      metaclusters=funForward(inMat, pars[[method]], colors=colors)
      runOutput[rowsKeep]=metaclusters[match(colors, rownames(inMat))]
    } else {
      runOutput[rowsKeep]=funForward(inMat, pars[[method]], colors=colors)
    }
  } else {
    load(inMatFile)
    rowsKeep=match(rownames(inMat), ids)
	
    columnNames=match(colnames(inMat), colnames(inData))

    runOutput=rep("Excluded", nrow(inData))
    if(file.exists(clFile)) {
      load(clFile)
      rowsKeep=rowsKeep[which(!is.na(clusters))]
      clusters=clusters[!is.na(clusters)]
	  cat(length(runOutput), length(rowsKeep), length(clusters), '\n')
      runOutput[rowsKeep] = clusters
    } else {
      runOutput[rowsKeep] = funForward(inMat, pars[[method]], colors=colors)
    }
  }
  cat('WARNING: range of transformed values is ', range(inMat),
  	'using magnitude ', pars$magnitude, '\n',
  	file = f('{runID}.log'), append = T)
  if(method %in% pars$visualisation_methods)
	  return(NULL)
  return(list('runOutput'=runOutput, 'columnNames'=columnNames, 'ids'=ids))
}

# CellAssign
run_cellassign <- function(inMat, p, colors)  {
  
  # Note: marker_gene_list requires an object of class matrix as input
  # marker_mat rows must match inmat columns
	outDir=f("{p$run_id}_aux")
  if(! dir.exists(outDir))
		dir.create(outDir)
  outFile=f("{outDir}/out.RData")
  inMatFile=f("{outDir}/inMat.RData")
  probsFile=f("{outDir}/probs.txt")
  clFile=f("{outDir}/clusters.RData")
  
  if(! file.exists(outFile)) {
	  
	marker_list=get_node_list(marker_gene_list[[p$markers]])
    markerMat=cellassign::marker_list_to_mat(marker_list)
    markerMat=markerMat[rownames(markerMat) %in% colnames(inMat), ]
    markerMat=markerMat[, apply(markerMat, 2, function(x) ! all(x == 0))]
    # if same signature for multiple cells because of missing markers, keep the higest in lineage
    cellcode=apply(markerMat, 2, paste, collapse="")
    exclude=unlist(sapply(unique(cellcode[duplicated(cellcode)]), function(x) {
      duplcates=names(cellcode)[cellcode==x]
      count=sapply(duplcates, function(celltype) {
        length(marker_gene_list[[p$markers]][[celltype]])
      })
      names(count)[count > min(count)]
    }))
    if(length(exclude) > 0)
      markerMat=markerMat[, ! colnames(markerMat) %in% exclude]
    write.table(markerMat, append = T, file=f("{p$run_id}.log"))
    print(markerMat)
    cat("Using", nrow(markerMat), "markers and ", ncol(markerMat),  "celltypes.\n")
    keep=colnames(inMat)[match(rownames(markerMat), colnames(inMat))]
    cat("From them,", length(keep), "are found in inMat.\n")
    
	# Constant library size
    s=rep(1, nrow(inMat))
    fit <- cellassign::cellassign(
		exprs_obj        = as.matrix(inMat[, ..keep]),
        marker_gene_info = markerMat,
        s                = s, 
		X                = p$X,
        learning_rate    = p$learning_rate,
        shrinkage        = p$shrinkage,
        verbose          = p$verbose)
    save(fit, file = outFile)
    save(inMat, file=inMatFile)
  } else {
    load(outFile)
    load(inMatFile)
  }
  probs=cellassign::cellprobs(fit)
  probs=cbind(imagename=rownames(inMat), probs)
  write.tab(probs, file  = probsFile)
  clusters=fit$cell_type
  save(clusters, file = clFile)
  return(clusters)
}

# Rphenograph
run_Rphenograph <- function(inMat, p, colors) {

  library(Rcpp)	
	
	outDir=f("{p$run_id}_aux")
  if(! dir.exists(outDir))
		dir.create(outDir)
		
  outFile=f("{outDir}/out.RData.gz")
  inMatFile=f("{outDir}/inMat.RData")
  probsFile=f("{outDir}/probs.txt")
  clFile=f("{outDir}/clusters.RData")
 
	
  pdfOut=paste0(p$run_id, ".heatmap.pdf")
  if(file.exists(clFile)) {
    load(clFile)
    return(clusters)
  }
  nrClusters=length(unique(colors[!is.na(colors)]))
  print(table(colors))
  if(! file.exists(outFile) | ! file.exists(inMatFile)) {
    if(nrClusters == 1) {
		# From the Rphenograph code: (k > nrow(data)-2)
	  if(nrow(inMat) - 2 < p[["k"]]) return(NA)
      Rout=Rphenograph::Rphenograph(inMat, k=p[["k"]])
    } else {
        Rout=tapply(1:nrow(inMat), colors, function(subset) {
          cluster=colors[subset[1]]
		  cat(cluster, ' has ', sum(colors == cluster), p[["k"]], '\n')
          outCellFile=f("{outDir}/{cluster}.cell.out.RData.gz")
          if(file.exists(outCellFile))  {
            load(outCellFile)
            return(out)
          }
          cat("\n", cluster, length(subset), "\n")
          # Minimum group size for clustering 
          # EDIT: if(length(subset) <= p[["k"]] * 150) return(NA)
		  if(length(subset) - 2 < p[["k"]]) return(NA)	  
          out=Rphenograph::Rphenograph(inMat[subset, ], k=p[["k"]])
          out$rows=subset
          save(out, file=outCellFile, compress="gzip")
          out
      })
    }
    save(Rout, file=outFile, compress="gzip")
    save(inMat, file=inMatFile)
  } else {
    load(outFile)
    load(inMatFile)
  }
  print(nrClusters)
  if(nrClusters == 1) {
    clusters=factor(igraph::membership(Rout[[2]]))
  } else {
    summary=sapply(names(Rout), function(cellType)  {
      print(cellType)
      out=Rout[[toString(cellType)]]
      if(length(out) >0 & is.na(out)) {
		  cat('No clustering info for ', cellType, '\n')
		  return(NULL)
	  }
      
      cluster=factor(igraph::membership(out[[2]]))
      if(cellType != "black")
        cluster=paste(cellType, cluster)
      cbind(cluster=cluster, row=out[[3]])
    })
    summary=data.frame(do.call(rbind, summary), stringsAsFactors = F)
    clusters=as.character(summary$cluster[match(1:nrow(inMat), summary$row)])
    clusters[is.na(clusters)] = colors[is.na(clusters)]
  }
  clusterSummary<-sapply(colnames(inMat), function(x) {
    tapply(t(inMat[[x]]), as.factor(clusters), mean)
  })
  clusterSummary[clusterSummary == 0]=min(clusterSummary[clusterSummary>0])
  pdf(pdfOut, useDingbats = F)
  pheatmap::pheatmap(apply(clusterSummary, 2, to_zscore), 
  	symm=F, scale = 'none', clustering_method = 'complete')
  dev.off()
  # igraph object
  # plot(Rout[[1]], cex=0.1)
  save(clusters,  file = clFile)
  return(clusters)

	
}

# FastPG
run_FastPG <- function(inMat, p, colors) {
	outDir=f("{p$run_id}_aux")
  if(! dir.exists(outDir))
		dir.create(outDir)
		
  outFile=f("{outDir}/out.RData.gz")
  inMatFile=f("{outDir}/inMat.RData")
  clFile=f("{outDir}/clusters.RData")
	
	pdfOut=paste0(p$run_id, ".heatmap.pdf")
	if(file.exists(clFile)) {
		load(clFile)
		return(clusters)
	}
	
	nrClusters=length(unique(colors[!is.na(colors)]))
	print(table(colors))
	print('Calling FastPG')
	if(! file.exists(outFile) | ! file.exists(inMatFile)) {
		if(nrClusters == 1) {
			if(nrow(inMat) - 2 < p[["k"]])
				return(NA)
			print('Calling FastPG on all data')
			print(class(inMat))
			Rout=FastPG::fastCluster(as.matrix(inMat), k=p[["k"]], num_threads=p[['num_threads']])
    } else {
        Rout=tapply(1:nrow(inMat), colors, function(subset) {
          cluster=colors[subset[1]]
          outCellFile=f("{outDir}/{cluster}.cell.out.RData.gz")
          if(file.exists(outCellFile))  {
            load(outCellFile)
            return(out)
          }
          cat("\n", cluster, length(subset), "\n")
		  
          # Minimum group size for clustering
          # EDIT: if(length(subset) <= p[["k"]] * 150 & nrow(inMat) > p[["k"]] * 500) return(NA)
		  if(length(subset) - 2 < p[["k"]]) return(NA)
		  out=FastPG::fastCluster(as.matrix(inMat[subset, ]), 
		  	k=p[["k"]], 
			num_threads=p[['num_threads']])
          out$rows=subset
          save(out, file=outCellFile, compress="gzip")
          out
      })
    }
    save(Rout, file=outFile, compress="gzip")
    save(inMat, file=inMatFile)
  } else {
    load(outFile)
    load(inMatFile)
  }
  
  if(nrClusters == 1) {
    clusters=Rout$communities #factor(igraph::membership(Rout[[2]]))
  } else {
    summary=sapply(names(Rout), function(cellType)  {
      # cat(cellType, sum(clusters == cellType), '\n')
      out=Rout[[toString(cellType)]]
	  
      if(length(out) >0 & is.na(out)) {
		  cat('No clustering info for ', cellType, '\n')
          return(NULL)
	  }
      cluster=out$communities # factor(igraph::membership(out[[2]]))
	  
      if(cellType != "black")
        cluster=paste(cellType, cluster)
      cbind(cluster=cluster, row=out[[3]])
    })
    summary=data.frame(do.call(rbind, summary), stringsAsFactors = F)
    clusters=as.character(summary$cluster[match(1:nrow(inMat), summary$row)])
    clusters[is.na(clusters)] = colors[is.na(clusters)]
  }
  clusterSummary<-sapply(colnames(inMat), function(x) {
    tapply(t(inMat[[x]]), as.factor(clusters), mean)
  })
  clusterSummary[clusterSummary == 0]=min(clusterSummary[clusterSummary > 0])
  pdf(pdfOut, useDingbats = F)
  pheatmap::pheatmap(apply(clusterSummary, 2, to_zscore), symm=F)
  dev.off()
  # igraph object
  # plot(Rout[[1]], cex=0.1)
  save(clusters,  file = clFile)
  return(clusters)
}

# Rphenograph metaclusters
run_MC <- function(inMat, p, colors) {
  
	
	outDir=f("{p$run_id}_aux")
  if(! dir.exists(outDir))
		dir.create(outDir)
		
  outFile=f("{outDir}/out.RData.gz")
  inMatFile=f("{outDir}/inMat.RData")
  clFile=f("{outDir}/clusters.RData")
	
  if(!file.exists(outFile)) {
    Rout=Rphenograph::Rphenograph(inMat, k=p[["k"]])
    save(Rout, file=outFile, compress="gzip")
    save(inMat, clusters, file=inMatFile)
  } else {
    load(outFile)
    load(inMatFile)
  }
  pheatmap::pheatmap(to_magnitude(inMat, p$magnitude))
  clusters=factor(membership(Rout[[2]]))
  # plot(Rout[[1]], cex=0.1)
  save(clusters, file = clFile)
  return(clusters)
}

# flowSOM
run_flowSOM <- function(inMat, p, colors) {
  
  # p         - is a list of all parameters for the functions: 
  # Functions - ReadInput, BuildSOM, BuildMST, and Metaclutersing
  #             all embedded in the wrapper function FlowSOM
  # colsToUse - the columns to use will be preselected in the preselectedious step
  # nClus     - Metaclustering option
  # seed      - 42 (Seed for reproducible results)
  # metaclustering method options: metaClustering_consensus,
  #         metaClustering_hclust,metaClustering_kmeans,metaClustering_som
  
	
	outDir=f("{p$run_id}_aux")
  if(! dir.exists(outDir))
		dir.create(outDir)
		
  outFile=f("{outDir}/out.RData")
  inMatFile=f("{outDir}/inMat.RData")
  clFile=f("{outDir}/clusters.RData")
	
  nrClusters=length(unique(colors[!is.na(colors)]))
  print(table(colors))
  if(!file.exists(outFile) | !file.exists(inMatFile)) {
    if(nrClusters == 1) {
      obj=flowCore::flowFrame(as.matrix(inMat))
	  print('ReadInput')
      input=FlowSOM::ReadInput(obj, # transform=as.logical(p[["transform"]]),
                               # compensate=as.logical(p[["compensate"]]),
                               scale=as.logical(p[["scale"]]))
	  print('FlowSOM')
      fSOM=FlowSOM::BuildSOM(input, xdim=p[["xdim"]],  ydim=p[["ydim"]], rlen=p[["rlen"]])
	  print('Build')
      fMST=FlowSOM::BuildMST(fSOM, tSNE=p[["tSNE"]])
    } else {
      fMST=tapply(1:nrow(inMat), colors, function(subset) {
        cluster=colors[subset[1]]
        outCellFile=f("{outDir}/{cluster}.cell.out.RData.gz")
        print(outCellFile)
        if(file.exists(outCellFile)) {
          load(outCellFile)
          return(out)
        }
        cat("\n", cluster, length(subset), "\n")
        if(length(subset) <= p[["xdim"]] * p[["ydim"]]) return(NA)
        obj=flowCore::flowFrame(as.matrix(inMat[subset, ]))
        input=FlowSOM::ReadInput(obj, scale=as.logical(p[["scale"]]))
        fSOM=FlowSOM::BuildSOM(input, xdim=p[["xdim"]],  ydim=p[["ydim"]], rlen=p[["rlen"]])
        out=FlowSOM::BuildMST(fSOM, tSNE=p[["tSNE"]])
        out$rows=subset
        save(out, file=outCellFile, compress="gzip")
        out
      })
    }
    save(fMST, file=paste0(p$run_id, ".out.RData"))
    save(inMat, file = inMatFile)
  } else {
    load(outFile)
    load(inMatFile)
  } 

  print('Exporting clusters')
  print(nrClusters)
  if(nrClusters==1) {

	print('Meta')
    runOutput<-FlowSOM::MetaClustering(fMST$map$codes,  max=p[["maxClust"]],
            method="metaClustering_consensus")
	print('Mapping')
    clusters=runOutput[fMST$map$mapping[,1]]
	print('heatmap')
    pheatmap::pheatmap(apply(fMST$map$codes, 2, to_percentile))
    FlowSOM::PlotStars(fMST) #, backgroundValues=as.factor(clusters))

  } else {

    summary=sapply(names(fMST), function(cellType)  {
      print(cellType)
      out=fMST[[toString(cellType)]]
      if(length(out) >  0 & is.na(out))
        return(NULL)
	  print(cellType)
      runOutput<-FlowSOM::MetaClustering(out$map$codes, max=p[["maxClust"]],
                                         method="metaClustering_consensus")
	  print('clustered')
      cluster=runOutput[out$map$mapping[,1]]
      annotation_row=cbind(cluster=runOutput)
      # pheatmap::pheatmap(apply(out$map$codes, 2, to_percentile), main = cellType)
                         # annotation_row = cluster)
	  print('pheatmap')
      pheatmap::pheatmap(apply(out$map$codes, 2, to_min_max), main = cellType,
                         annotation_row = annotation_row)
     
      if(cellType != "black")
        cluster=paste(cellType, cluster)
      cbind(cluster=cluster, row=out$rows)
    })
    summary=data.frame(do.call(rbind, summary), stringsAsFactors = F)
    clusters=as.character(summary$cluster[match(1:nrow(inMat), summary$row)])
    clusters[is.na(clusters)] = "Excluded"
  }
  save(clusters, file=clFile)
  return(clusters)
}

# rTSNE
run_rtsne <- function(inMat, p, colors)  {
  
  # check_duplicates: cells with the same value after PCA reduction
  # (best to make sure there are no duplicates, esp. for large datasets)
  # the colors/clusters can be provided under pars$rtsne$clusters
  
	outDir=f("{p$run_id}_aux")
  if(! dir.exists(outDir))
		dir.create(outDir)
		
  outFile=f("{outDir}/out.RData.gz")
  inMatFile=f("{outDir}/inMat.RData")
  probsFile=f("{outDir}/probs.txt")
  clFile=f("{outDir}/clusters.RData")
	
  print(p[["nthread"]])
  if(!file.exists(outFile)) {
    rtsne_out<-Rtsne::Rtsne(as.matrix(inMat),
                       pca=as.logical(p[["pca"]]),
                       verbose=as.logical(p[["verbose"]]),
                       check_duplicates= p[["check_duplicates"]],
                       num_threads=p[["nthread"]])
    save(rtsne_out, p, file=outFile)
    save(inMat, p, file=inMatFile)
  } else {
    load(outFile)
    load(inMatFile)
  }
  clusterColors=sapply(colors, function(x) {
    if(! x %in% names(cellTypeColors)) return('grey')
    cellTypeColors[[x]]
  })
  tsnePdf=paste0(p$run_id, ".tsne.pdf")
  pdf(tsnePdf, width=7, height=7)
  plot(rtsne_out$Y, pch=20, col=clusterColors, xlab="", ylab="", 
       xaxt='n', yaxt='n', cex=0.5)
  plot=plot_scatter(x=rtsne_out$Y[,1], y=rtsne_out$Y[,2],
                   color=clusterColors, xlab="tSNE1", ylab="tSNE2")
  print(plot)
  dev.off()
  # points(rtsne_out$Y, col="darkgrey", pch=1)
  return(colors)
}

# UMAP
run_umap <- function(inMat, p, colors)  {
  outFile=paste0(p[["run_id"]], ".umap.RData")
  if(file.exists(outFile)) {
    load(outFile)
  } else {
    out=umap::umap(data.table::setDF(inMat), random_state=555)
    save(out, colors, p, file = outFile)
  }
  # plot only subsampled from each cluster
  umapPdf=paste0(p$run_id, ".umap.pdf")
  pdf(file=umapPdf, width=7, height = 7)
  names(clusterColors)=sort(unique(colors))
  plotUmap(x=out, labels=colors, colors=palette$cellTypeColors, legend=T, cex=0.1)
  dev.off()
} 

# kmeans
run_kmeans <- function(inMat, p, colors) {
  out=stats::kmeans(inMat, centers=p[["centers"]])
  pheatmap(apply(out$centers, 2, to_percentile))
  return(out$cluster)
}

# RClusterpp - not tested
run_rclusterpp <- function(inMat, p, colors) {
  Rclusterpp.hclust(inMat, method=p[["method"]])
}

# immunoClust - not tested
run_immunoClust <- function(inMat, p, colors) {
  obj<-convert_mat2fcs(inMat)
  out<-immunoClust::cell.process(obj, classify.all=p[["classify.all"]])
  return(out@label)
}
