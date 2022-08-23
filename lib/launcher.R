
# Generic function
run_method <- function(inData, method, pars, runID, wDir, regFile, nfDir,
                     feature="MeanIntensity", colors=NULL)  {
  
  # inData     - a data frame object with cells in rows and markers in columns
  # method     - accepts the value of Rphenograph, X-shift, CellAssign, kmeans FlowSOM
  # feature    - the value that summarizes the intensity signal (mean, integrated, median)
  # parameters - method-specific parameters + contains markers to be excluded/included
  
  # If file.pattern=Cells.csv
  # columnNames=get_column_names(colnames=colnames(inData),
  #                                features=feature, markers_select=group,
  #                                channels_exclude=pars$channels_exclude)
  
  ids=with(inData, paste(ObjectNumber, basename(file)))
  colnames(inData)=markers_format(colnames(inData))
  funForward=match.fun(paste0("run_", method))
  pars[[method]]$run_id=runID
  pars[[method]]$markers=pars$markers
  pars[[method]]$magnitude=pars$magnitude
  
  # Load info
  areaData=load_files(file.path(nfDir, "AreaShape"), run=pars$run)
  areaMatch=match(ids, with(areaData, paste(ObjectNumber, basename(file))))
  
  # if the analysis has been done
  inMatFile=f("{runID}.inMat.RData")
  clFile=f("{runID}.clusters.RData")
  
  if(! file.exists(inMatFile))  {
    # Keep only numeric columns
    columnNames=setdiff(colnames(inData), c(pars$channels_exclude, "ObjectNumber", "file"))
    rowsKeep=rep(TRUE, nrow(inData))
  
    if(!is.null(pars$markers) & pars$markers != "all") {
      indices=match(unique(unlist(marker_gene_list[[pars$markers]])), colnames(inData))
      keep=colnames(inData)[indices[!is.na(indices)]]
      cat("Column names matching the markers list:", keep, "\n", file=f("{runID}.log"), append=T)
      columnNames=intersect(columnNames, keep)
    }
    colnames(inData)=gsub(paste0(".*", feature, "_([^_]+).*"), "\\1", colnames(inData))
    
    # Exlude cells with an area smaller or equal to pars$area_exclude
    if(!is.null(pars$area_exclude) & pars$area_exclude > 0) {
      cat("Number of cells with area smaller than", pars$area_exclude, "is",
          sum(areaData$AreaShape_Area[areaMatch] <= pars$area_exclude), "\n",
          file=f("{runID}.log"), append=T)
      rowsKeep=rowsKeep & areaData$AreaShape_Area[areaMatch] > pars$area_exclude
    }
    
    # Exclude cells with intensity smaller or equal to pars$min_expression for all markers
    rowsKeep=apply(inData[, ..columnNames], 1, function(x) any(x > pars$min_expresssion))
    cat("Number of cells with zero expression for all markers of interest",
        sum(!rowsKeep), "\n", file=f("{runID}.log"), append=T)
    
    # Include batches for cellassign; # rows must match cov matrix
    if(method == "cellassign") {
      # add batch effects
      tma_block=get_region_info(panel=pars$panel, cellIDs=inData$file, featurePattern="TMA$", regFile=regFile)
      tma_side=get_region_info(panel=pars$panel, cellIDs=inData$file, featurePattern="duplicate_core$", regFile=regFile)
      cohort=get_region_info(panel=pars$panel, cellIDs=inData$file, featurePattern="^cohort$", regFile=regFile)
      tma_values=paste(cohort, tma_block, tma_side)
      if(any(is.na(tma_values))) {
        stop("Missing TMA values")
      }
      cat("Number of cells without TMA & tissue values is",
          sum(is.na(tma_values)), "\n",
          file=file.path(paste0(runID, ".log")), append=T)
      rowsKeep=rowsKeep & !is.na(tma_values)  # & !is.na(tissue_values)
      print(table(tma_values[rowsKeep]))
      if(length(unique(tma_values[rowsKeep])) > 1) {
        pars[[method]]$X=model.matrix(~ 0 + tma_values[rowsKeep])   
      } else {
        pars[[method]]$X=NULL
      }
      # create_covariate_matrix(list(tma_values[rowsKeep], tissue_values[rowsKeep]))
    }
    
    if(! dir.exists(dirname(runID)))
      dir.create(dirname(runID), recursive=T)
     
    inMat=inData[rowsKeep, ..columnNames]
	print(dim(inMat))
    colors=colors[rowsKeep] # toString
    if(pars[[method]]$transformation != "none") {
      transformInput=match.fun(paste0("to_", pars[[method]]$transformation))
      inMat=transformInput(inMat)
    }
    
    runOutput=rep("Excluded", nrow(inData))
    rownames(inMat) = ids[rowsKeep]
    cat("Markers for analysis:\n", colnames(inMat), "\n", 
        file=f("{runID}.log"), append=T)
    if(method == "MC")  {
      print("Calculating centroid for each marker")
      inMat=apply(inMat, 2, function(x) tapply(x, as.factor(colors), median))
      metaclusters=funForward(inMat, pars[[method]], colors=colors)
      runOutput[rowsKeep]=metaclusters[match(colors, rownames(inMat))]
    } else {
      runOutput[rowsKeep]=funForward(inMat, pars[[method]], colors=colors)
    }
  } else {
    load(inMatFile)
    rowsKeep=match(rownames(inMat), ids)
	print(head(rownames(inMat)))
	print(sum(is.na(rowsKeep)))
	print(head(ids[is.na(rowsKeep)]))
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
  if(method %in% pars$visualisation_methods) return(NULL)
  return(list('runOutput'=runOutput, 'columnNames'=columnNames, 'ids'=ids))
}

run_Rphenograph <- function(inMat, p, colors) {

  library(Rcpp)	
  outFile=paste0(p$run_id, ".out.RData.gz")
  inMatFile=paste0(p$run_id, ".inMat.RData")
  clFile=paste0(p$run_id, ".clusters.RData")
  pdfOut=paste0(p$run_id, ".heatmap.pdf")
  if(file.exists(clFile)) {
    load(clFile)
    return(clusters)
  }
  nrClusters=length(unique(colors[!is.na(colors)]))
  print(table(colors))
  if(!file.exists(outFile) | !file.exists(inMatFile)) {
    if(nrClusters == 1) {
      Rout=Rphenograph::Rphenograph(inMat, k=p[["k"]])
    } else {
        Rout=tapply(1:nrow(inMat), colors, function(subset) {
          cluster=colors[subset[1]]
          outCellFile=paste0(p$run_id, ".", cluster, ".cell.out.RData.gz")
          if(file.exists(outCellFile))  {
            load(outCellFile)
            return(out)
          }
          cat("\n", cluster, length(subset), "\n")
          # Minimum group size for clustering 
          if(length(subset) <= p[["k"]] * 150) return(NA)
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
  if(nrClusters == 1) {
    clusters=factor(igraph::membership(Rout[[2]]))
  } else {
    summary=sapply(names(Rout), function(cellType)  {
      print(cellType)
      out=Rout[[toString(cellType)]]
      if(length(out) >0 & is.na(out))
        return(NULL)
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
  pheatmap::pheatmap(apply(clusterSummary, 2, to_zscore), symm=F)
  dev.off()
  # igraph object
  # plot(Rout[[1]], cex=0.1)
  save(clusters,  file = clFile)
  return(clusters)

	
}
# Rphenograph
run_FastPG<-function(inMat, p, colors) {
  
  outFile=paste0(p$run_id, ".out.RData.gz")
  inMatFile=paste0(p$run_id, ".inMat.RData")
  clFile=paste0(p$run_id, ".clusters.RData")
  pdfOut=paste0(p$run_id, ".heatmap.pdf")
  if(file.exists(clFile)) {
    load(clFile)
    return(clusters)
  }
  nrClusters=length(unique(colors[!is.na(colors)]))
  print(table(colors))
  print('Calling FastPG')
  if(!file.exists(outFile) | !file.exists(inMatFile)) {
    if(nrClusters == 1) {
	  print('Calling FastPG on all data')
	  print(class(inMat))
	  print(dim(inMat))
      Rout=FastPG::fastCluster(as.matrix(inMat), k=p[["k"]], num_threads=p[['num_threads']])
    } else {
	   print(dim(inMat))
        Rout=tapply(1:nrow(inMat), colors, function(subset) {
          cluster=colors[subset[1]]
          outCellFile=paste0(p$run_id, ".", cluster, ".cell.out.RData.gz")
          if(file.exists(outCellFile))  {
            load(outCellFile)
            return(out)
          }
          cat("\n", cluster, length(subset), "\n")
          # Minimum group size for clustering 
          if(length(subset) <= p[["k"]] * 150) return(NA)
          out=FastPG::fastCluster(as.matrix(inMat[subset, ]), k=p[["k"]], num_threads=p[['num_threads']])
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
	print(head(clusters))
  } else {
    summary=sapply(names(Rout), function(cellType)  {
      print(cellType)
      out=Rout[[toString(cellType)]]
      if(length(out) >0 & is.na(out))
        return(NULL)
      cluster=out$communities #factor(igraph::membership(out[[2]]))
      if(cellType != "black")
        cluster=paste(cellType, cluster)
      cbind(cluster=cluster, row=out[[3]])
    })
    summary=data.frame(do.call(rbind, summary), stringsAsFactors = F)
    clusters=as.character(summary$cluster[match(1:nrow(inMat), summary$row)])
	cat('Number of clusters with na:', sum(is.na(clusters)), '\n')
    clusters[is.na(clusters)] = colors[is.na(clusters)]
  }
  
  clusterSummary<-sapply(colnames(inMat), function(x) {
    tapply(t(inMat[[x]]), as.factor(clusters), mean)
  })
  clusterSummary[clusterSummary == 0]=min(clusterSummary[clusterSummary>0])
  pdf(pdfOut, useDingbats = F)
  pheatmap::pheatmap(apply(clusterSummary, 2, to_zscore), symm=F)
  dev.off()
  # igraph object
  # plot(Rout[[1]], cex=0.1)
  save(clusters,  file = clFile)
  return(clusters)
}

# Rphenograph metaclusters
run_MC<-function(inMat, p, colors) {
  
  outFile=paste0(p$run_id, ".out.RData.gz")
  inMatFile=paste0(p$run_id, ".inMat.RData")
  clFile=paste0(p$run_id, ".clusters.RData")
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

# CellAssign
run_cellassign<-function(inMat, p, colors)  {
  
  # Note: marker_gene_info requires an object of class matrix as input
  # marker_mat rows must match inmat columns
  
  outFile=paste0(p$run_id, ".out.RData")
  inMatFile=paste0(p$run_id, ".inMat.RData")
  probsFile=paste0(p$run_id, ".probs.txt")
  clFile=paste0(p$run_id, ".clusters.RData")
  if(!file.exists(outFile)) {
    markerMat=cellassign::marker_list_to_mat(marker_gene_list[[p$markers]])
    markerMat=markerMat[rownames(markerMat) %in% colnames(inMat), ]
    markerMat=markerMat[, apply(markerMat, 2, function(x) !all(x == 0))]
    # if same signature for multiple cells because of missing markers, keep the higest in lineage
    cellcode=apply(markerMat, 2, paste, collapse="")
    exclude=unlist(sapply(unique(cellcode[duplicated(cellcode)]), function(x) {
      duplcates=names(cellcode)[cellcode==x]
      count=sapply(duplcates, function(celltype) {
        length(marker_gene_list[[p$markers]][[celltype]])
      })
      names(count)[count > min(count)]
    }))
    if(length(exclude) > 0) {
      markerMat=markerMat[, ! colnames(markerMat) %in% exclude]
    }
    write.table(markerMat, append = T, file=f("{p$run_id}.log"))
    print(markerMat)
    cat("Using", nrow(markerMat), "markers and ", ncol(markerMat),  "celltypes.\n")
    keep=colnames(inMat)[match(rownames(markerMat), colnames(inMat))]
    cat("From them,", length(keep), "are found in inMat.\n")
    
    s=rep(1, nrow(inMat))
    fit<-cellassign::cellassign(exprs_obj=as.matrix(inMat[, ..keep]),
                    marker_gene_info=markerMat,
                    s=s, X=p$X,
                    learning_rate=p$learning_rate,
                    shrinkage=p$shrinkage,
                    verbose=p$verbose)
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

# flowSOM
run_flowSOM<-function(inMat, p, colors) {
  
  # p         - is a list of all parameters for the functions: 
  # Functions - ReadInput, BuildSOM, BuildMST, and Metaclutersing
  #             all embedded in the wrapper function FlowSOM
  # colsToUse - the columns to use will be preselected in the previous step
  # nClus     - Metaclustering option
  # seed      - 42 (Seed for reproducible results)
  # metaclustering method options: metaClustering_consensus,
  #         metaClustering_hclust,metaClustering_kmeans,metaClustering_som
  
  outFile=paste0(p$run_id, ".out.RData")
  inMatFile=paste0(p$run_id, ".inMat.RData")
  clFile=paste0(p$run_id, ".clusters.RData")
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
        outCellFile=paste0(p$run_id, ".", cluster, ".cell.out.RData.gz")
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
run_rtsne<-function(inMat, p, colors)  {
  
  # check_duplicates: cells with the same value after PCA reduction
  # (best to make sure there are no duplicates, esp. for large datasets)
  # the colors/clusters can be provided under pars$rtsne$clusters
  
  print(p[["nthread"]])
  outFile=paste0(p[["run_id"]], ".out.RData")
  inMatFile=paste0(p[["run_id"]], ".inMat.RData")
  if(!file.exists(outFile)) {
    rtsne_out<-Rtsne(as.matrix(inMat),
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

# immunoClust
run_immunoClust<-function(inMat, p, colors) {
  obj<-convert_mat2fcs(inMat)
  out<-immunoClust::cell.process(obj, classify.all=p[["classify.all"]])
  return(out@label)
}

# kmeans
run_kmeans<-function(inMat, p, colors) {
  out=stats::kmeans(inMat, centers=p[["centers"]])
  pheatmap(apply(out$centers, 2, to_percentile))
  return(out$cluster)
}

# RClusterpp
run_rclusterpp<-function(inMat, p, colors) {
  Rclusterpp.hclust(inMat, method=p[["method"]])
}

run_umap<-function(inMat, p, colors)  {
  outFile=paste0(p[["run_id"]], ".umap.RData")
  if(file.exists(outFile)) {
    load(outFile)
  } else {
    out=umap::umap(data.table::setDF(inMat), random_state=555)
    save(out, colors, p, file = outFile)
  }
  umapPdf=paste0(p$run_id, ".umap.pdf")
  pdf(file=umapPdf, width=7, height = 7)
  names(clusterColors)=sort(unique(colors))
  plotUmap(x=out, labels=colors, colors=clusterColors, legend=F, cex=0.2)
  dev.off()
}

get_probabilities <- function(cellIDs, clusters, probFile) {
  
  if(!file.exists(probFile))
     return(rep(NA, nrow(inData)))
  probs=data.table::fread(probFile, sep = "\t")
  probMatch=match(cellIDs, probs$imagename)
  probs=probs[probMatch, ]
  export=tapply(1:length(cellIDs), clusters, function(subset) {
    cluster=clusters[subset[1]]
    print(cluster)
    if(! cluster %in% colnames(probs)) {
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
  names=clusterNames$positive[match(clusters, clusterNames$cluster)]
  names=gsub("up:(.*) down:.*", "\\1", names)
  names=gsub("pos:(.*) neg:.*", "\\1", names)
  names[is.na(names)] = ''
  getExpSummary=match.fun(fun)
  summary=tapply(1:nrow(dfExp), names, function(subset) {
    cluster=names[subset[1]]
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
