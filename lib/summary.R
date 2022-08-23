# Generic function

summarise_output <- function(inData, method, pars, runID, runOutput, columnNames, nfDir, colors=NULL) {
  
		ids=with(inData, paste(ObjectNumber, basename(imagename)))

	    colnames(inData)=colnames(inData) %>%
	  		gsub(paste0(".*", feature, "_([^_]+).*"), "\\1", .)
			
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

		clusterNames=assign_cluster_positivity(
			dfExp = inData[, ..columnNames], 
			clusters=runOutput,
			run=pars$run,
			runID=runID, 
			panel=pars$panel, 
			magnitude=pars$magnitude,
			fromList=ifelse(method == "cellassign", pars$markers, NA)
		)

		majorTypes=NULL
		if(method != 'cellassign' & ! all(is.null(colors))) {
			if(!all(colors=='black')) 
				majorTypes=gsub(" [0-9]+", "", colors)
		}
		clusterNames$positive=sapply(clusterNames$positive, toString)
		positive=clusterNames$positive[match(runOutput, clusterNames$cluster)]
		positive[is.na(positive)] = 'pos: neg:'

		print("Getting tissue category")
		tissue_categs=get_tissue_category(
			cellIDs = gsub(".txt", "", ids),
			panel = pars$panel,
			tissAreaDir=pars$tissAreaDir,
			class="regional"
		)
		vessel_categs=get_tissue_category(
			cellIDs = gsub(".txt", "", ids),
			panel = pars$panel,
			tissAreaDir=pars$tissAreaDir,
			class="vessels"
		)
		colnames(vessel_categs) = gsub('region', 'vessels_region', colnames(vessel_categs))
			
		print('Major type assignment')
		majorTypes=assign_celltype(
			names = positive, 
			majorTypes=majorTypes,
			markers = marker_gene_list[[pars$ref_markers]],
			major = T, 
			region=tissue_categs$Tumour
		)

		print("Getting marker expression")
		meanFeature=get_markers_expression(
			inData[, ..columnNames],
			clusters=runOutput,
            clusterNames = clusterNames, 
			fun="mean", 
            magnitude = pars$magnitude
		)

		print("Getting spatial coordinates")
		spatData=load_files(file.path(nfDir, "LocationCenter"), run=pars$run)
		spatMatch=match(with(inData, paste(ObjectNumber, basename(imagename))),
				  with(spatData, paste(ObjectNumber, basename(imagename))))

		print('Getting area')
		areaData=load_files(file.path(nfDir, "AreaShape"), run=pars$run)
		areaMatch=match(ids, with(areaData, paste(ObjectNumber, basename(imagename))))

		print("Getting probabilities")
		probabilities=get_probabilities(cellIDs = ids, clusters=runOutput,
								  	probFile = f("{runID}.probs.txt"))

		if(method == "cellassign")	{
			celltypes=runOutput
			celltypes=gsub('CD([48])$', 'CD\\1 T cells', celltypes)
			clusterNames$cluster=gsub('CD8$', 'CD8 T cells', clusterNames$cluster)
			clusterNames$major=clusterNames$cluster
		} else {
			print('Cell subtype assignment')
			# from nested list to list
			celltypes=assign_celltype(
				names = positive, 
				markers = markersList, 
				majorTypes = majorTypes,
				region=tissue_categs$Tumour
			)
			if(pars$subset == subtypesDir) {
				clusterNames$major = 
					gsub("[0-9]+$|Excluded$", "", clusterNames$cluster) %>%
					gsub(' $', '', .)
			} else {
				clusterNames$major=NULL
			}
		}

		clusterNames$majorType=assign_celltype(
			names = clusterNames$positive,
			majorTypes=clusterNames$major,
			markers = marker_gene_list[[pars$ref_markers]],
			major = T
		)
		clusterNames$cellType=assign_celltype(
			names = clusterNames$cellType, 
			markers = markersList,
			majorTypes = clusterNames$majorType
		)
			
		write.tab(clusterNames, file=clusterNameFile)
		nrClusters=length(unique(runOutput))
        print('Marker expression heatmap')
		heatmapOrder=plot_heatmap(dfExp = inData[, ..columnNames],
		                                clusters = runOutput, runID = runID,
		                                labels=clusterNames)


		cat('Marker expression heatmap for ', nrClusters, ' clusters\n')

		print('Reviewing cell assignments for plotting only')
		print(pars$cellAssignFile)
		clusterNames$cellType_rev=review_cellType_by_major(
			cellTypes      = clusterNames$cellType,
			majorTypes	   = clusterNames$majorType,
			positivity	   = clusterNames$positive, 
			cellAssignFile = pars$cellAssignFile
		)
		clusterNames$majorType_rev=review_major_by_cellType(
			cellTypes	   = clusterNames$cellType,
			majorTypes	   = clusterNames$majorType,
			positivity	   = clusterNames$positive,
			cellAssignFile = pars$cellAssignFile
		)

	  print("Plotting expression")
	  plot_expression(
		  dfExp = inData[, ..columnNames],
          clusters = runOutput,
          clusterNames, pars[[method]],
          magnitude = pars$magnitude
	  )

		print("Wrapping up in data frame")
		outDF=data.frame(
			imagename     = gsub(".txt", "", basename(inData$imagename)),
			object        = sapply(inData$ObjectNumber, toString),
			centerX       = spatData$LocationCenter_X[spatMatch],
			centerY       = spatData$LocationCenter_Y[spatMatch],
			runID         = gsub(paste0(pars$inDir, "/?"), "", runID),
			cluster       = sapply(runOutput, toString),
			cellType      = celltypes,
			majorType     = majorTypes,
			positive      = positive,
			meanIntensity = meanFeature,
			area          = areaData$AreaShape_Area[areaMatch],
			probability   = probabilities,
			vessel_categs,
			tissue_categs
		)
		write.tab(outDF, file=f("{runID}.clusters.txt"))
		fst::write_fst(outDF, path=f("{runID}.clusters.fst"), compress = 75)

		return(outDF)
}
  
get_probabilities <- function(cellIDs, clusters, probFile) {
  
		if(!file.exists(probFile))
		return(rep(NA, nrow(inData)))
		probs=data.table::fread(probFile, sep = "\t")
		probMatch=match(cellIDs, probs$imagename)
		probs=probs[probMatch, ]
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

		if(! is.null(magnitude)) 
			dfExp=to_magnitude(dfExp, magnitude)

		clusterNames$cluster=gsub('\\.', ' ', clusterNames$cluster)
		names=clusterNames$positive[match(clusters, clusterNames$cluster)]
		names=gsub("up:(.*) down:.*", "\\1", names)
		names=gsub("pos:(.*) neg:.*", "\\1", names)
		names[is.na(names)] = ''
		getExpSummary=match.fun(fun)
		summary=tapply(1:nrow(dfExp), names, function(subset) {
				cluster=names[subset[1]]
				# cat('Expression values for', cluster, '\n')
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
