
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
		
		#overlap=sapply(cellTypeNames, function(celltype) {
		#	sum(cellMarkers %in% get_celltype_markers(markers, celltype))
		#})
		
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