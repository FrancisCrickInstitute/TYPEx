get_subtree <- function(tree, value) {
	
	if(! is.list(tree))
		return(NULL)
	if(value %in% names(tree))
		return(tree[[value]])
	indices = lapply(names(tree), function(label)
		get_subtree(tree[[label]], value)
	)
	if(is.null(unlist(indices))) return(NULL)
	return(rlist::list.clean(indices, fun=is.null, recursive=T)[[1]])
}

get_node_depth <- function(tree, value, level = 1) {
	
	if(length(tree) == 1) {
		if(tree == value) 
			return(level)
		return(NULL)
	}
	if(value %in% names(tree)) {
		names(level) = value
		return(level)
	}
	indices = sapply(names(tree), function(label)
		get_node_depth(tree[[label]], value, level + 1)
	)
	indices=unlist(indices)
	return(indices[! is.na(indices)])
}

isChild <- function(tree, parent, child) {
	
	if(parent == child)
		return(TRUE)
	
	parentSubtree=get_subtree(tree, parent)
	if(! is.list(parentSubtree) | is.null(parentSubtree))
		return(FALSE)
	childIndex=get_node_depth(parentSubtree, child)
	return(! is.null(childIndex))
}

get_node_breadth <- function(tree, value, tree_path) {
	
	levels=strsplit(tree_path, split = '\\.')[[1]]
	values=tree[levels[1]]
	if(value %in% names(values))
		return(which(names(values) == value))
	for(i in 1:(length(levels)+1)) {
		values=values[[levels[i]]]
		if(value %in% names(values))
			return(which(names(values) == value))
	}
	return(NA)
}

remove_node <- function(tree, value) {
	lapply(tree, function(x)  {
		if( is.list(x) ) { 
			if(! is.null(names(x))) {
				remove_node(x[names(x) != value], value)
			} else {
				remove_node(x, value)
			}
		} else x
	})
}

get_celltype_markers <- function(tree, celltype) {
	if(! is.list(tree))
		return(NULL)
	indices = sapply(names(tree), function(label) {
		if(label == celltype | celltype %in% gsub(' - Other', '', label)) {
			if(is.list(tree[[label]])) {
				get_celltype_markers(tree[[label]], f("{celltype} - Other")) %>%
					return
			}
			return(tree[[label]])
		}
		get_celltype_markers(tree[[label]], celltype)
	})
	indices=indices %>% unlist %>% unname
	return(indices[! is.na(indices)])
}
  
get_celltypes <- function(tree, label = NULL) {
	
	# return only the cell types for which major_markers have been defined
	if(! is.list(tree) & is.null(label))
		return(NULL)
	
	if(! is.list(tree))
		return(toString(label))
	cellTypes = sapply(names(tree), function(label) {
		get_celltypes(tree[[label]], label)
	})
	cellTypes=cellTypes %>% unlist %>% unname
	return(cellTypes[! is.na(cellTypes)])
}

# convert nested to a list of the final nodes
get_node_list <- function(tree) {

	if(! is.list(tree))
		return(tree)
	cellTypes = get_celltypes(tree)
	sapply(cellTypes, 
		function(x) get_celltype_markers(tree, x))
}
