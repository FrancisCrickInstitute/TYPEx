# NUMERICAL TRANSFORMATIONS

# TRANSFORMATIONS
to_magnitude <- function(inMat, dynamic_range=10 ** 6) {
  dynamic_range * inMat
}

to_center <- function(inMat) {
  apply(inMat, 2, function(x) scale(x, scale =F))
}

to_max <- function(inMat) {
  apply(inMat, 2, function(x) x/max(x, na.rm=TRUE))
}

to_min_max <- function(values) {
  values[is.nan(values)]=NA
  if(max(values, na.rm=TRUE) == 0) return(values)
  (values - min(values, na.rm=T))/(max(values, na.rm=TRUE) - min(values, na.rm=TRUE))
}

to_mean_center <- function(values) {
  values[is.nan(values)]=NA
  values - mean(values, na.rm = T)
}

to_zscore <- function(values) {
  (values - mean(values, na.rm = T))/sd(values, na.rm = T)
}

to_robust_zscore <- function(values) {
  (values - median(values, na.rm = T))/mad(values, na.rm = T)
}
to_percentile <- function(values) {
  fun=ecdf(values)
  sapply(values, fun)
}

to_asin <- function(values) {
  asin(values)
}

to_asinh <- function(values, magnitude = 10 ** 6) {
  if(all(values < 1))
		values = values * magnitude
  asinh(values)
}

# Shortcuts
f<-glue::glue 

paste1 <- function(...) {
  paste(..., sep="_", collapse="_")
}

list.rec <- function(...) {
  list.dirs(..., full.names = T, recursive = T)
}

list.frec <- function(...) {
  list.files(..., full.names=T, recursive = T)  
}

write.tab <- function(...) {
  write.table(..., sep = "\t", row.names = F, quote = F)
}

grepv <- function(...) {
  grep(..., value = T)
}

get.pvalue_label = function(pvalues) {
  if ( all(is.na(pvalues)) ) return( rep( "", length(pvalues) ) )
  if(!is.numeric(pvalues) ) stop( "pvalues are not numeric")
  labels = rep( "", length( pvalues) )
  labels[ pvalues <= 0.1 ] = "."
  labels[ pvalues <= 0.05 ] = "*"
  labels[ pvalues <= 0.01 ] = "**"
  labels[ pvalues <= 0.001 ] = "***"
  return(labels)
}

get_combos <- function(array) {
  if(length(array) == 1) return(array)
  combos =c(sapply(1:length(array), function(inx) {
      get_combos(array[-inx])
  }), paste(array, collapse=":"))
  return(unique(combos))
}

# MODEL MATRIX / COVARIATE
create_covariate_matrix <- function(covList) {
  X=do.call(cbind, lapply(covList, function(values) {
    do.call(cbind, lapply(names(table(values)), function(x) {
      sapply(values == x, sum)
    }))
  }))
  return(as.matrix(X))
}

# FILE UTILS
convert_fst2csv <- function(fst_file) {
  data=fst::read.fst(fst_file)
  data.table::fwrite(data, file  = gsub("fst", "csv", fst_file))
}

load_files <- function(inDir, pattern="*.txt", 
                       resample=FALSE, resample_frac=0.67,
                       full.names=T, recursive=T,
                       col.names=NULL, run=NULL,
                       col.exclude=NULL,
                       resampleBy="file")  {
  if(!dir.exists(inDir)) {
    stop("Dir does not exists: ", inDir,
         ", file pattern: ", pattern, " (sshfscamp)")
  }
  
  inFiles=list.files(path      =inDir,
                     pattern   =pattern,
                     full.names=full.names,
                     recursive =recursive)

  if(!is.null(run))
    inFiles=grep(paste0(run, "\\/"), inFiles, value=T)
  
  if(length(inFiles) == 0)
    stop("No input files at the given in ", inDir, " run:", run)
  fullSet=vector(mode="list")
  start=Sys.time()
  cat("Loading files", length(inFiles), "\n")
  for(inFile in inFiles) { # sample(1:length(inFiles), function(x))
    cat(".")
    extension=gsub('.*\\.([^.]+)$', '\\1', basename(inFile))
    if(extension %in% c('csv') | pattern == "*.txt") {
      data=data.table::fread(inFile, sep=",", header=T, check.names=F)
    } else if(extension %in% c('fst')) {
      data=read.fst(inFile)
      data=data.table::setDT(data)
    } else  if(extension=='txt') {
      data=data.table::fread(inFile, sep="\t", header=T, check.names=F)
    }
    if(!is.null(col.names))
      data=data[, colnames(data) %in% col.names, with = F]
    if(!is.null(col.exclude))
      data=data[, ! colnames(data) %in% col.exclude, with = F]
    data$file=inFile
    if(resample) {
      if(!is.null(resampleBy) && !is.na(resampleBy)) {
        data=data[unlist(tapply(1:nrow(data), data[[resampleBy]],
                                function(x) sample(x, replace=F, size=resample_frac*length(x)))), ]
      } else {
        data=data[sample(1:nrow(data), replace=F, size=resample_frac*nrow(data)), ]
      }
    }
    fullSet[[inFile]]=data
  }
  cat("\n")
  end=Sys.time()
  cat("Total number of loaded files", length(fullSet), "\n",
      end - start, "loading time", "\n", append=T)
  inData=data.table::rbindlist(fullSet, fill = T)
  cat("Total number of loaded rows",
      nrow(inData), "\n", append=T)
  cat("Total number of selected rows",
      nrow(inData), "\n", append=T)
  rownames(inData)=gsub(inDir, "", rownames(inData))
  if(any(! col.names %in% colnames(inData)))  {
    print("Warning:", col.names[! col.names %in% colnames(inData)], 
          "not found in data.\n")
  }
  colnames(inData)=markers_format(colnames(inData))
  return(inData)
}


plotUmap <- function(x, labels, main="",
                     colors=c("1"="#ff7f00", "2"="#e377c2", "3"="#17becf"),
                     pad=0.1, cex=0.2, pch=19, add=FALSE, legend.suffix="",
                     cex.main=1, cex.legend=1, legend = F, subsample=NULL) {

  layout = x
  if (is(x, "umap")) layout = x$layout

  if(!is.null(subsample)) {
	indices = tapply(1:nrow(layout), labels, function(x)
		sample(x, size = subsample * length(x)))
	layout = do.call(rbind, lapply(indices, function(x) layout[x, ]))
	labels = labels[ unlist(indices) ]
  }
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
   par(mar=c(0.2,0.7,1.2,0.7), ps=10)
   plot(xylim, xylim, type="n", axes=F, frame=F)
   rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=sapply(labels, function(x) {
    if(!x %in% names(colors)) return('grey')
    colors[[toString(x)]]
    }),
        cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  labels.u = unique(labels)
  legend.pos = "topright"
  legend.text = as.character(labels.u)
  if (add) {
   legend.pos = "bottomright"
   legend.text = paste(as.character(labels.u), legend.suffix)
  }
  if(legend)
    legend(legend.pos, legend=legend.text,
           col=sapply(labels.u, function(x) {
             if(!x %in% names(colors)) return('grey')
             colors[[x]]
             }), bty="n", pch=pch, cex=cex.legend)
}
