#!/usr/bin/env Rscript
source(glue::glue(Sys.getenv("BASE_DIR"), '/conf/settings.R'))


arg_parser=argparser::arg_parser("Summarize typing results")
add=argparser::add_argument

arg_parser=add(arg_parser, arg="--inDir", default='output', help='Where typing results are')
arg_parser=add(arg_parser, arg="--subset", default="sampled",
               help=paste0("The output typing directory [", majorDir, sampledDir, "]"))
arg_parser=add(arg_parser, arg="--method", default="Rphenograph",
               help=paste("Methods to be compared"))
arg_parser=add(arg_parser, "--run", default="final", help="NextFlow run")
arg_parser=add(arg_parser, "--panel", default="p2", help="Panel of markers")
arg_parser=add(arg_parser, "--markers", default="major_markers",
               help="Marker lists defined in TME_settings.R")
arg_parser=add(arg_parser, "--subtype_markers", default="major_markers",
               help="Marker lists defined in TME_settings.R")
arg_parser=add(arg_parser, arg="--ntasks", default=8,
               help="Number of threads")
## REFERENCE if sampled or full
arg_parser=add(arg_parser, "--reference_subset", default="major",
               help=paste0("The output typing directory [", majorDir, sampledDir, "]"))
arg_parser=add(arg_parser, "--reference_method", default="Rphenograph", help="Panel of markers")
arg_parser=add(arg_parser, "--reference_markers", default="major_markers",
               help="Marker lists defined in TME_settings.R")
args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))


downsample = 500000
transformation = "all"
runIDPattern = ".clusters.txt" #".full.RData"

outDir=with(args, file.path(inDir, subset, "robustness")) 

ref=with(args, list(subset=reference_subset, method=reference_method,
      markers=reference_markers, subtype_markers=reference_markers, run=run))
ids=c('subset', 'method', 'markers', 'run')
uniqueID=sapply(ids, function(x) ref[[x]] != args[x])
analysis_id=with(args, paste(c(args[unlist(ids)], ref[ids[uniqueID]]), collapse = "_"))


load(file.path(outDir,
               paste("matchedLabels", analysis_id, "RData", sep = ".")))

print(ls())
print(dim(clusters))
# Show down-sampled data to represent each cluster
downsizedIndices = unlist(tapply(1:nrow(clusters), clusters[,1], function(x) 
  sample(x, size = min(downsample, length(x)), replace = F )))
clusterMat = do.call(cbind, lapply(results, function(x)  x$cluster[downsizedIndices]))
print("Clustering")
dendro = fastcluster::hclust.vector(clusterMat, method = "centroid")

cols = c(tol21rainbow, color_clusters, rainbow(n = length(unique(clusterMat[,1])) - 51))
colorMat = WGCNA::labels2colors(clusters, colorSeq = cols)[, 1:length(runID_files)]
pdf(file.path(outDir, 
              paste0("Comparison.,",analysis_id, ".pdf", collapse = ".")))
plotDendroAndColors(dendro = dendro, 
                    colorMat[downsizedIndices, ],
                    names(results),
                    main = "",
                    autoColorHeight = FALSE, 
                    colorHeight = 0.65,
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    guideHang = 0.05,
                    addGuide = FALSE,
                    guideAll = FALSE, 
                    cex.main = 1, cex.lab = 1, 
                    cex.colorLabels = 0.6, 
                    marAll = c(8, 3, 3, 5))
# ddply(results[[1]], .(cluster, names), summarize, count = length(cluster))
labels = tapply(colorMat[, 1], results[[1]]$cluster, unique)
par(mar = c(0,0,0,0))
plot(1:length(labels), 1:length(labels), 
     type = "n", axes = F, xlab = "", ylab = "")
labelNamesFull = sapply(names(labels), 
                        function(x) results[[1]]$names[results[[1]]$cluster== x][1])
labelNames = gsub("up:(.*) down.*", "\\1", labelNamesFull)
# labelNames[labelNames == ""] =  gsub("up:.* (down.*)", "\\1", labelNamesFull)[labelNames == ""]
legend("top", legend = labelNames, cex = 0.5, bty = "n",
       ncol = 3,horiz = F, pt.cex = 1.5, x.intersp = 1, y.intersp = 1,
       pch = 20, col = labels, xpd=NA, title = names(results)[1])
run = 3 #length(runID_files)
labels = tapply(colorMat[, run], clusters[, run], unique)
labelNamesFull = sapply(names(labels), function(x) {
  origCluster = results[[run]]$cluster[which(clusters[, run] == x)[1]]
  results[[run]]$clus[results[[run]]$cluster == toString(origCluster)][1]
})
labelNames = gsub("up:(.*) down.*", "\\1", labelNamesFull)
# labelNames[labelNames == ""] = gsub("up:.* (down.*)", "\\1",
# labelNamesFull)[labelNames == ""]
legend("bottom", legend = labelNames, cex = 0.5, bty = "n",
       ncol = 3,horiz = F, pt.cex = 1.5, x.intersp = 1, y.intersp = 1,
       pch = 20, col = labels, xpd=NA, title = names(results)[run])
dev.off()

# Plot number of clusters
pdf(file.path(outDir, paste0("Cluster_stats", analysis_id, "pdf", sep = ".")))
# Plot number of cells per cluster (using the same colors), 
resultMrg = do.call(rbind, lapply(results, function(x) cbind(x$cluster, x$uniqueID)))
g <- ggplot(resultMrg, aes(as.factor(cluster), fill = uniqueID))
plot = g + geom_bar(position = "dodge", color = "darkgrey") + 
  theme_classic() +
  xlab("Cluster") + ylab("# Cells") +
  theme(legend.position = "top", legend.title = element_blank()) +
  scale_fill_manual(values = brewer.pal(length(runID_files), "Set1"))
print(plot)
# using matched cluster names

# calculate F1 score
colnames(clusters) = names(results)
precision = data.frame(t(sapply(colnames(clusters)[-1], function(col) {
  tp = sum(clusters[,col] == clusters[, 1], na.rm = T)
  p = tp / sum(!is.na(clusters[,1]))
  return(c(p = p, id = col))
})), check.names = F)
precision$p = as.numeric(as.character(precision$p))
g <- ggplot(precision, aes( as.factor(id), y = p))
plot = g + geom_bar(stat = "identity", 
                    fill = "transparent", col = "black") +
  theme_classic(base_size = 18) +
  theme(axis.text = element_text(color = "black")) +
  xlab("") + ylab("Precision score")
print(plot)
dev.off()
