source("scripts/comparison/settings.R")

downsample = 500000
to_compare = "Rphenograph"
transformation = "all"
sampling = "FALSE"
runIDPattern = ".clusters.txt" #".full.RData"

analysis_id = paste(to_compare, transformation, sampling, sep = "_")
if(sampling) {
  outputDir = "/camp/home/angelom/labwd/analyses/sampled/" # compare paramters
} else {
  outputDir = "/camp/home/angelom/labwd/analyses/full/" # compare paramters
}# outputDir = "/Users/angelom/Documents/projects/IMC/full/testing_runs" # compare paramters

load(file.path(outputDir,
               paste("matchedLabels", analysis_id, "RData", sep = ".")))

# Show down-sampled data to represent each cluster
downsizedIndices = unlist(tapply(1:nrow(clusters), clusters[,1], function(x) 
  sample(x, size = min(downsample, length(x)), replace = F )))
clusterMat = do.call(cbind, lapply(results, function(x)  x$cluster[downsizedIndices]))
print("Clustering")
dendro = fastcluster::hclust.vector(clusterMat, method = "centroid")

cols = c(tol21rainbow, color_clusters, rainbow(n = length(unique(clusterMat[,1])) - 51))
colorMat = WGCNA::labels2colors(clusters, colorSeq = cols)[, 1:length(runID_files)]
pdf(file.path(outputDir, 
              paste0("Comparison.,",analysis_id, ".pdf", collapse = ".")))
plotDendroAndColors(dendro = dendro, 
                    colorMat[downsizedIndices, ],
                    names(results),
                    main = "",
                    # main = "Gene dendrogram and module labels from resampled data sets",
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