#!/usr/bin/env Rscript

library(tidyr)
library(ggplot2)
library(data.table)
library(plyr)

source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/utilities.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/imc_utils.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/conf/settings.R'))

# Overlay function

# TODO
# get marker names from input params and infer the name of the model with one celltype excluded
# 

# MeanIntensity 0 for some cells (Mesenchymal)
arg_parser=argparser::arg_parser("Summarize typing results")
add=argparser::add_argument
arg_parser=add(arg_parser, "--celltypeReviewFile", 
		help='File with image annotations')
arg_parser=add(arg_parser, arg="--wDir", help="output")
arg_parser=add(arg_parser, arg="--cellReviewDir", 
	help=paste("review"))
arg_parser=add(arg_parser, arg="--subset", 
	default="major", help='major/sampled')
arg_parser=add(arg_parser, arg="--panel",
	 default='p1', help="Panel of markers")
arg_parser=add(arg_parser, arg="--major_method",
	 default="cellassign", help="Panel of markers")
arg_parser=add(arg_parser, arg="--markers", 
	default="major_markers", help="Panel of markers")
arg_parser=add(arg_parser, "--stratify_label", 
	default=NULL, help="file or NULL")

arg_parser=add(arg_parser, arg="--subtype_markers", 
	default="subtype_markers", help="Marker lists defined in TME_settings.R")
args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))

methods=c(args$major_method ,"csm")
features=c('area', 'meanIntensity', 'probability')
cellTypePattern=".review.fst"

subsetReview=read.delim(args$celltypeReviewFile, stringsAsFactors = F)
subsetReview=subset(subsetReview, control=="review")

if(! args$stratify_label %in% names(pars_stratify)) 
	stop("ERROR: parameters label ", args$stratify_label,
		"does not exist in conf/stratification.json")

args=c(args, pars_stratify[[args$stratify_label]])
csmCol=f("cellType_csm_{args$markers}")
regionCol=f("region_cellassign_{args$markers}")
markerSets=c(args$markers, args$excludedCellType)
cellAsCols=sapply(markerSets, 
	function(marker) f("cellType_cellassign_{marker}"))

abrv=strsplit(args$stratify_label, split = '')[[1]][1:3] %>% 
	tolower %>% paste0(., collapse = "")
ref_markers_list=paste1(pars$markers, abrv)
	
# IO
out=f("{args$wDir}/review/{args$panel}_{args$markers}")
visDir=f("{out}/overlay")
analysisID=with(args, f("{panel}_{markers}"))
modelOut=with(args, f("{cellReviewDir}/{analysisID}.RData"))

if(! dir.exists(out))
	dir.create(out, recursive = T)
if(! dir.exists(visDir)) 
	dir.create(visDir, recursive =T)

print('Merging files')
data=vector(mode="list")
for(markers in markerSets) {
  for(method in methods) {
    wdir=f("{args$wDir}/", with(args, f(analysisPath)))
    if(method == "csm" & length(grep(method, names(data)))) 
		next
    revFile=list.files(wdir, pattern=cellTypePattern, full.names = T)
	print(wdir)
	print(cellTypePattern)
    if(! length(revFile)) 
		next
    tmp=fst::read_fst(revFile)
    cat(markers, method, length(unique(tmp$imagename)), '\n')
	if(method == "csm")
		tmp$positive=tmp$names
	
	# Only keep CD3+CD4+ T cells as a positive control for CD4 T cells
    tmp$CD3=get_marker_frequency(marker = 'CD3', data = tmp, column='positive')
    tmp$cellType[tmp$cellType == 'CD4 T cells' & grepl('\\-', tmp$CD3)] = NA
	tmp$cellType=gsub(" - Other", "", .) %>% 
		gsub('Smooth muscle cells', "Myofibroblasts", .)
   
    data[[paste1(method, markers)]]=data.table::as.data.table(cbind(tmp, method=method, markers=markers))
  }
}
if(! length(data)) {
	print('WARNING: Skipping stratification by confidence. 
		No images for manual review selected.')
} else {
	mrg=data.table::rbindlist(data, use.names = T, fill = T)
	# Exclude unassigned, esp from csm, so that there are not considered as false negatives
	mrg=mrg[!mrg$cellType %in% c('Excluded', 'Unassigned'), ]
	mrg$cellID=with(mrg, paste(object, imagename))
	mrg$centerX=round(mrg$centerX, 2)
	mrg$centerY=round(mrg$centerY, 2)

	for(feature in grepv("Intensity", features))
	  mrg[, (feature) := lapply(.SD, function(x) log2(to_magnitude(x) + 1)), .SDcols=feature]
	smoothMuscle=which(mrg$cellType == "Smooth muscle cells")
	if(length(smoothMuscle)) mrg$cellType[smoothMuscle] = "Myofibroblasts"
	mrg=mrg[mrg$imagename %in% subsetReview$imagename, ]
	recast<-data.table::dcast(mrg, cellID + panel + imagename + 
	                            centerX + centerY + object ~ method + markers,
	                          value.var=c("cluster", "cellType", "positive"))
	print(head(recast))
	# Those in csm that do not have any associated raw masks
	recast$cellType_csm_mcsa[is.na(recast$cellType_csm_mcsa)] = 'None'

	imagenames=unique(recast$imagename)
	# Exclude the cells that have been correctly called by CSM in the negative control
	negative=setdiff(which(recast[[csmCol]] == recast[[cellAsCols[1]]] &
	                         recast[[cellAsCols[1]]] != recast[[cellAsCols[2]]] |
	                         recast[[csmCol]] == 'None' & 
	                         recast[[cellAsCols[1]]] %in% args$mostLikelyCellType_full & 
	                         recast[[cellAsCols[2]]] %in% args$mostLikelyCellType_ref),
	                 which(recast[[cellAsCols[1]]] %in% c(args$excludedCellType,args$dependentCell) &
	                         recast[[cellAsCols[2]]] %in% args$dependentCell))
	if(toupper(args$panel)=='P2')
	  negative = union(negative,
	                   which(recast[[csmCol]] != args$mostLikelyCellType_full & 
	                           recast[[cellAsCols[1]]] == args$mostLikelyCellType_full &
	                           recast[[cellAsCols[2]]] == args$mostLikelyCellType_ref & #recast[[cellAsCols[2]]]
	                           recast$imagename %in% review$imagename[subsetReview$description=='Carcinosarcoma']))
						   
	# Exclude low-intensity lymphocytes that might have been assigned to Mesenchymal
	negative=recast$cellID[negative]
	positive=recast$cellID[which(recast[[csmCol]] == recast[[cellAsCols[2]]] &
	                             recast[[csmCol]] == recast[[cellAsCols[1]]])]
	if(regionCol %in% colnames(recast))
	  positive=union(positive, recast$cellID[which(recast[[regionCol]] == "Tumour" & recast[[cellAsCols[2]]] == "Epithelial cells")])
	positiveExcluded=recast$cellID[which(recast[[cellAsCols[1]]] %in% args$excludedCellType &
	                                     recast[[cellAsCols[2]]] %in% c(args$dependentCell, args$mostLikelyCellType_ref) &
	                                     recast[[csmCol]] %in% c(args$dependentCell, args$excludedCellType))]
	negativeExcluded=recast$cellID[which(recast[[csmCol]] != recast[[cellAsCols[1]]] &
	                                     recast[[cellAsCols[1]]] %in% args$excludedCellType &
	                                   ! recast[[cellAsCols[2]]] %in% c(args$mostLikelyCellType_ref, args$dependentCell))]
	#recast$cellType_cellassign_mcsamy[recast$cellID %in% grepv('P2_TMA001_L_20190508-roi_11', negativeExcluded)] %>% table
	totalargs$excludedCellType=mrg$markers==args$markers & mrg$method == "cellassign" &
	  !is.na(mrg$probability) & mrg$cellType %in% args$excludedCellType
	posExcludedIndices=mrg$cellID %in% positiveExcluded & totalargs$excludedCellType
	negExcludedIndices=mrg$cellID %in% negativeExcluded & totalargs$excludedCellType

	print('Negative')
	paste(recast[[cellAsCols[1]]][recast$cellID %in% negative],
	                       recast[[cellAsCols[2]]][recast$cellID %in% negative]) %>%
						   		table %>% sort %>% print
	print('Positive')
	paste(recast[[cellAsCols[1]]][recast$cellID %in% positive],
	                       recast[[cellAsCols[2]]][recast$cellID %in% positive]) %>% 
						   table %>% sort %>% print
	print('positiveExcluded')
	recast$names_csm_mcsa[recast[[cellAsCols[1]]] %in% args$excludedCellType &
						  recast$cellID %in% positiveExcluded] %>% table %>% 
						  	sort(., decreasing = T) %>% head %>% print 
	print('Negative excluded')
	recast$names_csm_mcsa[recast[[cellAsCols[1]]] %in% args$excludedCellType &
						  recast$cellID %in% negativeExcluded] %>% table %>%
						 	 sort(., decreasing = T) %>% head %>% print

	mrg$control="undetermined"
	mrg$control[which(mrg$cellID %in% negative)] = "negative"
	mrg$control[which(mrg$cellID %in% positive)] = "positive"
	mrg$control[totalargs$excludedCellType]="undetermined"
	mrg$control[negExcludedIndices]="negative"
	mrg$control[posExcludedIndices]="positive"
	print(table(mrg$control))
	# mrg$control[mrg$cellType == "Excluded"] = "negative"

	# mrg$cellType[smoothMuscle]="Smooth muscle cells"
	# Predict positive and negative cell objects
	controlCellIDs=mrg$control %in% c("positive", "negative") &
	  mrg$markers==ref_markers_list & mrg$method == "cellassign" &
	  !is.na(mrg$probability)
	print('Probability')
	print(table(mrg$cellType[is.na(mrg$probability & mrg$method == "cellassign")]))

	dfFlt=mrg[controlCellIDs | posExcludedIndices | negExcludedIndices, ]
	dfFlt$control=sapply(dfFlt$control == "positive", sum)

	# Model on selected covariates
	print('Creating model with')
	print(table(dfFlt$control, dfFlt$cellType))
	model<-lme4::glmer(as.formula(f("control ~ probability + meanIntensity +
	                                (1|cellType)")), data=dfFlt, family=binomial)
	predict=predict(model, newdata=dfFlt, type="response", allow.new.levels=T)
	ROCRpred <- ROCR::prediction(predict[!is.na(predict)], dfFlt$control[!is.na(predict)])
	sensit=ROCR::performance(ROCRpred, 'sens')
	specif=ROCR::performance(ROCRpred, 'spec')
	# cat(sensit, specif, "\n")
	modelCutoff=sensit@x.values[[1]][which.min(abs(sensit@y.values[[1]] - specif@y.values[[1]]))]

	save(modelCutoff, model, file = modelOut)
	sensDf <- data.frame(x=unlist(sensit@x.values), y=unlist(sensit@y.values))
	specDf <- data.frame(x=unlist(specif@x.values), y=unlist(specif@y.values))
	pdf(f("{out}/Model.{analysisID}.pdf"), height=3, width=4)
	plot=sensDf %>% ggplot(aes(x,y)) +
	  geom_line() +
	  geom_line(data=specDf, aes(x,y,col="red")) +
	  scale_y_continuous(sec.axis = sec_axis(~., name = "Specificity")) +
	  labs(x='Cutoff', y="Sensitivity") +
	  theme(axis.title.y.right = element_text(colour = "red"), legend.position="none") +
	  geom_vline(aes(xintercept=modelCutoff), linetype="dashed") +
	  guides(color='none') +
	  cowplot::theme_cowplot()
	print(plot)
	dev.off()

	mrg$predicted=predictProb > modelCutoff
	mrg$predicted=ifelse(mrg$predicted, "positive", "negative")
	dfMrg=droplevels(mrg[mrg$method == "cellassign" & mrg$markers == ref_markers_list |
				posExcludedIndices | negExcludedIndices, ])
	print(table(dfMrg$cellType, dfMrg$predicted))

	stats=data.frame(table(dfMrg$control))
	stats$Freq=round(stats$Freq/sum(stats$Freq), 3)
	stats$Var1=factor(stats$Var1, levels=sort(stats$Var1))

	pdf(f("{out}/Control_stats.{analysisID}.pdf"), height=5, width=5)
	g <- ggplot(stats, aes(x="", y=Freq, fill=Var1))
	plot=g + geom_bar(stat="identity") + theme_classic(base_size=18) +
	  scale_fill_manual(name="Cell type controls", 
	                    values=palette$cellTypingStatusCols,
	                    labels=paste0(levels(stats$Var1), " (", 
	                                  stats$Freq[match(levels(stats$Var1), stats$Var1)] * 100, "%)")) +
	  coord_polar(theta="y", start=0, direction=-1) +
	  theme(axis.line=element_blank(),
	        axis.text= element_blank(),
	        axis.ticks=element_blank()) +
	  xlab("") + ylab("")
	  # scale_y_continuous(labels=function(x) paste0(x * 100, "%"))
	print(plot)

	predStats=data.table(table(dfMrg$control, dfMrg$predicted))
	predStats$Freq=sapply(1:nrow(predStats),  function(x) 
	  round(predStats$N[x]/sum(predStats$N[predStats$V1 == predStats$V1[x]]), 3))
	predStats$Estimated=ifelse(predStats$V2 == predStats$V1, "Predicted", NA)
	predStats=predStats[predStats$Estimated != "Inc",]
	predStats=predStats[predStats$N > 0, ]
	predStats$V1=factor(predStats$V1, levels=unique(sort(predStats$V1)))
	g <- ggplot(predStats, aes(x="", y=Freq))
	plot=g + geom_bar(stat="identity",fill = "lightblue", color="black" ) +
	  theme_classic(base_size=18) +
	  facet_grid(. ~ V1) +
	  coord_polar(theta="y", start=0, direction=-1) +
	  theme(axis.line=element_blank(), axis.ticks=element_blank(), 
	        axis.text=element_blank()) +
	  xlab("") + ylab("") +
	  scale_y_continuous(limits = c(0, 1))
	print(plot)

	cellStats=data.table(table(dfMrg$control, dfMrg$cellType))
	# cellStats$Freq=round(cellStats$N/sum(cellStats$N), 3)
	cellStats$Freq=sapply(1:nrow(cellStats), function(x) with(cellStats, round(N[x]/sum(N[V1 == V1[x]]), 3)))
	cellStats$V1=factor(cellStats$V1, levels=unique(sort(cellStats$V1)))
	cellStats$V2=factor(cellStats$V2, levels=names(palette$cellTypeColors)[names(palette$cellTypeColors) %in% cellStats$V2])
	g <- ggplot(cellStats, aes(x=V1, y=Freq, fill=V2))
	plot=g + geom_bar(stat="identity", color="black" ) +
	  theme_classic(base_size=18) +
	  theme(axis.line=element_blank(), axis.ticks=element_blank(), 
	        axis.text.x=element_text(angle=45)) +
	  xlab("") + ylab("") +
	  scale_fill_manual(values=palette$cellTypeColors) +
	  guides(fill=guide_legend(title='Cell type'))
	# scale_y_continuous(limits = c(0, 1))
	print(plot)
	dev.off()
	# browseURL(f("{out}/Control_stats.{analysisID}.pdf"))

	pdf(f("{out}/Control_comparison.{analysisID}.pdf"), height=4)
	for(feature in features) {

	  dfMrg=droplevels(mrg[mrg$method == "cellassign" & 
	  					   mrg$markers == ref_markers_list |
	                       mrg$cellType %in% args$excludedCellType &
						   mrg$markers== args$markers, ]
	  )
	  # dfMrg=dfMrg[dfMrg$cellType!= args$excludedCellType, ]
	  cat(feature, nrow(dfMrg), "\n")
	  cellTypeStats=table(dfMrg$cellType)
	  print(cellTypeStats)
	  g <- ggplot(dfMrg)
	  plot = g + 
	    cowplot::theme_cowplot(font_size = 14) + # ggtitle(feature) +
	    theme(axis.text.x=element_text(angle=45, hjust = 1),
	          legend.position="top") +
	    guides(fill=guide_legend(nrow=length(unique(mrg$method)))) +
	    xlab(NULL)
	    # facet_grid(method ~ markers, scale="free_x") 
	    if(feature=="Probability")
	        plot = plot + geom_hline(aes(yintercept=0.5), linetype="dashed", color="grey")
	    if(feature == "area")
	      plot = plot + scale_y_log10()

	  plot1=plot + geom_boxplot(aes_string("cellType", feature, fill="control"), outlier.size = 0.05) +
	    scale_fill_manual(values = palette$cellTypingStatusCols)
	  print(plot1)
	  plot2=plot + geom_boxplot(data=droplevels(subset(dfMrg, !is.na(predicted))),
	                            aes_string("cellType", feature, fill="predicted")) +
	    scale_fill_manual(values=palette$cellTypingStatusCols) +
	    facet_grid( . ~ control)
	  print(plot2)
	  plot3=plot + geom_density(aes_string(feature, color="control"), fill="grey", alpha=0.2) +
	    cowplot::theme_cowplot(font_size = 22) +
	    scale_color_manual(values=palette$cellTypingStatusCols) +
	    ylab("Density")
	  print(plot3)
	}
	dev.off()
	# browseURL(f("{out}/Control_comparison.{analysisID}.pdf"))
	for(imagename in rev(imagenames)) {

	  pdfOut=f("{visDir}/{imagename}.{analysisID}.pdf")
	  # pdf(pdfOut, width=10, height=11)
	  for(analysis in grepv("cellType_", colnames(recast))) {
		  
	    method=gsub("cellType_([^_]+)_.*", "\\1", analysis)
	    if(method=="csm") next
	    markers=gsub("cellType_[^_]+_(.*)", "\\1", analysis)
	    rowSel=mrg$imagename==imagename & mrg$method == method & mrg$markers == markers
	    if(!sum(rowSel)) next
    
	    for(region in c("negative", "positive")) {
		
	      for(column in c("predicted", "control")) {
			  
	        pngOut=f("{visDir}/{region}_{column}.{method}.{imagename}.{analysisID}.{markers}.png")
	        png(pngOut,  res = 400, units = "in", width=5.73, height=6)
	        regionSel=rowSel & mrg[[column]] == region
	        overlay_cell_types(imagename, imageDf=mrg[regionSel, ],
	                           cellTypeCol="cellType", cellTypes=column,
	                           title=f("{imagename} {region} {column} {markers}"), 
	                           legend = T, ptx=0.1)
	        dev.off()
	      }
	    }
	  }
	}
}
	

