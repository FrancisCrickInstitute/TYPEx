#source("/Users/angelom/Documents/projects/IMC/scripts/settings/settings.R")
#source("/Users/angelom/Documents/projects/IMC/scripts/settings/color_settings.R")

library(ComplexHeatmap)
library(beeswarm)
library(scales)

paletteLen=256

regFile="/Volumes/rubicon/tracerx/master_files/metadata/metadata.tracerx.txt"
typingDir='/Volumes/lab-swantonc/working/angelom/analyses/imc/publication'
rubiconDir='/Volumes/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/'

options(stringsAsFactors=F)
arg_parser=argparser::arg_parser("Summarize typing results")
add=argparser::add_argument
arg_parser=add(arg_parser, "--analysis1", default = "p1,mcsa", 
               help = "NextFlow run")
arg_parser=add(arg_parser, "--analysis2", default = "p2,mcsa", 
               help = "Panel of markers")
args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))

imgData=read.csv(regFile, sep="\t", check.names = F)
imgData = subset(imgData, tx100_manuscript_filter &
                  ! tissue_type %in% c("Benign tumour-adjacent", 'LN'))

outDir = with(args, file.path('compare_intensity_by_poscell'))
if(!dir.exists(outDir))
  dir.create(outDir, recursive = T)

args$panel = strsplit(args$analysis1, split = ',')[[1]][1]
resultPattern = f("tx100_cell_objects.*_{args$panel}.fst")
analysis1=load_files(inDir = rubiconDir, 
                     pattern = f(resultPattern), 
                     recursive = F)
inDir1 = with(args, f("{typingDir}/{panel}/output/", 
                      "features/MeanIntensity"))
exp1 = load_files(inDir = inDir1, pattern = "*fst", recursive = F)
exp1 = setDT(exp1)

clinMatch = match(exp1$imagename, imgData$imagename)
exp1$RegionID_T = imgData$RegionID_T[clinMatch]
exp1 = subset(exp1, imagename %in% imgData$imagename)
anaMatch1 = match(with(exp1, paste(imagename, ObjectNumber)),
                  with(analysis1, paste(imagename, object)))
exp1$positive = analysis1$positive[anaMatch1]
exp1$cellType = analysis1$cellType[anaMatch1]
exp1$majorType = analysis1$majorType[anaMatch1]
print(dim(exp1))
print(length(unique(exp1$imagename)))

args$panel = strsplit(args$analysis2, split = ',')[[1]][1]
resultPattern = f("tx100_cell_objects.*_{args$panel}.fst")
inDir2 = with(args, f("{typingDir}/{panel}/output/",
                      "features/MeanIntensity"))
analysis2 = load_files(inDir = rubiconDir,
                       pattern = f(resultPattern),
                       recursive = F)
exp2 = load_files(inDir = inDir2, 
                  pattern = "*fst", 
                  recursive = F)
exp2 = setDT(exp2)
clinMatch = match(exp2$imagename, imgData$imagename)
exp2[, RegionID_T := imgData$RegionID_T[clinMatch]]
exp2 = subset(exp2, imagename %in% imgData$imagename)
anaMatch2 = match(with(exp2, paste(imagename, ObjectNumber)),
                  with(analysis2, paste(imagename, object)))
exp2$positive = analysis2$positive[anaMatch2]
exp2$cellType = analysis2$cellType[anaMatch2]
exp2$majorType = analysis2$majorType[anaMatch2]

markersCommon = setdiff(union(colnames(exp1), colnames(exp2)), 
                      c("file", "imagename", "ObjectNumber", "RegionID_T",
                        'positive', 'majorType', 'markerPos', 'cellType'))

for(marker in markersCommon) {
    
  exp1$markerPos = get_marker_frequency(data = exp1, marker, 'positive')
  exp2$markerPos = get_marker_frequency(data = exp2, marker, 'positive')
  
  # get average intensity per image
  markerCols = intersect(colnames(exp1),markersCommon)

  avg1 = exp1[1:100, lapply(.SD, c(median, length)), .SDcols = markerCols,
              by = .(imagename, majorType, markerPos,)]
  avg1$panel = 'p1'
                            # majorType == 'Epithelial cells' &
                            # markerPos %in% c("GZMB+", "LAG3+", "TIM3+", "PD1+", "CD103+")))
  print(head(avg1))
  # get average intensity per image
  markerCols = intersect(colnames(exp2), markersCommon)
  avg2 = exp2[, lapply(.SD, median), .SDcols = markerCols,
              by = .(imagename, majorType, markerPos)]
  avg2 = subset(avg2, ! ( majorType %in% c('T cells - Other', 'Ambiguous', 'Leukocytes - Other') ))
                            # majorType == 'Epithelial cells' &
                            # markerPos %in% c("GZMB+", "LAG3+", "TIM3+", "PD1+", "CD103+", 'GITR', 'FOXP3')))
  avg2$panel = 'p2'
  print(head(avg2))
  
  # What is the LAG3 intensity in images in P1 with many LAG3+
  cmbDF = rbindlist(list(avg1, avg2), use.names = T, fill = T)
  cmbDF = reshape2::melt(cmbDF, id.vars = c('panel', "imagename",
                                            "majorType", "markerPos"))
  cmbDF = dcast(cmbDF, RegionID_T + majorType + markerPos+ variable ~ panel)
  # cmbDF = na.omit(cmbDF)
  # cmbDF$variable = droplevels(cmbDF$variable)
  cmbDF$p1 = cmbDF$p1 * 10 ** 5
  cmbDF$p2 = cmbDF$p2 * 10 ** 5
  cmbDF$p1[cmbDF$p1 < 1e-2] = 1e-2
  cmbDF$p2[cmbDF$p2 < 1e-2] = 1e-2
  cmbDF= subset(cmbDF, majorType != 'Ambiguous')
  
  OtherSel= ! sapply(unique(cmbDF$majorType), function(x) any(grepl('\\+', cmbDF$markerPos[cmbDF$majorType == x])))
  # cmbDF$majorType[! cmbDF$majorType %in% names(OtherSel)[OtherSel]] = "Other"
  nrCelltypes = unique(cmbDF$variable)
  pdfOut = f("{outDir}/median_intensities_per_positive_type.log10.{marker}.pdf")
  pdf(pdfOut, height = 8, width = 8)
  sub = subset(cmbDF, variable == marker)
  # sub = na.omit(sub)
  sub = droplevels(sub)
  for(panel in c("p1", "p2")) {
    if(! marker %in% sub$variable) next
    if(all(is.na(sub[sub$variable == marker, panel])))
      next
    g <- ggplot(sub,  aes(majorType, get(panel)))
    plot = g + # geom_smooth(method='lm', color='grey20', se = F) +
     geom_violin(aes(fill = markerPos), # varwidth = T, 
                 position = position_dodge(width = 1), 
                  color='black') +
      geom_quasirandom(aes(fill = markerPos), dodge.width = 1, color='black', alpha = .5, pch = 21) +
      # geom_line(aes(group = RegionID_T))
      stat_summary(geom = 'pointrange', aes(fill = markerPos), # varwidth = T,
                   position = position_dodge(width = 1),
                   color='black', fun = median) +
      theme_classic(base_size = 16) + 
      scale_fill_manual(values = c('blue', 'red')) +
      theme(legend.position="top",
            axis.text.x=element_text(color="black", angle = 90, hjust = 1),
            axis.text.y=element_text(color="black"),
            strip.background=element_blank(),
            legend.title = element_blank()) +
      scale_y_log10() +
      # guides(fill = FALSE) + #guide_legend(title="Cell Types")) + 
      # scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
      # scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
      facet_wrap(. ~ confidence,  scale = 'free', nrow=4) + #
      ylab(panel) + xlab("") +
      ggtitle(marker) +
      # xlab(expression("Median Cell Intensity"  ~ "[ Pan-Immune panel" ~ "cells /" ~ mm^2 ~ "]")) + 
      # ylab(expression("Median Cell Intensity"  ~ "[ T cell-Stroma panel" ~ "cells /" ~ mm^2 ~ "]")) + 
      theme(strip.text=element_text(angle=0))
    print(plot)
  }
  dev.off()
}
