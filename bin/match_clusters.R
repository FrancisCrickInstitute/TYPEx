#!/usr/bin/env Rscript

source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/utilities.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/conf/settings.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/imc_utils.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/matchLabels.R'))

resultPattern=".clusters.fst"

arg_parser=argparser::arg_parser("Summarize typing results")
add=argparser::add_argument

arg_parser=add(arg_parser, arg="--inDir", help='Where typing results are')
arg_parser=add(arg_parser, arg="--subset", default="sampled", help=paste0("The output typing directory [", majorDir, sampledDir, "]"))
arg_parser=add(arg_parser, arg="--method", default="FastPG", help=paste("Methods to be compared"))
arg_parser=add(arg_parser, "--run", default="publication", help="NextFlow run")
arg_parser=add(arg_parser, "--panel", default="p2", help="Panel of markers")
arg_parser=add(arg_parser, "--markers", default="major_markers", help="Marker lists defined in TME_settings.R")
arg_parser=add(arg_parser, "--subtype_markers", default="major_markers", help="Marker lists defined in TME_settings.R")
arg_parser=add(arg_parser, arg="--ntasks", default=8, help="Number of threads")
## REFERENCE if sampled or full
arg_parser=add(arg_parser, "--reference_subset", default="major",
               help=paste0("The output typing directory [", majorDir, sampledDir, "]"))
arg_parser=add(arg_parser, "--reference_method", default="FastPG", help="Panel of markers")
arg_parser=add(arg_parser, "--reference_markers", default="major_markers",
               help="Marker lists defined in TME_settings.R")
args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))

ref=with(args, list(subset=reference_subset, method=reference_method, 
      markers=reference_markers, subtype_markers=reference_markers, run=run))

inDir=file.path(args$inDir, with(args, f(analysisPath)))
inFiles=list.frec(inDir, pattern=resultPattern)

refDir=file.path(args$inDir, with(ref, f(analysisPath)))
refFile=list.frec(refDir, pattern=resultPattern)

outDir=file.path(args$inDir, args$subset, "robustness")
if(! dir.exists(outDir))  dir.create(outDir)
ids=c('subset', 'method', 'markers', 'run')
uniqueID=sapply(ids, function(x) ref[[x]] != args[x])
analysisID=with(args, paste(c(args[unlist(ids)], ref[ids[uniqueID]]), collapse = "_"))
cat('Reading ', refFile, "in", refDir, '\n')
refData=fst::read.fst(refFile)
cat("Matching labels", analysisID, '\n')
clusters=vector(mode='list')
clusters[[1]]=sapply(refData$cluster, toString)
pind=WGCNA::initProgInd()
results=vector(mode="list")
print(head(refData))
results[[refFile]]=refData[, c('imagename', 'object', 'cluster', 'cellType', 'positive')]

for(r in 1:length(inFiles)) {
  print(r)
  fileOut=f("{outDir}/matched.{analysisID}.{r}.txt")
  cat('Reading ', inFiles[r], '\n')
  iter=fst::read.fst(inFiles[r])
  iter$cellType[iter$cluster  == 'Excluded'] = NA
  iter$cluster[iter$cluster == 'Excluded'] = NA
  ind=match(with(refData, paste(imagename, object)),
            with(iter, paste(imagename, object)))
  results[[inFiles[r]]] = iter[ind, c('imagename', 'object', 'cluster', 'cellType', 'positive')]
  if(all(is.na(ind))) {
    print("No matching cells found")
    next
  }
  res=matchLabels(source=iter$cluster[ind],
                              reference=refData$cluster,
                              nThreads=args$ntasks,
                              ignoreLabels = "NA")
  clusters[[r+1]]=res
  pind=WGCNA::updateProgInd((r-1)/length(inFiles), pind)
  write.tab(clusters[[r+1]], file = fileOut)
}
save(clusters, results, inFiles, file=file.path(outDir, paste("matchedLabels", analysisID, "RData", sep=".")))
cat('Output saved in ', outDir, '\n')
