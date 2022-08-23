#!/usr/bin/env Rscript

# Collate files by feature
# source("scripts/settings/settings.R")
source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/utilities.R'))
source(glue::glue(Sys.getenv("BASE_DIR"), '/lib/imc_utils.R'))


arg_parser=argparser::arg_parser("Collate image info by feature")
add=argparser::add_argument
arg_parser=add(arg_parser, "--inDir", default="tracerx", help="Panel of markers")
arg_parser=add(arg_parser, "--outDir", default="tx100", help="Panel of markers")
arg_parser=add(arg_parser, "--panel", default="p1", help="Panel of markers")
arg_parser=add(arg_parser, "--run", default="final", help="NextFlow run")
arg_parser=add(arg_parser, "--feature", default="LocationCenter", help="Panel of markers")

args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))
args[["channels_exclude"]]= c("Argon", "131Xe", "134Xe", "DNA1", "DNA2", "processed")

# outDir=f("{labwd}/analyses/imc/cell_objects/{args$study}/{args$cohort}")

if(!dir.exists(args$outDir)) dir.create(args$outDir, recursive=T)

for(feature in args$feature) {

  inData=load_files(inDir=file.path(args$inDir, feature),
                    run=args$run, col.exclude = args$channels_exclude)
  fileOut=file.path(args$outDir, paste(args$panel, feature, args$run, sep = "_"))
  inData$file=basename(inData$file)
  inData$file=gsub(".txt", "", inData$file)
  colnames(inData)[colnames(inData) == "file"] = "imagename"
  write.table(inData, file  = f("{fileOut}.csv"), sep = ",", quote  = F, row.names = F)
  fst::write.fst(inData, path = f("{fileOut}.fst"), compress = 75)

}

