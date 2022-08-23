#!/usr/bin/env Rscript
library(plyr)
source(file.path(Sys.getenv('BASE_DIR'), 'conf/settings.R'))
source(file.path(Sys.getenv('BASE_DIR'), 'lib/utilities.R'))

arg_parser=argparser::arg_parser("Run a cell typing method")
add=argparser::add_argument
arg_parser=add(arg_parser, "--outDir", help="NextFlow run")
arg_parser=add(arg_parser, "--tissAreaDir", help="Panel of markers")
arg_parser=add(arg_parser, "--panel", default='p2', help="Panel of markers")

args=argparser::parse_args(arg_parser, argv=commandArgs(trailingOnly=TRUE))

areaDir=f('{args$tissAreaDir}/area_info')
areaInfoFiles=list.files(areaDir)
summary=vector(mode='list')
for(file in areaInfoFiles) {
  cat(which(areaInfoFiles == file), file, '\n')
  data=read.csv(f('{areaDir}/{file}'))
  summary[[file]]=ddply(data, .(Label), summarise, 
        Count=length(Label),
        "Total Area"=sum(Area),
        AverageSize=mean(Area),
        X=mean(X), Y=mean(Y),
        Perim.=mean(Perim.), Minor=mean(Minor),
        Major=mean(Major),
        Angle=mean(Angle), Circ.=mean(Circ.), 
        AR=mean(AR), 
        Round=mean(Round), Solidity=mean(Solidity))
}

sumOut=do.call(rbind, summary)
sumOut$imagename=gsub('^[^_]+_(.*).tiff', '\\1', sumOut$Label)
sumOut$panel=tolower(gsub('_.*', '', sumOut$imagename))
sumOut$TissueType=gsub("([^_]+).*", "\\1", sumOut$Label)

fileOut=with(args, f("{tissAreaDir}/{panel}_tissue_area.csv"))
write.table(sumOut, file=fileOut, row.names = F, quote = F, sep = ',')



