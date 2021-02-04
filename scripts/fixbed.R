#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


bed <- read.table(args[1], header = F, sep = '\t')

for(i in 1:nrow(bed)){
 bed$V4[i] <- min(bed$V2[i], bed$V3[i])
 }
for(i in 1:nrow(bed)){
 bed$V5[i] <- max(bed$V2[i], bed$V3[i])
 }

sortbed <- bed[c(1,4:5)]
write.table(sortbed, file=args[2], sep = '\t', col.names = F, row.names=F, quote = F)
