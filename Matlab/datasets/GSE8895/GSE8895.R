library(GEOquery)
gset <- getGEO("GSE8895", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL90", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
ex <- exprs(gset)
#symbols = fData(gset)[,'Gene symbol']
#rownames(ex) <-  fData(gset)[,'Platform_ORF']

library("data.table")
t = data.table(ex)
t = cbind(fData(gset)[,'Platform_ORF'], t)
t = t[!(is.na(t$V1) | t$V1==""), ]

library("xlsx")
setwd("D:/MyWorks/ssGSEA/Matlab/Yeast/datasets/GSE8895")
write.xlsx(t, 'GSE8895.xlsx', row.names = FALSE)