library(wavelets)
args<-commandArgs(TRUE)
dat = read.csv("data/testdata.csv",header=F)
temp = dwt(dat[1:512,1], filter=args[1],n.levels = as.integer(args[2]), boundary = args[3])
write.table(do.call(cbind, lapply(temp@W, function(i) c(i, rep(NA, 512/2-length(i))))), file="W.csv", row.names=F, col.names=F, sep=", ", na = "nan")
write.table(do.call(cbind, lapply(temp@V, function(i) c(i, rep(NA, 512/2-length(i))))), file="V.csv", row.names=F, col.names=F, sep=", ", na = "nan")