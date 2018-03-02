args<-commandArgs(trailingOnly=TRUE)
q<-read.table(args[1],header=FALSE,skip=1)
n=ncol(q)
q_new<-q[][5:n]
pdf(args[2])
barplot(t(as.matrix(q_new)),col=rainbow(args[3]),xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()
