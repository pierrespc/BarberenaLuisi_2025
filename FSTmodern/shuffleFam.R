#!/bin/Rscript


fam=commandArgs(trailingOnly=T)
a<-read.table(fam,stringsAsFactors=F,header=F)


a$V1<-sample(a$V1,nrow(a))

write.table(a,fam,col.names=F,row.names=F,sep=" ",quote=F)
