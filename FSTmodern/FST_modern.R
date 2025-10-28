#!/bin/Rscript
require(scales)
setwd("~/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/AdmixtureMasked/DiffAF/")

require(plotrix)


change<-function(string){
  return(strsplit(string,split="/")[[1]][2])
}
listDir<-sapply(list.dirs()[-1],change,USE.NAMES = F)

for(pop in listDir){
  tmp<-read.table(paste(pop,"/",pop,".frq",sep=""),stringsAsFactors =F, header=T)
  if(pop==listDir[1]){
    AFusp<-tmp
  }else{
    AFusp<-cbind(AFusp,tmp[,c(5,6)])
  }
  names(AFusp)[length(AFusp)+c(-1,0)]<-paste(pop,c("","_N"),sep="")
  
}


IndPops<-listDir[ grepl("delaFuente",listDir) | grepl("Reich",listDir)]
prop<-matrix(NA,length(IndPops),length(IndPops))
prop<-matrix(NA,length(IndPops),length(IndPops))
colnames(prop)<-rownames(prop)<-IndPops

DIFF<-ALL<-prop
for(p1 in c(1:(nrow(prop)-1))){
  pop1<-colnames(prop)[p1]
  for(p2 in c((p1+1):nrow(prop))){
    pop2<-colnames(prop)[p2]
    count<-sum(abs(AFusp[,pop2]-AFusp[,pop1])>0.2,na.rm = T)
    noNA<-sum(!is.na(AFusp[,pop2]-AFusp[,pop1]))
    
    prop[pop1,pop2]<-prop[pop2,pop1]<-count/noNA
    DIFF[pop1,pop2]<-DIFF[pop2,pop1]<-count
    ALL[pop1,pop2]<-ALL[pop2,pop1]<-noNA
  }
}



IndPops<-c("Aymara-Reich","Cabecar-Reich","Chono-Reich","Diaguita-Reich","Guarani-Reich","Huilliche-delaFuente","Kaweskar-delaFuente")
forBarplot<-matrix(NA,choose(length(IndPops),2),length(seq(0,0.9,0.1)))
colnames(forBarplot)<-paste(seq(0,0.9,0.1),"-",seq(0.1,1,0.1),sep="")

iter=0
listComp<-c()
listMore3<-c()
listBelow1<-c()
listPerc<-c()
for(p1 in c(1:(length(IndPops)-1))){
  for(p2 in c((p1+1):length(IndPops))){
    iter=iter+1
    pop1<-IndPops[p1]
    pop2<-IndPops[p2]
    comp<-paste(pop1,pop2,sep=" vs ")
  
    vector<-abs(na.omit(AFusp[,pop2]-AFusp[,pop1]))
    noNA<-length(vector)
    perc<-round(quantile(vector,0.95),digits = 4)
    listPerc<-c(listPerc,perc)
    sum<-0
    for(af in seq(0,0.9,0.1)){
      lower<-af
      upper<-ifelse(af==0.9,1.01,af+0.1)
      #print(c(lower,upper))
      tmp<-sum((vector<upper & vector >= lower ) | vector==lower)
      vector<-vector[ ! ((vector<upper & vector>=lower) | vector==lower)]
      forBarplot[iter,paste(af,"-",af+0.1,sep="")]<-tmp/noNA
      sum<-sum+tmp
    }
    
    more3<-round(sum(forBarplot[iter,paste(seq(0.3,0.9,0.1),seq(0.4,1,0.1),sep="-")]),digits=3)
    below1<-round(forBarplot[iter,"0-0.1"],digits=3)
    listMore3<-c(listMore3,more3)
    listBelow1<-c(listBelow1,below1)
    #print(vector)
    listComp<-c(listComp,paste(comp," (",below1," / ",perc," / ",more3,")",sep=""))
    forBarplot[iter,"0.6-0.7"]<-forBarplot[iter,"0.6-0.7"]+length(vector)/noNA
    print(c(sum(forBarplot[iter,]),noNA,sum))
  }
}

listComp<-str_remove_all(listComp,"-Reich")
listComp<-str_remove_all(listComp,"-delaFuente")
row.names(forBarplot)<-listComp


pdf("DiffAF_Indexes.pdf",width=10,height=5)
par(mfrow=c(1,3))
boxplot(listBelow1,main="B\nProportion of SNPs with |f(A in pop1) -  f(A in pop2)| < 0.1\nDistribution from all pairs of populations",cex.main=0.8)
boxplot(listPerc,main="C\n95% percentile fro |f(A in pop1) -  f(A in pop2)| \nDistribution from all pairs of populations",cex.main=0.8)
boxplot(listMore3,main="D\nProportion of SNPs with |f(A in pop1) -  f(A in pop2)| > 0.3\nDistribution from all pairs of populations",cex.main=0.8)
dev.off()
pdf("DiffAF.pdf")

par(mar=c(15,4,4,1))
barplot(t(forBarplot),las=2,cex.names = 0.6,legend.text = F,xlab = "",col=paste("grey",seq(100,10,length.out=10)),cex.lab=0.8,
        ylab="Proportion of SNPs with |f(A in pop1) - f(A in pop2)|",
        cex.axis=0.6,cex.main=0.8,
        main="A\nProportion of SNPs with\na given difference in Allele frequency between two populations")
plot(0,0,"n",axes=F,ann=F)
legend("center",pch=15,col=paste("grey",seq(100,10,length.out=10)),title="Observed |f(A in pop1) - f(A in pop2)|",
       #legend = paste("|f(A in pop1) - f(A in pop2)| in [",seq(0,0.9,0.1),"-",seq(1,1,0.1),c(rep("[",8),"]"))
       legend = paste("[",seq(0,0.9,0.1),"-",seq(1,1,0.1),c(rep("[",8),"]"),sep="")
       )
dev.off()




IndPops<-c("Aymara-Reich","Cabecar-Reich","Chono-Reich","Diaguita-Reich","Guarani-Reich","Huilliche-delaFuente","Kaweskar-delaFuente")

FSTbetweenModern<-list()
wilcoxTest<-c()

diffFST<-list()
forBarplot<-matrix(NA,choose(length(IndPops),2),length(seq(0,0.9,0.1)))
colnames(forBarplot)<-paste(seq(0,0.9,0.1),"-",seq(0.1,1,0.1),sep="")
iter=0
for(p1 in c(1:(length(IndPops)-1))){
  for(p2 in c((p1+1):length(IndPops))){
    iter=iter+1
    pop1<-IndPops[p1]
    pop2<-IndPops[p2]
    comp<-paste(pop1,pop2,sep=" vs ")
    comp<-str_remove_all(comp,"-Reich")
    comp<-str_remove_all(comp,"-delaFuente")
    tmp1<-read.table(paste(pop1,"/",pop1,".frq",sep=""),stringsAsFactors =F, header=T)
    tmp2<-read.table(paste(pop2,"/",pop2,".frq",sep=""),stringsAsFactors =F, header=T)
    
    merged<-merge(tmp1,tmp2,by="SNP")
    numHudson<-(merged$MAF.x-merged$MAF.y)^2 - 
      (merged$MAF.x)*(1-merged$MAF.x)/(merged$NCHROBS.x-1) -
      (merged$MAF.y)*(1-merged$MAF.y)/(merged$NCHROBS.y-1)
    denomHudson=merged$MAF.x+merged$MAF.y-2*merged$MAF.x*merged$MAF.y
    fstHudson<-ifelse(merged$NCHROBS.x > 2 & merged$NCHROBS.y>2,numHudson/denomHudson,NA)
    fstHudson[fstHudson<0]<-0
    FSTbetweenModern[[comp]]<-fstHudson
    
    sum<-0
    vector=na.omit(fstHudson)
    noNA<-length(vector)
    for(af in seq(0,0.9,0.1)){
      lower<-af
      upper<-ifelse(af==0.9,1.01,af+0.1)
      #print(c(lower,upper))
      tmp<-sum((vector<upper & vector >= lower ) | vector==lower)
      vector<-vector[ ! ((vector<upper & vector>=lower) | vector==lower)]
      forBarplot[iter,paste(af,"-",af+0.1,sep="")]<-tmp/noNA
      sum<-sum+tmp
    }
    
    
    
    tmp1<-read.table(paste("../DiffAFshuffled/",pop1,"___",pop2,"/",pop1,".frq",sep=""),stringsAsFactors =F, header=T)
    tmp2<-read.table(paste("../DiffAFshuffled/",pop1,"___",pop2,"/",pop2,".frq",sep=""),stringsAsFactors =F, header=T)
    merged<-merge(tmp1,tmp2,by="SNP")
    numHudson<-(merged$MAF.x-merged$MAF.y)^2 - 
      (merged$MAF.x)*(1-merged$MAF.x)/(merged$NCHROBS.x-1) -
      (merged$MAF.y)*(1-merged$MAF.y)/(merged$NCHROBS.y-1)
    denomHudson=merged$MAF.x+merged$MAF.y-2*merged$MAF.x*merged$MAF.y
    fstHudsonShuffled<-ifelse(merged$NCHROBS.x > 2 & merged$NCHROBS.y>2,numHudson/denomHudson,NA)
    fstHudsonShuffled[fstHudsonShuffled<0]<-0
    FSTbetweenModern[[paste(comp," shuffled")]]<-fstHudsonShuffled
    wilcoxTest[comp]<-wilcox.test(fstHudson,fstHudsonShuffled,paired=TRUE,alternative = "greater")$p.value
    
    diffFST[[comp]]<-fstHudson-fstHudsonShuffled
  }
}

pdf("FstDifferences.pdf")
par(mar=c(10,4,4,1)+0.1)
boxplot(FSTbetweenModern,las=2,col=rep(c("red","grey"),length(FSTbetweenModern)/2),cex=0.6,ylab="Fst",main="Fst between pairs of modern groups (red)\nas compared to Fst when shuffling group labels (grey)",axes=F)  
axis(1,at=seq(1.5,length(FSTbetweenModern)+0.5,2),names(wilcoxTest),las=2,col="grey")
axis(2)
box()
abline(v=seq(0.5,length(FSTbetweenModern)+1.5,2),lty=1)
boxplot(FSTbetweenModern,las=2,col=rep(c("red","grey"),length(FSTbetweenModern)/2),cex=0.6,outline=F,ylab="Fst",main="Fst between pairs of modern groups (red)\nas compared to Fst when shuffling group labels (grey)",axes=F)  
axis(1,at=seq(1.5,length(FSTbetweenModern)+0.5,2),names(wilcoxTest),las=1,col="grey")
axis(2)
box()
abline(v=seq(0.5,length(FSTbetweenModern)+1.5,2),lty=2)
boxplot(diffFST,las=2,ylab="Fst - Suffled Fst",main="Difference between\nFst and Fst when shuffling group labels",cex=0.6)  
#axis(3,at=c(1:length(wilcoxTest)),scientific(wilcoxTest),las=2)
boxplot(diffFST,outline = F,las=2,ylab="Fst - Suffled Fst",main="Difference between\nFst and Fst when shuffling group labels",cex=0.6)  


plot(0,0,"n",xlim=c(0,length(FSTbetweenModern)),ylim=c(0,1),ylab="Fst - Suffled Fst",main="Difference between\nFst and Fst when shuffling group labels",cex=0.6,axes=F)
for(i in c(1:length(FSTbetweenModern))){
  violin_plot(na.omit(FSTbetweenModern[[i]]),at=i,add=T,violin_width = 1,col=ifelse(i%%2==1,"red","grey"),box_width = 0.03)

}
axis(1,at=seq(1.5,length(FSTbetweenModern)+0.5,2),names(wilcoxTest),las=1,col="grey",las=2)
axis(2)
box()

par(mar=c(15,4,4,1))
barplot(t(forBarplot),las=2,cex.names = 0.6,legend.text = F,xlab = "",col=paste("grey",seq(100,10,length.out=10)),cex.lab=0.8,
        ylab="Proportion of SNPs with FST",
        cex.axis=0.6,cex.main=0.8,
        main="A\nProportion of SNPs with\na given FST between two populations")
plot(0,0,"n",axes=F,ann=F)
legend("center",pch=15,col=paste("grey",seq(100,10,length.out=10)),title="Observed FST",
       #legend = paste("|f(A in pop1) - f(A in pop2)| in [",seq(0,0.9,0.1),"-",seq(1,1,0.1),c(rep("[",8),"]"))
       legend = paste("[",seq(0,0.9,0.1),"-",seq(1,1,0.1),c(rep("[",8),"]"),sep="")
)

dev.off()



