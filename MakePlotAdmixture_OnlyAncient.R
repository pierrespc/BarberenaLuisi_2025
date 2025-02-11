#!/usr/bin/Rscript

library(RColorBrewer)
require(stringr)

setwd("/Users/pierrespc/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/AdmixtureOnlyAncient/")



setInd="Lab_with_Compendium"
setSnp="1240K"
th=30000
Ref="../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv"
Lab="../../Uspallata_Annotation.tsv"
#OrderIND=paste("WithoutProjection/",setInd,".",setSnp,"/BestRUNperK/orderInds.txt",sep="")
OrderIND=paste("WithoutProjection/",setInd,".",setSnp,"/BestRUNperK/orderInds.ForK15.txt",sep="")
#OrderIND="None"
KMAX=15
pref=paste("WithoutProjection/",setInd,".",setSnp,"/BestRUNperK/ForAdmixture.TH",th,".MAF0.01.pruned",sep="")

fam<-read.table(paste("../DataSets/",setInd,".",setSnp,"/ForAdmixture.TH",th,".MAF0.01.pruned.fam",sep=""),stringsAsFactors = F,header=F)


ColorsRef<-read.table(Ref,stringsAsFactors = F,header=T,sep="\t",comment.char = "@")
ColorsRef$PopulationLeg<-str_remove(str_remove(str_remove(ColorsRef$Population,".Capt"),".Mix"),".SG")
ColorsRef$cex=0.7

ColorsLab<-read.table(Lab,stringsAsFactors = F,header=T,sep="\t",comment.char = "@",quote="\"",fileEncoding="latin1")
ColorsLab<-ColorsLab[,c("MainRegion","MainRegion","Individual","Color","Point")]
names(ColorsLab)<-c("Region","PopulationLeg","Population","Color","Point")
ColorsLab$Set<-"PresentStudy"
ColorsLab$cex=1.2
checkDUP<-unique(ColorsRef[,c("Point","Color")])
if(sum(duplicated(checkDUP))){
  stop("duplicated Color/Point in Ref")
}


Colors<-rbind(unique(ColorsRef[,names(ColorsLab)]),ColorsLab)


orderRegions<-read.table("listRegionOrdered.txt",stringsAsFactors = F,header=F)$V1

for(KNUM in c(2:17)){
  #a<-read.table(paste("chr1.1KG.PopArg.pruned.",KNUM,".Q",sep=""),stringsAsFactors=F,header=F)
  getFile=paste(pref,".",KNUM,".Q",sep="")
  print(getFile)
	a<-read.table(getFile,stringsAsFactors=F,header=F)
#	outfile<-paste("chr1.1KG.PopArg.pruned.K",KNUM,".pdf",sep="")
	outfile<-paste(pref,".",KNUM,sep="")


  numberK=dim(a)[2]
  numberInd=dim(a)[1]

  if(KNUM==2){
    MyColors<-c("blue1","indianred2")
  }
  
  if(KNUM==3){
    MyColors<-c("indianred2","goldenrod1","blue1")
  }
  
  if(KNUM==4){
    MyColors<-c("goldenrod1","blue1","orange4","indianred2")
  }
  
  if(KNUM==5){
    MyColors<-c("orange4","blue1","indianred2","goldenrod1","green4")
  }
  
  if(KNUM==6){
    MyColors<-c("goldenrod1","green4","blue1","indianred2","red2","orange4")
    #MyColors<-c("goldenrod1","green4","blue1","grey1" ,"goldenrod1","blue1")
  }
  
  if(KNUM==7){
    MyColors<-c("red2","deeppink","indianred2","orange4","blue1","green4","goldenrod1")
  }
  
  if(KNUM==8){
    MyColors<-c("goldenrod1","green4","indianred2","red2","orange4","violet","deeppink","blue1")
    
  }
  
  
  if(KNUM==9){
    MyColors<-c("goldenrod1","blue1","green4","indianred2","deeppink","orange4","grey1","red2","grey80")
    MyColors<-c("deeppink","violet","goldenrod1","orange4","red2","green4","indianred2","blue1","darkviolet")
  }
  
  if(KNUM==10){
    MyColors<-c("goldenrod1","indianred2","darkviolet","darkorange1","violet","deeppink","red2","blue1","orange4","green4")
  }
  
  if(KNUM==11){
    
    MyColors<-c("palegreen1","violet","goldenrod1","deeppink","blue1","indianred2","orange4","green4","red2","darkviolet","darkorange1")
    
  }
  
  if(KNUM==12){
    MyColors<-c("violetred","orange4","red2","violet","darkviolet","indianred2","palegreen1","darkorange1","green4","blue1","goldenrod1","deeppink")
    
  }
  
  if(KNUM==13){
    MyColors<-c("orange4","darkviolet","orange3","blue4","indianred2","green4","darkorange1","goldenrod1","red2","violet","deeppink","blue1","palegreen1")
  }
  if(KNUM==14){
    MyColors<-c("green4","orange4","goldenrod1","blue1","lightblue4","blue4","deeppink","darkorange1","red2","violetred","indianred2","violet","orange3","palegreen1")
  }
  if(KNUM==15){
    MyColors<-c("darkorange1","indianred2","blue1","darkviolet","green4","orange4","violet","orange3","red2","blue4","deeppink","palegreen3","goldenrod1","violetred","palegreen1")
  }
  if(KNUM==16){
    MyColors<-c("deeppink","red2","violet","grey1","grey40","orange4","grey60","orange3","palegreen1","darkviolet","darkorange1","indianred2","goldenrod1","green4","blue1","brown3")
  }
  if(KNUM==17){
    MyColors<-c("blue3","blue1","darkviolet","green4","grey80","deeppink","orange4","red2","palegreen1","brown3","orange3","goldenrod1","indianred2","grey1","violet","darkorange3","darkorange1")
  }
  
  if(length(unique(MyColors))!=KNUM){
    print(MyColors[duplicated(MyColors)])
    stop("pb num col")
  }
  a$Ind<-fam$V2
  a$Population<-fam$V1
  
  a<-merge(a,Colors,by="Population")

  

  #numPops<-max(table(unique(a[,c("Population","Region")])$Region))
  #numRegions<-1
  #numInds<-max(table(unique(a[,c("Ind","Region")])$Region))
  
  numInds<-nrow(a)
  numPops<-length(unique(a$PopulationLeg))
  numRegions<-length(unique(a$Region))
  
  if(dim(a)[1] != numberInd){
	  stop("your regions file and admixture output do not coincide: do not have the same number of Pops")
  }

  out<-c()
  for(reg in orderRegions){
    print(reg)
    temp<-a[ a$Region == reg,]
    if(dim(temp)[1]==0){
      print(paste(reg,"skipped"))
      
      next
    }	
	  meanOverRegion<-vector(length=numberK)
	  names(meanOverRegion)<-paste("V",c(1:numberK),sep="")
	
	
	  meanOverPop<-data.frame(matrix(NA,length(unique(temp$PopulationLeg)),numberK))
	  names(meanOverPop)<-paste("V",c(1:numberK),sep="")
	  rownames(meanOverPop)<-unique(temp$PopulationLeg)
	
	  for(K in c(1:numberK)){
		  meanOverRegion[paste("V",K,sep="")]<-mean(temp[,paste("V",K,sep="")])
		  for(pop in unique(temp$PopulationLeg)){
			  meanOverPop[pop,paste("V",K,sep="")]<-mean(temp[temp$PopulationLeg==pop,paste("V",K,sep="")])
		  }
	  }
	
	  meanOverRegion<-meanOverRegion[order(as.numeric(meanOverRegion),decreasing=T)]
	  meanOverPop<-meanOverPop[,names(meanOverRegion)]
	  meanOverPop<-meanOverPop[order(meanOverPop[,1],decreasing=T),]
	  #print(head(meanOverPop))
	  for(pop in rownames(meanOverPop)){
		  temp2<-a[ a$PopulationLeg == pop,]
		  temp2<-temp2[order(temp2[,names(meanOverRegion)[1]],decreasing=T),]
		  out<-rbind(out,temp2)
	  }
  }

  if(OrderIND !="None"){
    orderInd<-read.table(OrderIND,stringsAsFactors=F,header=F)
    rownames(out)<-out$Ind
    out<-out[ orderInd$V1,]
    
  }
  

  #separ<-ceiling(numberInd/500)
  #separPop<-ceiling(numberInd/1000)
  separPop=1
  separ=1
  pdf(paste(outfile,".pdf",sep=""), width=400,height=20)
  #par(mar=c(1, 2, 40, 2) + 0.1)
  par(mar=c(0,0,0,0))

  plot(0,0,"n",ylim=c(0,1),xlim=c(0,(numInds+separPop*(numPops-numRegions)+numRegions*(separ)))*1.1,main=paste("K =", numberK),ann=F,axes=F)
  xleft=0
  atPop<-c()
  dimPrevRegion<- 0
  atRegion<-c()
  for(reg in orderRegions){
	  temp<-out[ out$Region == reg,]
	  Population=temp$PopulationLeg[1]
	
	  Population2=Population
	  axis(3,at=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),label=Population2,cex.axis=6,las=2,tick=T)
	  
	  for(ind in unique(temp$Ind)){
		  ybottom=0
		  if(temp$PopulationLeg[ temp$Ind==ind] != Population){
			  Population=temp$PopulationLeg[ temp$Ind==ind]
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
			  xleft=xleft+separPop

			  Population2=Population
			  axis(3,at=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),label=Population2,cex.axis=6,las=2,tick=T)
			  dimPrevRegion=dimPrevRegion+separPop
		  }
		  for(k in c(1:numberK)){
			  ytop=ybottom+temp[ temp$Ind==ind,paste("V",k,sep="")]
			  rect(xleft=xleft,xright=xleft+1,ybottom=ybottom,ytop=ytop,col=MyColors[k],border=NA)
			  ybottom=ytop
		  }
		  xleft=xleft+1
	  }
	  xleft=xleft+separ
	  atRegion[reg]<-mean(c(0,sum(temp$Region==reg)))+dimPrevRegion
	  dimPrevRegion=atRegion[reg]+mean(c(0,sum(temp$Region==reg)))+separ
  }

  
  forLeg<-unique(out[,c("PopulationLeg","Color")])
  #plot(0,0,"n",pch=0,ann=F,axes=F)
  #legend("topleft",legend=paste(c(1:dim(forLeg)[1]),": ",forLeg$Population,sep=""),ncol=10,cex=8,text.col=forLeg$Color)

  
  dev.off()
  
  meanByPop=data.frame(matrix(NA,0,KNUM))
  names(meanByPop)=paste(c(1:KNUM),sep="")
  iterP=0
  for(region in unique(out$Region)){
    iterP=iterP+1
    meanByPop[iterP,]<-apply(out[out$Region==region,paste("V",c(1:KNUM),sep="")],2,mean)
    rownames(meanByPop)[iterP]<-region
    if(length(unique(out$PopulationLeg[out$Region==region]))>1){
      for(pop in unique(out$PopulationLeg[out$Region==region])){
        print(pop)
        iterP=iterP+1
        meanByPop[iterP,]<-apply(out[out$PopulationLeg==pop,paste("V",c(1:KNUM),sep="")],2,mean)
        rownames(meanByPop)[iterP]<-pop
      }
    }else{
      if(unique(out$PopulationLeg[out$Region==region])!=region){
        for(pop in unique(out$PopulationLeg[out$Region==region])){
          print(pop)
          iterP=iterP+1
          meanByPop[iterP,]<-apply(out[out$PopulationLeg==pop,paste("V",c(1:KNUM),sep="")],2,mean)
          rownames(meanByPop)[iterP]<-pop
        }
      }
    }
  }
  #write.table(meanByPop,paste("chr1.1KG.PopArg.pruned.",KNUM,".MeanByGroup.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
  names(meanByPop)<-MyColors
  write.table(meanByPop,paste(outfile,".MeanByGroup.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
  #  write.table(out,paste("chr1.1KG.PopArg.pruned.",KNUM,".AncestryComponentByIndividual.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)  
  names(out)[names(out) %in% paste("V",c(1:KNUM),sep="")]<-MyColors
  write.table(out,paste(outfile,".AncestryComponentByIndividual.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
}


if(OrderIND!="None"){
  
  separPop=1
  separ=1
  outfile<-paste(pref,".OrderINDs",sep="")
  
  pdf(paste(outfile,".pdf",sep=""), width=400,height=30)
  #par(mar=c(1, 2, 40, 2) + 0.1)
  par(mar=c(0,0,0,0))
  
  plot(0,0,"n",ylim=c(0,1),xlim=c(0,(numInds+separPop*(numPops-numRegions)+numRegions*(separ)))*1.1,main=paste("K =", numberK),ann=F,axes=F)
  xleft=0
  atPop<-c()
  dimPrevRegion<- 0
  atRegion<-c()
  for(reg in orderRegions){
    temp<-out[ out$Region == reg,]
    Population=temp$PopulationLeg[1]
    point=temp$Point[1]
    col=temp$Color[1]
    Population2=Population
    #text(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0,labels=paste(Population2,": ",reg,sep=""),cex=4.5,srt=90,adj=c(0,0.5))
    text(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0.1,labels=Population2,cex=10,srt=90,adj=c(0,0.5))
    points(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0,cex=25,pch=point,col=ifelse(point<21,col,"black"),bg=col,lwd=15)
    #text(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0,labels=reg,cex=9,srt=90,adj=c(0,0.5))
    
    for(ind in unique(temp$Ind)){
      ybottom=0
      if(temp$PopulationLeg[ temp$Ind==ind] != Population){
        Population=temp$PopulationLeg[ temp$Ind==ind]
        col=temp$Color[ temp$Ind==ind]
        point=temp$Point[ temp$Ind==ind]
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
        xleft=xleft+separPop
        
        Population2=Population
        #text(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0,labels=paste(Population2,": ",reg,sep=""),cex=4.5,srt=90,adj=c(0,0.5))
        text(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0.1,labels=Population2,cex=10,srt=90,adj=c(0,0.5))
        points(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0,cex=25,pch=point,col=ifelse(point<21,col,"black"),bg=col,lwd=15)
        
        #text(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0,labels=reg,cex=9,srt=90,adj=c(0,0.5))
        dimPrevRegion=dimPrevRegion+separPop
      }
      xleft=xleft+1
    }
    xleft=xleft+separ
    atRegion[reg]<-mean(c(0,sum(temp$Region==reg)))+dimPrevRegion
    dimPrevRegion=atRegion[reg]+mean(c(0,sum(temp$Region==reg)))+separ
  }
  
  
  dev.off()
  
}
 pdf(paste(pref,".CrossValidation.pdf",sep=""))
 
a<-read.table(paste(pref,".CV.BESTruns",sep=""),stringsAsFactors = F,header=T)
a<-a[ a$K >2,]
plot(CV~K,data=a,type="l",main="Cross-Validation Score",axes=F)
axis(2,at=seq(round(min(a$CV),digits = 4),round(max(a$CV),digits = 4),by = 2e-4))
axis(1,at=a$K)
points(CV~K,data=a,pch=4)
points(a$K[ a$CV==min(a$CV)],a$CV[ a$CV==min(a$CV)],col="red",pch=1,cex=3,lwd=2)
dev.off()
