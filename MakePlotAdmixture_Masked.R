#!/usr/bin/Rscript

library(RColorBrewer)
require(stringr)

setwd("/Users/pierrespc/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/AdmixtureMasked/")



setInd="Lab_with_Compendium_GEHmodern_GenoModern"
setSnp="1240K"
th=30000
Ref="../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv"
Lab="../../Uspallata_Annotation.tsv"
OrderIND=paste("WithoutProjection/",setInd,".",setSnp,"/BestRUNperK/OrderInd.OnK15.txt",sep="")
#OrderIND="None"

KMAX=15
folderOUT=paste("WithoutProjection/",setInd,".",setSnp,"/BestRUNperK/",sep="")
fam=read.table(paste("../DataSets/",setInd,".",setSnp,"/ForAdmixture.TH30000.GENO0.5.MAF0.01.pruned.fam",sep=""),stringsAsFactors = F,header=F)


change<-function(string,split,pos){
  return(strsplit(string,split=split)[[1]][pos])
}

fam$V6<-fam$V2
fam$V1[ grepl(":",fam$V6)]<-sapply(fam$V6[ grepl(":",fam$V6)],change,split=":",pos=1)
fam$V2[ grepl(":",fam$V6)]<-sapply(fam$V6[ grepl(":",fam$V6)],change,split=":",pos=2)
fam$V1[str_ends(fam$V1,".Mix")]<-str_replace(fam$V1[str_ends(fam$V1,".Mix")],".Mix",".Capt")
fam$V1<-str_replace(fam$V1,"_Luisi","-Luisi")
fam$V1[ fam$V1 %in% c("25deMayo-Luisi","ColoniaAurora-Luisi","SanVicente-Luisi")]<-"Posadas-Luisi"
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


orderRegions<-read.table("listRegionOrdered_withGenoModern.txt",stringsAsFactors = F,header=F)$V1

for(KNUM in c(2:15)){
  #a<-read.table(paste("chr1.1KG.PopArg.pruned.",KNUM,".Q",sep=""),stringsAsFactors=F,header=F)
  a<-read.table(paste(folderOUT,"/ForAdmixture.TH30000.GENO0.5.MAF0.01.pruned.",KNUM,".Q",sep=""),stringsAsFactors=F,header=F)
  outfile<-paste(folderOUT,"/ForAdmixture.TH30000.GENO0.5.MAF0.01.pruned.",KNUM,sep="")


  numberK=dim(a)[2]
  numberInd=dim(a)[1]

  if(KNUM==2){
    MyColors<-c("deeppink","indianred2")
  }
  
  if(KNUM==3){
    MyColors<-c("indianred2","deeppink","mediumseagreen")
  }
  
  if(KNUM==4){
    MyColors<-c("indianred2","deeppink","mediumseagreen","violet")
  }
  
  if(KNUM==5){
    
    MyColors<-c("deeppink","violet","blue1","mediumseagreen","indianred2")
  }
  
  if(KNUM==6){
    MyColors<-c("orange4","mediumseagreen","blue1","indianred2","deeppink","violet")
  }
  
  if(KNUM==7){
    MyColors<-c("indianred2","mediumseagreen","blue1","violet","red2","deeppink","orange4")
    
  }
  
  if(KNUM==8){
    MyColors<-c("violet","orange4","indianred2","darkseagreen","mediumseagreen","blue1","red2","deeppink")
    
  }
  
  
  if(KNUM==9){
    MyColors<-c("deeppink","blue1","darkseagreen","darkorange1","mediumseagreen","orange4","violet","indianred2","red2")
  }
  
  if(KNUM==10){
    MyColors<-c("goldenrod1","blue1","violet","darkseagreen","darkorange1","red2","orange4","indianred2","mediumseagreen","deeppink")
  }
  
  if(KNUM==11){
    MyColors<-c("violet","red2","indianred2","darkseagreen","goldenrod1","blue1","mediumseagreen","deeppink","palegreen1","darkorange1","orange4")
  }
  
  if(KNUM==12){
    MyColors<-c("goldenrod1","darkseagreen","palegreen1","blue1","orange3","indianred2","red2","orange4","mediumseagreen","deeppink","violet","darkorange1")
    
  }
  
  if(KNUM==13){
    MyColors<-c("violet","green4","indianred2","orange4","darkorange1","palegreen1","orange3","mediumseagreen","goldenrod1","blue1","red2","deeppink","darkseagreen")
  }
  if(KNUM==14){
    MyColors<-c("darkseagreen","palegreen1","orange4","blue1","mediumseagreen","orange3","goldenrod1","darkorange1","darkviolet","deeppink","violet","indianred2","red2","green4")
    
  }
  if(KNUM==15){
    MyColors<-c("orange3","mediumseagreen","violet","palegreen1","brown3","red2","green4","deeppink","darkseagreen","indianred2","goldenrod1","darkorange1","blue1","darkviolet","orange4")
    
  }
  if(KNUM==16){
    MyColors<-c("violet","red2","deeppink","grey1","grey40","orange4","grey60","orange3","palegreen1","darkdeeppink","darkorange1","indianred2","goldenrod1","green4","blue1","brown3")
  }
  if(KNUM==17){
    MyColors<-c("blue3","blue1","darkdeeppink","green4","grey80","violet","orange4","red2","palegreen1","brown3","orange3","goldenrod1","indianred2","grey1","deeppink","darkorange3","darkorange1")
  }
  
  if(length(unique(MyColors))!=KNUM){
    print(MyColors[duplicated(MyColors)])
    stop("pb num col")
  }
  a$Ind<-fam$V2
  a$Population<-fam$V1
  
  if(sum(!a$Population %in% Colors$Population)>0){
    print(unique(a[ ! a$Population %in% Colors$Population,"Population"]))
    stop("WTF")
  }
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
  if(sum(! orderRegions %in% a$Region)>0){
    print("order regions not represented")
    print(orderRegions[! orderRegions %in% a$Region])
  }
  if(sum(!  a$Region %in% orderRegions)>0){
    print("order regions does not include")
    print(table(a$Region[! a$Region %in% orderRegions]))
  }
  out<-c()
  for(reg in orderRegions){
    #print(reg)
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
  for(reg in orderRegions[orderRegions %in% out$Region]){
	  temp<-out[ out$Region == reg,]
	  Population=temp$PopulationLeg[1]
	
	  Population2=Population
	  #axis(3,at=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),label=Population2,cex.axis=6,las=2,tick=T)
	  
	  for(ind in unique(temp$Ind)){
		  ybottom=0
		  if(temp$PopulationLeg[ temp$Ind==ind] != Population){
			  Population=temp$PopulationLeg[ temp$Ind==ind]
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
			  xleft=xleft+separPop

			  Population2=Population
			  #axis(3,at=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),label=Population2,cex.axis=6,las=2,tick=T)
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
        #print(pop)
        iterP=iterP+1
        meanByPop[iterP,]<-apply(out[out$PopulationLeg==pop,paste("V",c(1:KNUM),sep="")],2,mean)
        rownames(meanByPop)[iterP]<-pop
      }
    }else{
      if(unique(out$PopulationLeg[out$Region==region])!=region){
        for(pop in unique(out$PopulationLeg[out$Region==region])){
          #print(pop)
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
  out<-apply(out,2,as.character)
  write.table(out,file=paste(outfile,".AncestryComponentByIndividual.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
}


out<-as.data.frame(out,stringsAsFactors = F)
if(OrderIND!="None"){
  
  separPop=1
  separ=1
  outfile<-paste("WithoutProjection/",setInd,".",setSnp,"/BestRUNperK/ForAdmixture.TH30000.GENO0.5.MAF0.01.pruned.OrderINDs",sep="")
  
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
    point=as.numeric(temp$Point[1])
    col=temp$Color[1]
    Population2=Population
    #text(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0,labels=paste(Population2,": ",reg,sep=""),cex=4.5,srt=90,adj=c(0,0.5))
    text(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0.1,labels=Population2,cex=7,srt=90,adj=c(0,0.5))
    points(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0,cex=8,pch=point,col=ifelse(point<21,col,"black"),bg=col,lwd=10)
    #text(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0,labels=reg,cex=9,srt=90,adj=c(0,0.5))
    
    for(ind in unique(temp$Ind)){
      ybottom=0
      if(temp$PopulationLeg[ temp$Ind==ind] != Population){
        Population=temp$PopulationLeg[ temp$Ind==ind]
        col=temp$Color[ temp$Ind==ind]
        point=as.numeric(temp$Point[ temp$Ind==ind])
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
        xleft=xleft+separPop
        
        Population2=Population
        #text(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0,labels=paste(Population2,": ",reg,sep=""),cex=4.5,srt=90,adj=c(0,0.5))
        text(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0.1,labels=Population2,cex=7,srt=90,adj=c(0,0.5))
        points(x=xleft+mean(c(0,sum(temp$PopulationLeg==Population))),y=0,cex=8,pch=point,col=ifelse(point<21,col,"black"),bg=col,lwd=10)
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
pdf(paste("ForAdmixture.TH30000.GENO0.5.MAF0.01.pruned.CrossValidation.pdf",sep=""))
 
a<-read.table(paste(pref,".CV.BESTruns",sep=""),stringsAsFactors = F,header=T)
a<-a[ a$K >2,]
plot(CV~K,data=a,type="l",main="Cross-Validation Score",axes=F)
axis(2,at=seq(round(min(a$CV),digits = 4),round(max(a$CV),digits = 4),by = 2e-4))
axis(1,at=a$K)
points(CV~K,data=a,pch=4)
points(a$K[ a$CV==min(a$CV)],a$CV[ a$CV==min(a$CV)],col="red",pch=1,cex=3,lwd=2)
dev.off()
