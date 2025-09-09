#!/bin/bash
require(scales)
require(stringr)
setIND="Lab_with_Compendium_GEHmodern"



for(snpTested in c("all","onlyTs")[1]){
  for(TEMPORALITY in c("","_noLH")[1]){
    if(snpTested=="all"){
      listSNPs<-c("1240K","1240K.TVs","SG","SG.TVs")
    }else{
      listSNPs<-c("1240K","SG")
    }
    
    th="50000"
    folder="~/Documents/PostDocPasteur/aDNA//2024-09-01_Uspallata_noHighCov//Analyses/qpWaves_Admixtools2_Uspallata/"
    setwd(folder)
    pdf(paste("Summary_Doublets_6Sets_4SNPpanels.",snpTested,TEMPORALITY,".pdf",sep=""),height=ifelse(TEMPORALITY=="",15,10),width=8)
    par(mar=c(5,15,4,11))
    
    target="Uspallata"
    info<-read.table("../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv",stringsAsFactors = F,header=T,sep="\t")
    listRightPops<-read.table("RightPops_Sets.tsv",stringsAsFactors=F,header=T)
    
    listSETS<-names(listRightPops[-c(1:6)])
    print(nrow(listRightPops))
    info<-rbind(info[,c("Population","Region","Point","Color")],cbind("Population"=c("Late","Early"),"Region"="Uspallata","Point"=21,"Color"="goldenrod1"))
    info$Point=as.numeric(info$Point)
    listRightPops<-merge(listRightPops,info[,c("Population","Point","Color")],by.x="Pop",by.y="Population")
    
    listOrder<-read.table("listRegionOrdered.txt")
    
    
    tmp<-unique(listRightPops$Region[! listRightPops$Region %in% listOrder$V1])
    if(length(tmp)>1){
      print("not in order")
      print(tmp)
      stop()
    }
    
    listRightPopsBU<-listRightPops
    listRightPops<-c()
    for(reg in listOrder$V1[c(length(listOrder$V1):1)]){
      tmp<-listRightPopsBU[listRightPopsBU$Region ==reg,]
      listRightPops<-rbind(listRightPops,tmp[order(tmp$Pop),])
    }
    print(nrow(listRightPops))
    
    for(setSNP in listSNPs){
      a<-read.table(paste("Summary_Doublets_",target,"_",length(listSETS),"Sets_",setSNP,".TH",th,".tsv",sep=""),sep="\t",header=T,stringsAsFactors = F)
      names(a)[grepl("Rank",names(a))]<-paste(names(a)[grepl("Rank",names(a))],setSNP,sep="_")
      if(setSNP=="1240K"){
        merged<-a
      }else{
        merged<-merge(merged,a,by=c("LeftPops","Region","Point","Color"),all=T)
      }
    }
    
    
    if(TEMPORALITY!=""){
      merged<-merged[ ! (str_ends(merged$Region,"_LH") | str_ends(merged$Region,"_Unknown")),]
    }
    merged
    out<-c()
    for(reg in listOrder$V1[listOrder$V1 %in% merged$Region]){
      tmp<-merged[merged$Region ==reg,]
      out<-rbind(out,tmp[order(tmp$LeftPops),])
    }
    
    
    plot(0,0,"n",xlim=c(0,(length(listSNPs)*length(listSETS)+1)),ylim=c(0,(nrow(listRightPops)+nrow(out)+10)),
         axes=F,ann=F)
    
    #title(main=paste("Testing Rank 0 rejection with 2 left pops:\n1 South American Group (y-axis) & ",target,"\nUsing ",length(listSETS)," sets for right populations\n and 4 SNP panels",sep=""),cex.main=0.9)
    axis(3,at=seq(length(listSNPs)/2+.5,length(listSNPs)*length(listSETS)+0.5,length(listSNPs)),paste("Set",c(1:length(listSETS))),tick = F,line=ifelse(TEMPORALITY=="",-3,-2),cex.axis=0.7)
    
    y=0
    
    for(i in c(1:nrow(out))){
      y=y+1
      snp=0
      for(setSNP in listSNPs){
        snp=snp+1
        for(set in c(1:length(listSETS))){
          rect(xleft = (set-1)*length(listSNPs)+snp-0.5,xright=(set-1)*length(listSNPs)+snp+0.5,
               ybottom = y-0.5,ytop=y+0.5,
               border="black",
               col=ifelse(is.na(out[i,paste("Rank1_",listSETS[set],"_",setSNP,sep="")]),"grey",
                          ifelse(out[i,paste("Rank1_",listSETS[set],"_",setSNP,sep="")]<0.001,"red",
                                 ifelse(out[i,paste("Rank1_",listSETS[set],"_",setSNP,sep="")]<0.005,"orange",
                                        ifelse(out[i,paste("Rank1_",listSETS[set],"_",setSNP,sep="")]<0.01,"yellow3",
                                               ifelse(out[i,paste("Rank1_",listSETS[set],"_",setSNP,sep="")]<0.05,"yellow1",
                                                      "white"))))))
          
        }
      }
      points(x=-0.4,y=y,bg=out$Color[i],pch=out$Point[i],col=ifelse(out$Point[i]<21,out$Color[i],"black"),cex=0.8)
      #axis(2,at=y,labels = paste(out$Region[i],": ",out$LeftPops[i],sep=""),line=0,las=2,tick=F,cex.axis=0.6)
      axis(2,at=y,labels = out$LeftPops[i],line=0,las=2,tick=F,cex.axis=0.6)
    }
    
    
    
    ##### ADD REF TO REGIONS
    segments(x0 = -2,x1=length(listSNPs)*length(listSETS)+0.5,y0=0.5,lwd=3)
    prevY=1
    for(i in c(2:nrow(out))){
      if(out$Region[i]!=out$Region[i-1]){
        segments(x0 = -2,x1=length(listSNPs)*length(listSETS)+0.5,y0=i-0.5,lwd=3)
        axis(4,at=mean(c(i,prevY))-0.5,labels = out$Region[i-1],las=2,cex.axis=0.6,tick = F,line=-1)
        prevY=i
      }
    }
    segments(x0 = -2,x1=length(listSNPs)*length(listSETS)+0.5,y0=i+0.5,lwd=3)
    axis(4,at=mean(c(i,prevY)),labels = out$Region[i],las=2,cex.axis=0.6,tick = F,line=-1)
    
    
    axis(1,at=seq(1,length(listSNPs)*length(listSETS),1),las=2,labels = rep(listSNPs,length(listSETS)),cex.axis=0.6,line = ifelse(TEMPORALITY=="",-2,-1.5))
    
    #segments(x0 = seq(0.5,24.5,4),y0=ybreak,y1=y+0.5,lwd=3)
    
    
    ybreak=y+0.5
    segments(x0 = seq(0.5,length(listSNPs)*length(listSETS)+0.5,length(listSNPs)),y0=0.5,y1=ybreak,lwd=3)
    #abline(v=seq(0.5,24.5,4),lwd=2)
    y=y+2
    text(x=seq(length(listSNPs)/2+0.5,length(listSNPs)*length(listSETS)+0.5,length(listSNPs)),y=y,labels = paste("Set",c(1:length(listSETS))),cex=0.7)
    
    y=y+4
    
    
    
    rect(xleft=c(0,4.25,8.5,12.75,17)/length(listSNPs)*4,xright=c(4,8.25,12.5,16.75,21)/length(listSNPs)*4,
         ybottom=rep(y-0.75,5),ytop=rep(y+0.75,5),
         col=c("red","orange","yellow3","yellow1","grey"))
    text(x=c(2,6.25,10.5,14.75,19)/length(listSNPs)*4,y=rep(y,3),labels = c(paste("P<",c(0.001,0.005,0.01,0.05),sep=""),"NA"),cex=0.7)
    
    rect(xleft=0,xright=10/length(listSNPs)*4,
         ybottom=y+1.5,ytop=y+3,
         col="lightblue")
    text(x=5/length(listSNPs)*4,y=y+2.25,labels = "included as right pop",cex=0.7)
    
    
    
    y=y+4
    ybreak=y+0.5
    for(i in c(1:nrow(listRightPops))){
      y=y+1
      for(set in c(1:length(listSETS))){
        rect(xleft = (set-1)*length(listSNPs)+0.5,xright=(set)*length(listSNPs)+0.5,
             ybottom = y-0.5,ytop=y+0.5,
             border="black",
             col=ifelse(listRightPops[i,listSETS[set]],"lightblue","white"))
        
      }
      points(x=-0.4,y=y,bg=listRightPops$Color[i],pch=listRightPops$Point[i],col=ifelse(listRightPops$Point[i]<21,listRightPops$Color[i],"black"),cex=0.8)
      #axis(2,at=y,labels = paste(listRightPops$Region[i],": ",listRightPops$Pop[i],sep=""),line=0,las=2,tick=F,cex.axis=0.6)
      axis(2,at=y,labels = listRightPops$Pop[i],line=0,las=2,tick=F,cex.axis=0.6)
    }
    
    
    ##### ADD REF TO REGIONS
    segments(x0 = -2,x1=length(listSETS)*length(listSNPs)+0.5,y0=ybreak,lwd=3)
    prevY=ybreak
    y=ybreak
    for(i in c(2:nrow(listRightPops))){
      y=y+1
      print(c(y,ybreak,i))
      if(listRightPops$Region[i]!=listRightPops$Region[i-1]){
        segments(x0 = -2,x1=length(listSETS)*length(listSNPs)+0.5,y0=y,lwd=3)
        axis(4,at=mean(c(y,prevY)),labels = listRightPops$Region[i-1],las=2,cex.axis=0.6,tick = F,line=-1)
        prevY=y
      }
    }
    segments(x0 = -2,x1=length(listSETS)*length(listSNPs)+0.5,y0=y+1,lwd=3)
    axis(4,at=mean(c(y,prevY))+0.5,labels = listRightPops$Region[i],las=2,cex.axis=0.6,tick = F,line=-1)
    
    
    segments(x0 = seq(0.5,length(listSETS)*length(listSNPs)+0.5,4),y0=y+1,y1=ybreak,lwd=3)
    dev.off()
    
  }
}
stop()
#####onkly SG

out<-out[ str_ends(out$LeftPops,".SG"),]


plot(0,0,"n",xlim=c(0,(4*length(listSETS)+1)),ylim=c(0,(nrow(listRightPops)+nrow(out)+10)),
     axes=F,ann=F)
title(main=paste("Testing Rank 0 rejection with doublets:\n1 South American Group (y-axis) & ",target,"\nUsing ",length(listSETS)," sets for right populations\n and 4 SNP panels",sep=""),cex.main=0.9)
axis(3,at=seq(2.5,4*length(listSETS),4),paste("Set",c(1:length(listSETS))),tick = F,line=-3,cex.axis=0.9)
y=0
for(i in c(1:nrow(out))){
  y=y+1
  snp=0
  for(setSNP in c("1240K","1240K.TVs","SG","SG.TVs")){
    snp=snp+1
    for(set in c(1:length(listSETS))){
      rect(xleft = (set-1)*4+snp-0.5,xright=(set-1)*4+snp+0.5,
           ybottom = y-0.5,ytop=y+0.5,
           border="black",
           col=ifelse(is.na(out[i,paste("Rank1_",listSETS[set],"_",setSNP,sep="")]),"grey",
                      ifelse(out[i,paste("Rank1_",listSETS[set],"_",setSNP,sep="")]<0.01,"red",
                             ifelse(out[i,paste("Rank1_",listSETS[set],"_",setSNP,sep="")]<0.05,"orange",
                                    ifelse(out[i,paste("Rank1_",listSETS[set],"_",setSNP,sep="")]<0.1,"yellow2","white")))))
      
    }
  }
  points(x=-0.4,y=y,bg=out$Color[i],pch=out$Point[i],col=ifelse(out$Point[i]<21,out$Color[i],"black"))
  axis(2,at=y,labels = paste(out$Region[i],": ",out$LeftPops[i],sep=""),line=0,las=2,tick=F,cex.axis=0.6)
}

axis(1,at=seq(1,24,1),las=2,labels = rep(c("1240K","1240K.TVs","SG","SG.TVs"),6), cex.axis=0.8)


y=y+2
text(x=seq(2.5,4*length(listSETS),4),y=y,labels = paste("Set",c(1:length(listSETS))),cex=0.9)

y=y+4

abline(v=seq(0.5,24.5,4),lwd=2)

rect(xleft=c(0,4.25,8.5,12.75),xright=c(4,8.25,12.5,16.75),
     ybottom=rep(y-0.75,4),ytop=rep(y+0.75,4),
     col=c("red","orange","yellow2","grey"))
text(x=c(2,6.25,10.5,14.75),y=rep(y,3),labels = c(paste("P <",c(0.01,0.05,0.1)),"NA"),cex=0.8)

y=y+4
for(i in c(1:nrow(listRightPops))){
  y=y+1
  for(set in c(1:length(listSETS))){
    rect(xleft = (set-1)*4+0.5,xright=(set)*4+0.5,
         ybottom = y-0.5,ytop=y+0.5,
         border="black",
         col=ifelse(listRightPops[i,listSETS[set]],"black","white"))
    
  }
  points(x=-0.4,y=y,bg=listRightPops$Color[i],pch=listRightPops$Point[i],col=ifelse(listRightPops$Point[i]<21,listRightPops$Color[i],"black"))
  axis(2,at=y,labels = paste(listRightPops$Region[i],": ",listRightPops$Pop[i],sep=""),line=0,las=2,tick=F,cex.axis=0.6)
  
}
}
dev.off()


