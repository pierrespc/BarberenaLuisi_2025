#!/bin/Rscript
require(stringr)


setwd("~/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F//")


annotfile<-"../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv"
annotstudy<-"../../Uspallata_Annotation.tsv"
annot<-read.table(annotfile,stringsAsFactors = F,header=T,sep="\t",comment.char="@")
annot<-unique(annot[,c("Region","Population","Color","Point")])
annot<-annot[ annot$Population !="Peru_SoroMikayaPatjxa_6800BP.SG",]
study<-read.table(annotstudy,stringsAsFactors = F,header=T,sep="\t",comment.char="@")[,c("MainRegion","MainRegion","Color","Point")]

names(study)[c(1,2)]<-c("Region","Population")
study$Point[study$Color=="darkorange4"]<-21
study<-unique(study)

annot<-rbind(annot,study)

listRegOrdered<-read.table("../F4_UspallataOneGroup_toPairsModern/listRegionOrdered.txt",stringsAsFactors = F,header=F)$V1

fam<-read.table("../../Analyses/DataSets/Lab_with_Compendium_GEHmodern_GenoModern.1240K/ForAdmixture.TH30000.GENO0.5.MAF0.01.pruned.fam",stringsAsFactors = F,header=F)
nInd<-table(fam$V1)

Kept<-names(nInd)[ nInd>0]
th="30000"
for( setsnp in c("1240K")){
  f4file<-paste("F4_UspallataOneGroup_toPairsModern/F4.UspallataOneGroup_toPairsModern.finalSet.Lab_with_Compendium_GEHmodern_GenoModern.",setsnp,".TH",th,".OUT",sep="")
  
  for(focus in c("Uspallata")){
    outfile<-str_replace(f4file,"OUT",paste(focus,".pdf",sep=""))
    pdf(outfile,width=12,height=12)
    par(mfrow=c(1,1),mar=c(4,4,14,6)+0.1)
    print(c(setsnp,focus))
    f4<-read.table(f4file,stringsAsFactors = F,header=T,sep="\t")
    f4<-f4[ f4$Pop2_X==focus,]
    f4bis=f4
    prout=f4bis$Pop3_Y
    f4bis$Pop3_Y<-f4bis$Pop4_Z
    f4bis$Pop4_Z<-prout
    f4bis$Dstat=-f4bis$Dstat
    f4bis$Z=-f4bis$Z
    f4<-rbind(f4,f4bis)
    print(paste("dim before merging annot",nrow(f4)))
    print(unique(f4$Pop3_Y[! f4$Pop3_Y %in% annot$Population]))
    f4<-merge(f4,annot,by.x="Pop3_Y",by.y="Population")
    print(paste("dim after merging annot Pop3_Y",nrow(f4)))
    print(unique(f4$Pop4_Z[! f4$Pop4_Z %in% annot$Population]))
    f4<-merge(f4,annot,by.x="Pop4_Z",by.y="Population")
    print(paste("dim after merging annot Pop4_Z",nrow(f4)))
    
    
    #listRegOrdered<-read.table("listRegionOrdered.txt",stringsAsFactors=F,header=F)$V1
    listRefReg<-listRegOrdered[  listRegOrdered %in% f4$Region.x & !(grepl("Unmasked",listRegOrdered))]
    check1<-unique(f4$Region.x[!f4$Region.x %in% listRegOrdered])
    if(length(check1)>1){
      warning(paste("those regions are missing in your orderRegion list: ",paste(check1,collapse=" | ")))
    }
    check2<-unique(f4$Region.y[!f4$Region.y %in% listRegOrdered])
    if(length(check2)>1){
      warning(paste("those regions are missing in your orderRegion list: ",paste(check2,collapse=" | ")))
    }
    
    
    
    if(length(listRefReg) != length(listRegOrdered)){
      print(paste("watch out some regions in your order not represented in F4:",paste(listRegOrdered[ ! listRegOrdered %in% listRefReg],collapse=" | ")))
      
    }
    
    f4<-f4[ f4$Region.x %in% listRefReg & f4$Region.y %in% listRefReg,]
  
    f4<-f4[ f4$Pop3_Y %in%Kept & f4$Pop4_Z %in%Kept,]
    ###order listRefReg in decreasing median Zvalue
    #medianZ<-c()
    meanZ<-c()
    for(reg in listRefReg){
      tmp<-f4[ f4$Region.x==reg & f4$Region.y!=reg ,]
      #medianZ[reg]<-median(tmp$Z)
      meanZ[reg]<-mean(tmp$Z)
    }
    
    #listRefReg<-names(sort(medianZ,decreasing = T))
    listRefReg<-names(sort(meanZ,decreasing = T))
    
    #order f4
    f4bis<-c()
    for(reg in listRefReg){
      tmp<-f4[ f4$Region.x==reg,]
      f4bis<-rbind(f4bis,tmp[order(tmp$Pop3_Y),])
    }
    
    print(nrow(f4bis))
    f4<-f4bis
    remove(f4bis)
    listPOPS<-unique(f4$Pop3_Y)
    plot(0,0,"n",xlim=c(-4,(length(listPOPS)+0.5)),ylim=c(-4,(length(listPOPS)+0.5)),ann=F,axes=F)
    title(main=paste(focus," closer to Pop1 or Pop2?"),line=12)
    x=0
    reg1="HA"
    for(pop1 in listPOPS){
      x=x+1
      y=0
      #print(x)
      for(pop2 in listPOPS){
        y=y+1
        
        if(pop1==pop2){
          next
        }
        tmp<-f4[ f4$Pop3_Y==pop1 & f4$Pop4_Z==pop2,]
        rect(xleft = x-0.5,xright=x+0.5,ybottom=y-0.5,ytop=y+0.5,
             col=ifelse(tmp$Z<= -5,"blue4",
                        ifelse(tmp$Z> -5 & tmp$Z <= -3,"blue2",
                               ifelse(tmp$Z> -3 & tmp$Z <= -2,"lightblue3",
                                      ifelse(tmp$Z> -2 & tmp$Z <= 0,"lightblue1",
                                             ifelse(tmp$Z> 0 & tmp$Z <= 2,"khaki1",
                                                    ifelse(tmp$Z> 2 & tmp$Z <= 3,"khaki3",
                                                           ifelse(tmp$Z> 3 & tmp$Z <= 5,"brown2",
                                                                  ifelse(tmp$Z> 5,"brown4",stop(tmp$Z))))))))),
             border=NA
        )
        
      }
      points(x=x,y=-1,
             pch=tmp$Point.x,
             col=ifelse(tmp$Point.x<21,tmp$Color.x,"black"),
             bg=tmp$Color.x,
             lwd=2,cex=1)
      points(y=x,x=-1,
             pch=tmp$Point.x,
             col=ifelse(tmp$Point.x<21,tmp$Color.x,"black"),
             bg=tmp$Color.x,
             lwd=2,cex=1)
      
      print(c(x,tmp$Pop3_Y))
      axis(1,at=x,tick = F,labels = tmp$Pop3_Y,las=2,cex.axis=0.6,line = -3.25)
      axis(2,at=x,tick = F,labels = tmp$Pop3_Y,las=2,cex.axis=0.6,line = -4)
      #axis(1,x=x,las=2,tick=F,)
      #if(!str_ends(pop1,".Capt")){
      #  points(x=-3,y=x,
      #         pch=16,
      #         col="black",
      #         lwd=2,cex=1)
      #  points(x=x,y=-3,
      #         pch=16,
      #         col="black",
      #         lwd=2,cex=1)
      #}
    }
    
    reg1="HA"
    x=0
    for(pop1 in listPOPS){
      x=x+1
      newreg1=unique(f4$Region.x[ f4$Pop3_Y==pop1])
      if(newreg1 != reg1){
        #abline(v=x-0.5,h=x-0.5)
        segments(x0 = -3,x1=length(listPOPS)+0.5,y0 = x-0.5,y1=x-0.5)
        segments(y0 = -3,y1=length(listPOPS)+0.5,x0 = x-0.5,x1=x-0.5)
        reg1=newreg1
        
      }
      
    }
    segments(x0 = -3,x1=length(listPOPS)+0.5,y0 = x+0.5,y1=x+0.5)
    segments(y0 = -3,y1=length(listPOPS)+0.5,x0 = x+0.5,x1=x+0.5)
    
    
    xleft=1
    for(reg in listRefReg){
      xright=xleft+length(unique(f4$Pop3_Y[ f4$Region.x==reg]))-1
      at=mean(c(xleft,xright))
      axis(3,at=at,labels = str_replace_all(reg,"_"," "),cex.axis=0.6,las=2,tick=F,line = -2)
      axis(4,at=at,labels = str_replace_all(reg,"_"," "),cex.axis=0.6,las=2,tick = F,line = -2)
      xleft=xright+1
    }
    write.table(f4,str_replace(outfile,".pdf",".tsv"),col.names=T,row.names=F,sep="\t",quote=F) 
    
  }
  plot(0,0,"n",axes=F,ann=F)
  legend("center",
         pch=c(rep(15,8),16),
         col=c("blue4","blue2","lightblue3","lightblue1","khaki1","khaki3","brown2","brown4","black"),
         legend = c("Z <= -5",
                    "-5 < Z <=-3",
                    "-3 < Z <= -2",
                    "-2 < Z <= 0",
                    "0 < Z <= 2",
                    "2 < Z <= 3",
                    "3 < Z <= 5",
                    "Z > 5",
                    "Shotgun sequencing"))
  dev.off()
}

