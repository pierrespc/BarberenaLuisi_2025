#!/bin/Rscript

require(stringr)
#library(reshape2)
#library(ggrepel)
library(ggplot2)
library(scales)
#require(viridis)

setwd("~/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F/")

###########################################################################
###########################################################################
###########################################################################
###########################################################################

plotMDS<-function(mdsFun,Maintitle,percentage,pdfName,Doggplot=T){
  #print(head(mdsFun))
  #print(sum(str_starts(names(mdsFun),"X")))
  #pdf(pdfName)
  #svg(pdfName)
  if(! Doggplot ){
    for(i in seq(1,sum(str_starts(names(mdsFun),"X")),2)){
      svg(paste(pdfName,"_Dim",i,".Dim",i+1,".svg",sep=""))
      plot(mdsFun[,paste("X",i,sep="")],mdsFun[,paste("X",i+1,sep="")],"n",
           main=paste(Maintitle,"\nDim. ",i+1," vs Dim. ",i,sep=""),
           xlab=paste("Dim. ",i, " (",percentage[i],"%)",sep=""),
           ylab=paste("Dim.",i+1," (",percentage[i+1],"%)",sep=""))
      points(mdsFun[ mdsFun$Study!="PresentStudy",paste("X",i,sep="")],mdsFun[mdsFun$Study!="PresentStudy",paste("X",i+1,sep="")],
             pch=mdsFun$Point[ mdsFun$Study!="PresentStudy"],
             col=ifelse(mdsFun$Point[mdsFun$Study!="PresentStudy"]<21,mdsFun$Color[mdsFun$Study!="PresentStudy"],"black"),
             bg=mdsFun$Color[mdsFun$Study!="PresentStudy"],
             cex=mdsFun$cex[mdsFun$Study!="PresentStudy"])
      points(mdsFun[ mdsFun$Study=="PresentStudy",paste("X",i,sep="")],mdsFun[mdsFun$Study=="PresentStudy",paste("X",i+1,sep="")],
             pch=mdsFun$Point[ mdsFun$Study=="PresentStudy"],
             col=ifelse(mdsFun$Point[mdsFun$Study=="PresentStudy"]<21,mdsFun$Color[mdsFun$Study=="PresentStudy"],"black"),
             bg=mdsFun$Color[mdsFun$Study=="PresentStudy"],
             cex=mdsFun$cex[mdsFun$Study=="PresentStudy"],
             lwd=2)
      dev.off()
    }
  }else{
    for(i in seq(1,sum(str_starts(names(mdsFun),"X")),2)){
      svg(paste(pdfName,"_Dim",i,".Dim",i+1,".svg",sep=""))
      
      gg<-ggplot()+
        geom_point(x = as.numeric(mdsFun[,paste("X",i,sep="")]), y = as.numeric(mdsFun[,paste("X",i+1,sep="")]),
                   bg=mdsFun$Color,colour=ifelse(mdsFun$Point<21,mdsFun$Color,"black"),shape=mdsFun$Point,data=mdsFun,size=mdsFun$cex*2,
                   stroke=ifelse(mdsFun$Study=="PresentStudy",1.5,1))+
        coord_cartesian(xlim=c(min(mdsFun[,paste("X",i,sep="")]),max(mdsFun[,paste("X",i,sep="")])),
                        ylim=c(min(mdsFun[,paste("X",i+1,sep="")]),max(mdsFun[,paste("X",i+1,sep="")])))+
        geom_label_repel(data=mdsFun,x =mdsFun[,paste("X",i,sep="")], y = mdsFun[,paste("X",i+1,sep="")],
                         label=ifelse(mdsFun$Study=="PresentStudy",mdsFun$Population,NA),
                         #colour=ifelse(toPlot$Publication=="This Study","black","gray10"),
                         colour=mdsFun$Color,
                         fill=NA,
                         #colour="white",
                         #fill=toPlot$Color,
                         label.size=NA)+
        labs(title = Maintitle,x=paste("Dim.",i),y=paste("Dim.",i+1))+
        #annotate(geom="text",x=xright-rangeX/20,y=ytop-rangeY/40,label= "Color Coding",col="black",size=7)+
        theme_bw(base_size=10)
      print(gg)
      dev.off()
    }
  }
  svg(paste(pdfName,"_Legend.svg",sep=""))
  
  rascoLeg<-unique(mdsFun[mdsFun$Study=="PresentStudy",c("Population","Region","Point","Color")])
  plot(0,0,"n",axes=F,ann=F)
  rascoLeg<-rascoLeg[ order(rascoLeg$Region,rascoLeg$Population),]
  
  legend("center",pch=rascoLeg$Point,col=ifelse(rascoLeg$Point<21,rascoLeg$Color,"black"),
         pt.bg=rascoLeg$Color,
         legend = paste(rascoLeg$Region,": ",rascoLeg$Population,sep=""),
         ncol=2,
         cex=0.3) 
  plot(0,0,"n",axes=F,ann=F)
  
  Modern<-unique(mdsFun[ grepl("Modern",mdsFun$Region),c("Population","Region","Point","Color")])
  Modern<-Modern[ order(Modern$Region),]
  legend("center",pch=Modern$Point,col=ifelse(Modern$Point<21,Modern$Color,"black"),
         legend = paste(Modern$Region,": ",Modern$Population,sep=""),
         ncol=2,
         cex=0.5) 
  
  plot(0,0,"n",axes=F,ann=F)
  Ancient<-unique(mdsFun[ ! grepl("Modern",mdsFun$Region) & mdsFun$Study!="PresentStudy",c("Population","Region","Point","Color")])
  Ancient<-Ancient[ order(Ancient$Region),]
  legend("center",pch=Ancient$Point,pt.bg=Ancient$Color,col=ifelse(Ancient$Point<21,Ancient$Color,"black"),
         legend = paste(Ancient$Region,": ",Ancient$Population,sep=""),
         ncol=2,
         cex=0.3)
  
  dev.off()
}

change<-function(string,split,pos){
  return(strsplit(string,split = split,fixed=T)[[1]][pos])
  
}

###########################################################################
###########################################################################
###########################################################################
###########################################################################

ColorsLab<-archeo<-read.csv("../../Uspallata_Annotation.tsv",
                            stringsAsFactors = F,header=T,sep="\t",comment.char = "@")
names(ColorsLab)[names(ColorsLab)=="MainRegion"]<-"Region"
names(ColorsLab)[names(ColorsLab)=="Site"]<-"Population"
ColorsLab$Unit=ColorsLab$Individual
ColorsLab$cex=1.5
ColorsLab$Study="PresentStudy"
ColorsLab$Color[ColorsLab$Population=="Potrero Las Colonias"]<-"goldenrod1"
ColorsLab$Point[ColorsLab$Population=="Potrero Las Colonias"]<-21
ColorsLab$Region="Uspallata"

listModern=c("SouthernNorthAmerica_UnmaskedModern","CentralAmerica_UnmaskedModern","Amazonia_UnmaskedModern","Amazonia_Modern","CentralAndes_UnmaskedModern","CentralAndes_Modern","CentralChile_Modern","SouthPatagonia_Modern","Russia_Modern")
dicoSETS<-list("All"=c("Africa_Modern","LabMembers"),
               "America_NorthAsia"=c("Africa_Modern","LabMembers","Europe_Modern","Saami_LH","Saami_Modern","Laos_Hoabinhian_MH","India_Andamese_LH","Indonesia_Sulawesi_MH","China_Longlin_EH","China_RedDeerCave_LP","CentralSouthAsia_Modern","Oceania_Modern"),
               "America_NorthAsiaAncient"=c("Africa_Modern","LabMembers","Europe_Modern","Saami_LH","Saami_Modern","Laos_Hoabinhian_MH","India_Andamese_LH","Indonesia_Sulawesi_MH","China_Longlin_EH","China_RedDeerCave_LP","CentralSouthAsia_Modern","Oceania_Modern",
                                            listModern),
               "NoNorthernNorthAmerica"=c("Africa_Modern","LabMembers","Europe_Modern","Saami_LH","Saami_Modern","Laos_Hoabinhian_MH","India_Andamese_LH","Indonesia_Sulawesi_MH","China_Longlin_EH","China_RedDeerCave_LP","CentralSouthAsia_Modern","Oceania_Modern",
                                          "Russia_LP","Russia_EH","Russia_MH","Russia_LH","Russia_Modern","Beringia_EH","Beringia_LH","NorthernNorthAmerica_MH","NorthernNorthAmerica_LH"),
               "NoNorthernNorthAmericaAncient"=c("Africa_Modern","LabMembers","Europe_Modern","Saami_LH","Saami_Modern","Laos_Hoabinhian_MH","India_Andamese_LH","Indonesia_Sulawesi_MH","China_Longlin_EH","China_RedDeerCave_LP","CentralSouthAsia_Modern","Oceania_Modern",
                                                 "Russia_LP","Russia_EH","Russia_MH","Russia_LH","Russia_Modern","Beringia_EH","Beringia_LH","NorthernNorthAmerica_MH","NorthernNorthAmerica_LH",
                                                 listModern),
               "NoNorthernNorthAmerica.NoCarribes.Ancient"=c("Africa_Modern","LabMembers","Europe_Modern","Saami_LH","Saami_Modern","Laos_Hoabinhian_MH","India_Andamese_LH","Indonesia_Sulawesi_MH","China_Longlin_EH","China_RedDeerCave_LP","CentralSouthAsia_Modern","Oceania_Modern",
                                                 "Russia_LP","Russia_EH","Russia_MH","Russia_LH","Russia_Modern","Beringia_EH","Beringia_LH","NorthernNorthAmerica_MH","NorthernNorthAmerica_LH",
                                                  "ArchaicCarribean_LH","CeramicCarribean_LH",
                                                 listModern),
               "OnlySouthAmerica"=c("Africa_Modern","LabMembers","Europe_Modern","Saami_LH","Saami_Modern","Laos_Hoabinhian_MH","India_Andamese_LH","Indonesia_Sulawesi_MH","China_Longlin_EH","China_RedDeerCave_LP","CentralSouthAsia_Modern","Oceania_Modern","CentralAmerica_UnmaskedModern",
                                    "Russia_LP","Russia_EH","Russia_MH","Russia_LH","Russia_Modern","Beringia_EH","Beringia_LH","NorthernNorthAmerica_MH","NorthernNorthAmerica_LH",
                                    "SouthernNorthAmerica_EH","SouthernNorthAmerica_LH",
                                    "California_SouthernMainland_MH","California_NorthernChannelIslands_MH","California_SouthernChannelIslands_MH","California_Baja_LH","California_Central_LH","California_NorthernChannelIslands_LH","California_SouthernChannelIslands_LH","California_SouthernMainland_LH",
                                    "CentralAmerica_EH","CentralAmerica_MH","CentralAmerica_LH","ArchaicCarribean_LH","CeramicCarribean_LH"),
               
               "OnlySouthAmericaAncient"=c("Africa_Modern","LabMembers","Europe_Modern","Saami_LH","Saami_Modern","Laos_Hoabinhian_MH","India_Andamese_LH","Indonesia_Sulawesi_MH","China_Longlin_EH","China_RedDeerCave_LP","CentralSouthAsia_Modern","Oceania_Modern","CentralAmerica_UnmaskedModern",
                                           "Russia_LP","Russia_EH","Russia_MH","Russia_LH","Russia_Modern","Beringia_EH","Beringia_LH","NorthernNorthAmerica_MH","NorthernNorthAmerica_LH",
                                           "SouthernNorthAmerica_EH","SouthernNorthAmerica_LH",
                                           "California_SouthernMainland_MH","California_NorthernChannelIslands_MH","California_SouthernChannelIslands_MH","California_Baja_LH","California_Central_LH","California_NorthernChannelIslands_LH","California_SouthernChannelIslands_LH","California_SouthernMainland_LH",
                                           "CentralAmerica_EH","CentralAmerica_MH","CentralAmerica_LH","ArchaicCarribean_LH","CeramicCarribean_LH",
                                           listModern),
               "OnlySouthAmericaAncient_NoCentralAndes"=c("Africa_Modern","LabMembers","Europe_Modern","Saami_LH","Saami_Modern","Laos_Hoabinhian_MH","India_Andamese_LH","Indonesia_Sulawesi_MH","China_Longlin_EH","China_RedDeerCave_LP","CentralSouthAsia_Modern","Oceania_Modern","CentralAmerica_UnmaskedModern",
                                                          "Russia_LP","Russia_EH","Russia_MH","Russia_LH","Russia_Modern","Beringia_EH","Beringia_LH","NorthernNorthAmerica_MH","NorthernNorthAmerica_LH",
                                                          "SouthernNorthAmerica_EH","SouthernNorthAmerica_LH",
                                                          "California_SouthernMainland_MH","California_NorthernChannelIslands_MH","California_SouthernChannelIslands_MH","California_Baja_LH","California_Central_LH","California_NorthernChannelIslands_LH","California_SouthernChannelIslands_LH","California_SouthernMainland_LH",
                                                          "CentralAmerica_EH","CentralAmerica_MH","CentralAmerica_LH","ArchaicCarribean_LH","CeramicCarribean_LH",
                                                          listModern,
                                                          "CentralAndes_EH","CentralAndes_MH","CentralAndes_EarlyIntermediate_LH","CentralAndes_LateIntermediate_LH","CentralAndes_MiddleHorizonLateIntermeditae_LH","CentralAndes_MiddleHorizon_LH","CentralAndes_LateHorizon_LH","CentralAndes_Unknown")
) 



listINDremove=read.table("../../Kinship/listToRemove_2ndDegree.txt",stringsAsFactors=F,header=T)
listINDremove<-c("Brazil-12",listINDremove$Individual)
listCarribRM<-read.table("../../..//ReferenceDataSet/CompendiumPostKin_Analyses_NoReichSG/ListCaribeanSubSample_TORM.txt",stringsAsFactors=F,header=F)$V2
listINDremove<-c(listCarribRM,listINDremove)
#print(listINDremove)



###########START!



outTABLE<-c()
outMDS<-c()
setind<-"Lab_with_Compendium_GEHmodern"
TH<-50000

ColorsALL<-c()
for(setsnp in c("1240K","1240K.TVs","SG","SG.TVs")[1]){
  f3<-read.table(paste("../F3_IND/",setind,".",setsnp,"/TH",TH,"/",setind,".",setsnp,".OUT",sep=""),stringsAsFactors = F,header=T)
  f3<-f3[ ! (f3$Source1 %in% listINDremove | f3$Source2 %in% listINDremove),]
  famind<-read.table(paste("../DataSets/",setind,".",setsnp,"/","/finalSet",ifelse(TH==0,"",paste(".TH",TH,sep="")),".ind.txt",sep=""),stringsAsFactors = F,header=F)
  famind<-famind[,c(3,1)]
  names(famind)<-c("Population","Individual")
  famind$Population<-str_remove(famind$Population,".DG")
  ColorsCompendium<-read.table("../../..//ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv",stringsAsFactors = F,header=T,sep="\t",comment.char = "@")
  ColorsCompendium$cex=1
  ColorsCompendium$Study="REF"
  #ColorsCompendium<-unique(ColorsCompendium[ColorsCompendium$Study!="DifferentProjects",c("Region","Population","Study","Color","Point","cex")])
  ColorsCompendium<-unique(ColorsCompendium[,c("Region","Population","Study","Color","Point","cex")])
  ColorsCompendium$Unit=ColorsCompendium$Population
  ColorsCompendium<-merge(ColorsCompendium,famind,by="Population")
  
  Colors<-rbind(ColorsLab[c("Individual","Region","Population","Study","Color","Point","cex")],
                ColorsCompendium[c("Individual","Region","Population","Study","Color","Point","cex")])
  Colors$Point<-as.numeric(Colors$Point)
  
  
  Colors<-rbind(ColorsLab[c("Individual","Region","Population","Study","Color","Point","cex")],
                ColorsCompendium[c("Individual","Region","Population","Study","Color","Point","cex")])
  
  Colors$Point<-as.numeric(Colors$Point)
  ColorsALL<-unique(rbind(Colors,ColorsALL))  
  
  
  readDM<-readLines(paste("../F3_IND/",setind,".",setsnp,"/TH",TH,"/",setind,".",setsnp,"_MDS.dist",sep=""))
  dm<-matrix(NA,length(readDM)-1,length(readDM)-1)
  colnames(dm)<-strsplit(readDM[1],split="\t")[[1]]
  rownamesDM<-c()
  for(i in c(2:length(readDM))){
    tmp<-strsplit(readDM[i],split="\t")[[1]]
    rownamesDM<-c(rownamesDM,tmp[1])
    dm[i-1,]<-as.numeric(tmp[-1])
  }
  rownames(dm)<-rownamesDM
  
  dm<-dm[ !rownames(dm) %in% listINDremove , ! colnames(dm) %in% listINDremove]
  
  #for(subset in c("PresentStudy",names(dicoSETS))){
  for(subset in c("NoNorthernNorthAmerica.NoCarribes.Ancient")){
    if(subset != "PresentStudy"){
      if(grepl("SG",setsnp)){
        listRemove<-Colors$Population[  str_ends(Colors$Population,".Capt") ] 
      }else{
        listRemove<-c()
      }
      for(remPattern in dicoSETS[[subset]]){
        listRemove<-c(listRemove,unique(Colors$Population[  grepl(remPattern,Colors$Region) ]))
      }
      #print(paste("removing those pops for ",subset))
      #print(listRemove)
    }else{
      listRemove<-unique(Colors$Population[  Colors$Study!="PresentStudy" ])
    }
    
    INDSub<-famind$Individual[ ! famind$Population %in% listRemove]
    #print(paste("those regions remain:"))
    #print(unique(sort(Colors$Region[ Colors$Individual %in% INDSub])))
    dmSub<-dm[ rownames(dm) %in% INDSub, colnames(dm) %in% INDSub ]
    
    mdsOUT<-cmdscale(dmSub,12,eig=T)
    mds<-data.frame(mdsOUT$points)
    perc<-round(100*mdsOUT$eig/sum(mdsOUT$eig),digits = 2)
    mds$Individual<-row.names(mds)
    print(famind[ famind$Individual %in% mds$Individual[ !mds$Individual %in% Colors$Individual],])
    mds<-merge(mds,Colors,by="Individual")
    
    plotMDS(mds,
            setsnp,
            perc,
            #paste("F3_IND/",setind,".",setsnp,"_MDS.",subset,".pdf",sep=""),F)
            paste("F3_IND/",setind,".",setsnp,"_MDS.",subset,".svg",sep=""),F)
    #write.table(mds,paste("F3_IND/",setind,".",setsnp,"_MDS.",subset,".tsv",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
  
    f3<-read.table(paste("../F3_IND/",setind,".",setsnp,"/TH",TH,"/",setind,".",setsnp,".OUT",sep=""),stringsAsFactors = F,header=T)
    f3<-f3[ ! (f3$Source1 %in% listINDremove | f3$Source2 %in% listINDremove),]
    
    names(f3)[4:7]<-paste(names(f3)[4:7],setsnp,sep="_")
    f3<-f3[ ! (f3$Source1 %in% famind$Individual[famind$Population %in% listRemove] | f3$Source2 %in% famind$Individual[famind$Population %in% listRemove]),]
    
    names(mds)[c(2:13)]<-paste("Dim",c(1:12),"_",setsnp,sep="")
    mds<-mds[, names(mds)!="Study"]
    if(length(outTABLE[[subset]])==0){
      outTABLE[[subset]]<-f3  
      outMDS[[subset]]<-mds
    }else{
      outTABLE[[subset]]<-merge(outTABLE[[subset]],f3,by=c("Source1","Source2","Target"),all=T)
      outMDS[[subset]]<-merge(outMDS[[subset]],mds,by=c("Individual","Region","Population","Color","Point","cex"),all=T)
    }
    
    if(length(mds$Individual)!=length(unique(c(f3$Source1,f3$Source2)))){
      stop()
    }
  
  }
}

for(subset in c("NoNorthernNorthAmerica.NoCarribes.Ancient")){
  print(nrow(outTABLE[[subset]]))
  outTABLE[[subset]]<-merge(outTABLE[[subset]],ColorsALL[,names(ColorsALL)!="Study"],by.x="Source1",by.y="Individual")
  names(outTABLE[[subset]])[length(outTABLE[[subset]])-c(4:0)]<-paste(names(outTABLE[[subset]])[length(outTABLE[[subset]])-c(4:0)],"1",sep="")
  outTABLE[[subset]]<-merge(outTABLE[[subset]],ColorsALL[,names(ColorsALL)!="Study"],by.x="Source2",by.y="Individual")
  names(outTABLE[[subset]])[length(outTABLE[[subset]])-c(4:0)]<-paste(names(outTABLE[[subset]])[length(outTABLE[[subset]])-c(4:0)],"2",sep="")
  print(nrow(outTABLE[[subset]]))
  
  write.table(outTABLE[[subset]],paste("F3_IND/",setind,"_F3byIND.",subset,".tsv",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
  write.table(outMDS[[subset]],paste("F3_IND/",setind,"_MDS.",subset,".tsv",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
}



