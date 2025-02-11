#!/bin/Rscript

Sys.setlocale( 'LC_ALL','C' ) 
require(stringr)

library(ggrepel)
library(ggplot2)


require(ape)

setwd("~/Documents/PostDocPasteur/aDNA//2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F/F3_IND")
setind<-"Lab_with_Compendium_GEHmodern"
TH<-"50000"
subset="NoNorthernNorthAmerica.NoCarribes.Ancient"

###########################################################################
###########################################################################
###########################################################################
###########################################################################
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

dicoOUT<-list("All"="DA236",
              "America_NorthAsia"="DA236",
              "America_NorthAsiaAncient"="DA236",
              "PresentStudy"="USR1",
              "NoNorthernNorthAmerica"="USR1",
              "NoNorthernNorthAmericaAncient"="USR1",
              "NoNorthernNorthAmerica.NoCarribes.Ancient"="USR1",
              "OnlySouthAmerica"="USR1",
              "OnlySouthAmerica_NoCentralAndes"="USR1",
              "OnlySouthAmericaAncient"="USR1",
              "OnlySouthAmericaAncient_NoCentralAndes"="USR1")



ColorsLab<-archeo<-read.csv("../../../Uspallata_Annotation.tsv",
                            stringsAsFactors = F,header=T,sep="\t",comment.char = "@")
names(ColorsLab)[names(ColorsLab)=="MainRegion"]<-"Region"
names(ColorsLab)[names(ColorsLab)=="Site"]<-"Population"
ColorsLab$Unit=ColorsLab$Individual
ColorsLab$cex=1.5
ColorsLab$Study="PresentStudy"




listINDremove=read.table("../../../Kinship/listToRemove_2ndDegree.txt",stringsAsFactors=F,header=T)
listINDremove<-c("Brazil-12",listINDremove$Individual[ listINDremove$Remove ])
listCarribRM<-read.table("../../../../ReferenceDataSet/CompendiumPostKin_Analyses_NoReichSG/ListCaribeanSubSample_TORM.txt",stringsAsFactors=F,header=F)$V2
listINDremove<-c(listCarribRM,listINDremove)
print(listINDremove)


for(setsnp in c("1240K","1240K.TVs","SG","SG.TVs")){
  f3<-read.table(paste("../../F3_IND/",setind,".",setsnp,"/TH",TH,"/",setind,".",setsnp,".OUT",sep=""),stringsAsFactors = F,header=T)
  
  famind<-read.table(paste("../../DataSets/",setind,".",setsnp,"/","/finalSet",ifelse(TH==0,"",paste(".TH",TH,sep="")),".ind.txt",sep=""),stringsAsFactors = F,header=F)
  famind<-famind[,c(3,1)]
  names(famind)<-c("Population","Individual")
  
  ColorsCompendium<-read.table("../../../..//ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv",stringsAsFactors = F,header=T,sep="\t",comment.char = "@")
  ColorsCompendium$cex=1
  #ColorsCompendium$Population=str_remove(str_remove(str_remove(ColorsCompendium$Population,".SG"),".Capt"),".Mix")
  ColorsCompendium$Unit=ColorsCompendium$Population
  ColorsCompendium<-unique(ColorsCompendium[,c("Region","Population","Study","Color","Point","cex")])
  #ColorsCompendium<-unique(ColorsCompendium[ColorsCompendium$Study!="DifferentProjects",c("Region","Population","Study","Color","Point","cex")])
  
  ColorsCompendium<-merge(ColorsCompendium,famind,by="Population")
  
  Colors<-rbind(ColorsLab[c("Individual","Region","Population","Study","Color","Point","cex")],
                ColorsCompendium[c("Individual","Region","Population","Study","Color","Point","cex")])
  print("here") 
  rownames(Colors)<-Colors$Individual
  print("there")
  f3<-f3[ ! f3$Source1 %in% famind$Individual[ famind$Population=="Ignore"] &
            ! f3$Source2 %in% famind$Individual[ famind$Population=="Ignore"],  ]
  listPop<-unique(c(f3$Source1,f3$Source2))
  
  if(! file.exists(paste("../../F3_IND/",setind,".",setsnp,"/TH",TH,"/",setind,".",setsnp,"_TREE.dist",sep=""))){
    dm<-data.frame(matrix(NA,length(listPop),length(listPop)))
    names(dm)<-listPop
    rownames(dm)<-listPop
    for(p1 in c(1:(length(listPop)-1))){
      print(p1)
      dm[p1,p1]<-0
      for(p2 in c((p1+1):length(listPop))){
        dm[p2,p1]<-dm[p1,p2]<-1 / na.omit(c(f3$f_3[ f3$Source1==listPop[p1] & f3$Source2==listPop[p2]],f3$f_3[ f3$Source1==listPop[p2] & f3$Source2==listPop[p1]]))
      }
      
    }
    p1=p1+1
    dm[p1,p1]<-0
    
    write.table(dm,paste(setind,".",setsnp,"/TH",TH,"/",setind,".",setsnp,"_TREE.dist",sep=""),sep="\t",quote=F,row.names=T,col.names=T)
  }else{
    readDM<-readLines(paste("../../F3_IND/",setind,".",setsnp,"/TH",TH,"/",setind,".",setsnp,"_TREE.dist",sep=""))
    dm<-matrix(NA,length(readDM)-1,length(readDM)-1)
    colnames(dm)<-strsplit(readDM[1],split="\t")[[1]]
    rownamesDM<-c()
    for(i in c(2:length(readDM))){
      tmp<-strsplit(readDM[i],split="\t")[[1]]
      rownamesDM<-c(rownamesDM,tmp[1])
      dm[i-1,]<-as.numeric(tmp[-1])
    }
    rownames(dm)<-rownamesDM
    
  }
  dm<-dm[ !rownames(dm) %in% listINDremove , ! colnames(dm) %in% listINDremove]
  
  if(subset != "PresentStudy"){
    listRemove<-c()
    for(remPattern in dicoSETS[[subset]]){
      listRemove<-c(listRemove,unique(Colors$Population[  grepl(remPattern,Colors$Region) ]))
    }
    print(paste("removing those pops for ",subset))
    print(listRemove)
  }else{
    listRemove<-unique(Colors$Population[  Colors$Study!="PresentStudy" ])
  }
  if(grepl("SG",setsnp)){
    listRemove<-unique(c(listRemove,Colors$Population[  str_ends(Colors$Population,".Capt")]))
  }

    
    INDSub<-famind$Individual[ ! famind$Population %in% listRemove | famind$Individual==dicoOUT[[subset]]]
    dmSub<-dm[ rownames(dm) %in% INDSub, colnames(dm) %in% INDSub ]
    
    dicoPLOT<-Colors[ rownames(dmSub),]
  
    
    pdf(paste(setind,".",setsnp,"_TREE.",subset,".pdf",sep=""),height=ifelse(grepl("SG",setsnp),15,30),width=20)
    dist<-as.dist(dmSub)
    nj<-bionj(dist)
    nj<-root(nj,outgroup = dicoOUT[[subset]],resolve.root = F)
    write.tree(nj,file=paste(setind,".",setsnp,"_TREE.",subset,".nwk",sep=""),append = F,tree.names = F)
    nj<-read.tree(paste(setind,".",setsnp,"_TREE.",subset,".nwk",sep=""))

    if(sum(! nj$tip.label %in% row.names(Colors) )>0){
      print(paste(nj$tip.label[! nj$tip.label %in% row.names(Colors)]))
      stop("pb individual not in your dico")
    }
    dicoPLOT<-Colors[nj$tip.label,]
    nj$tip.label=paste(dicoPLOT$Region,": ",dicoPLOT$Population," (",nj$tip.label,")",sep="")
    plotStats<-plot(nj,tip.color = dicoPLOT$Color,cex = 0.3,use.edge.length = F,show.tip.label=T,label.offset=10,align.tip.label=T)
    points(x=rep(nrow(dicoPLOT)+5,nrow(dicoPLOT)),y=c(1:nrow(dicoPLOT)),
       cex=0.7,pch=dicoPLOT$Point,bg=dicoPLOT$Color,col=ifelse(dicoPLOT$Point<21,dicoPLOT$Color,"black"))

    dev.off()
    outDISTname=paste(setind,".",setsnp,"_TREE.",subset,".NJdist",sep="")
    system(paste("echo \"\t\"", ncol(dmSub), " > ",outDISTname,sep=""))
    write.table(dmSub,outDISTname,col.names = F,row.names = T,append=T,sep="\t",quote=F)
    write.table(dicoPLOT,paste(outDISTname,".dico",sep=""),col.names = T,row.names = F,append=F,sep="\t",quote=F)
  }
