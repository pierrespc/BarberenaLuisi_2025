#!/bin/Rscript
#library(rnaturalearth)
require("ggplot2")
require("sf")
library(ggspatial)


PUTMODERN="NONE" #ONLY NONE BOTH

setwd("~/Documents/PostDocPasteur/aDNA//2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F/Maps/")
annot<-read.table("../../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv",stringsAsFactors = F,header=T,sep="\t")

annot<-annot[,c("Study","Region","Population","Latitude","Longitude")]


indLIST<-read.table("AnnotIndividualsIncluded.tsv",stringsAsFactors = F,header=T,sep="\t")


uspallata<-data.frame(cbind(lat=c(-28.852526,-28.852526,-34.304672,-34.304672),
                            long=c(-72.346136,-67.334717,-67.334717,-72.346136)))
  
coordSITES<-data.frame(rbind(c("Barrio Ramos I",-32.591333,-69.344217),
                               c("Potrero Las Colonias",-32.60925,-69.362067),
                               c("Tumulo II",-32.587167,-69.360167),
                               c("Barrancas B6",	-33.09325,-68.746417)))
names(coordSITES)<-c("Site","Latitude","Longitude")
coordSITES$Latitude<-as.numeric(coordSITES$Latitude)
coordSITES$Longitude<-as.numeric(coordSITES$Longitude)
AnnotStudy<-read.table("../../../Uspallata_Annotation.tsv",stringsAsFactors = F,header=T,sep="\t",comment.char = "@",quote="\"")
AnnotStudy<-merge(AnnotStudy[,c("MainRegion","Site")],coordSITES,by="Site")

#AnnotStudy<-AnnotStudy[,-1]
AnnotStudy$Region="Uspallata"
AnnotStudy$MainRegion<-ifelse(AnnotStudy$MainRegion=="Local_Early",paste(AnnotStudy$Site," (LH-HG|LH-HG)",sep=""),
                              ifelse(AnnotStudy$MainRegion=="Local_Late",paste(AnnotStudy$Site," (LH-LF|LH-F)",sep=""),
                                     ifelse(AnnotStudy$MainRegion=="Migrant",paste(AnnotStudy$Site," (LH-MF|LH-F)",sep=""),
                                            ifelse(AnnotStudy$MainRegion=="Migrant_Outlier",paste(AnnotStudy$Site," (LH-MFout|LH-F)",sep=""),
                                                   ifelse(AnnotStudy$MainRegion=="Local_PLC",paste(AnnotStudy$Site, " (LH-LFplc|LH-F)",sep=""),"WTF")))))

AnnotStudy<-AnnotStudy[,-1]
names(AnnotStudy)[names(AnnotStudy)=="MainRegion"]<-"Population"
AnnotStudy$Study="PresentStudy"
annot<-rbind(annot,unique(AnnotStudy[! duplicated(AnnotStudy$Population),names(annot)]))

annot<-merge(annot,indLIST,by.y=c("Study","MajorGroup","Group"),by.x=c("Study","Region","Population"))

if(nrow(annot)!=nrow(indLIST)){
  stop("issue")
}



annot<-unique(annot[,names(annot) %in% c("Study","Region","Population","Latitude","Longitude","Color","Point")])
annot$Latitude=as.numeric(annot$Latitude)
annot$Longitude=as.numeric(annot$Longitude)




annotSimplified<-c()
annotOverSimplified<-c()
for(reg in unique(annot$Region)){
  tmp<-annot[ annot$Region==reg,]
  tmp$Point<-tmp$Point[1]
  tmp$Color<-tmp$Color[1]
  tmp$LongitudeMEAN=mean(tmp$Longitude)
  tmp$LatitudeMEAN=mean(tmp$Latitude)
  tmp$Label=paste(tmp$Region," (",nrow(tmp)," Groups)",sep="")
  annotSimplified<-rbind(annotSimplified,unique(tmp[,c("Color","Point","Label","Latitude","Longitude","Region")]))  
  annotOverSimplified<-rbind(annotOverSimplified,unique(tmp[,c("Color","Point","Label","LatitudeMEAN","LongitudeMEAN","Region")]))  
}
annotSimplified$Point=as.numeric(annotSimplified$Point)
annotOverSimplified$Point=as.numeric(annotOverSimplified$Point)



world_map <- ne_countries(scale = "medium", returnclass = "sf")

#world_map <- world_map %>% dplyr::group_by(group)  %>% st_cast("LINESTRING")
  



annotSimplifiedTOWRITE<-merge(annot[,!names(annot) %in% c("Point","Color")],unique(annotSimplified[,c("Region","Point","Color")]),by="Region")

if(PUTMODERN =="NONE"){
    annot<-annot[ ! grepl("Modern",annot$Region),]
    system("mkdir OnlyAncient/")
    ncolLeg=3
}else{
  if(PUTMODERN =="BOTH"){
    system("mkdir WithModern/")
    ncolLeg=4
  }else{
    if(PUTMODERN =="ONLY"){
      annot<-annot[  grepl("Modern",annot$Region),]
      system("mkdir OnlyModern/")
      ncolLeg=4
    }
  }
}
annot<-annot[ ! is.na(annot$Latitude),]

listRegions<-read.table("listRegionOrdered_withGenoModern.txt",stringsAsFactors = F,header=F)$V1

###remove ancient in study region
if(PUTMODERN!="ONLY"){
  annotPLOT<-annot[ !(annot$Longitude> min(uspallata$long) & annot$Longitude< max(uspallata$long) &
                        annot$Latitude> min(uspallata$lat) & annot$Latitude< max(uspallata$lat)),  ]
  annotSimplifiedPLOT<-annotSimplified[ !(annotSimplified$Longitude> min(uspallata$long) & annotSimplified$Longitude< max(uspallata$long) &
                                            annotSimplified$Latitude> min(uspallata$lat) & annotSimplified$Latitude< max(uspallata$lat)),  ]
  annotOverSimplifiedPLOT<-annotOverSimplified[ !(annotOverSimplified$Longitude> min(uspallata$long) & annotOverSimplified$Longitude< max(uspallata$long) &
                                                    annotOverSimplified$Latitude> min(uspallata$lat) & annotOverSimplified$Latitude< max(uspallata$lat)),  ]
  
}else{
  annotPLOT<-annot  
    
}



for(ZONE in c("SouthAmerica","America","Global")[2]){
  if(ZONE == "SouthAmerica"){
    yrange=c(-56,0)
    xrange=c(-85,-35)
  }else{
    if(ZONE =="America"){
      yrange=c(-56,47)
      xrange=c(-123,-35)
    }else{
      yrange=c(-85.19218,83.59961)
      xrange=c(-180.0000,190.2708)
    }
  }
  
  #mapCroped<-st_crop(world_map, xmin = xrange[1], xmax = xrange[2], ymin = yrange[1], ymax = yrange[2])

  
  base<-ggplot()+
    labs(title="")+
    theme_classic()+
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank())

  
  Poly<-base+
    geom_sf(data=world_map,fill = "#D3D3D3", color = "white", size = 0.2)+
    coord_sf(xlim=xrange,ylim=yrange)
  
  
  if(sum(!annot$Region %in% listRegions)){
    print(annot$Region[!annot$Region %in% listRegions])
    stop()
  }
  colorScale<-c()
  fillScale<-c()
  shapeScale<-c()
  
  for(i in unique(annot$Population[annot$Latitude>yrange[1] & annot$Latitude<yrange[2] &
                                   annot$Longitude>xrange[1] & annot$Longitude<xrange[2]])){
    colorScale[i]<-ifelse(annot$Point[ annot$Population==i]<21,annot$Color[ annot$Population==i],"black")
    fillScale[i]<-annot$Color[ annot$Population==i]
    shapeScale[i]<-annot$Point[ annot$Population==i]
  }
  
  
  Points<-base+
    geom_point(data=annotPLOT,
               #geom_jitter(data=annotPLOT,width = 2,height=2,
               #size=ifelse(annotPLOT$Region=="This Study",5,3),stroke=1.3,
               size=ifelse(annotPLOT$Region=="This Study",2,2),stroke=1,
               mapping=aes(y=Latitude,x=Longitude,fill=Population,color=Population,shape=Population))+
    scale_color_manual(values=colorScale)+
    scale_fill_manual(values=fillScale)+
    scale_shape_manual(values=shapeScale)+
    geom_polygon(data=uspallata,fill="grey40",color="grey40",mapping=aes(x=long,y=lat))+
    coord_sf(xlim=xrange,ylim=yrange)
  
  pdf(paste(ifelse(PUTMODERN=="BOTH","WithModern/",
  #svg(paste(ifelse(PUTMODERN=="BOTH","WithModern/",
                   ifelse(PUTMODERN=="NONE","OnlyAncient/",
                          ifelse(PUTMODERN=="ONLY","OnlyModern/",stop()))),"MapComplex_",ZONE,".pdf",sep=""))
                          #ifelse(PUTMODERN=="ONLY","OnlyModern/",stop()))),"MapComplex_",ZONE,".svg",sep=""))
  print(base+
          geom_polygon(data=uspallata,fill="grey50",color="grey50",mapping=aes(x=long,y=lat))+
          geom_sf(data=world_map,fill = NA, color = "grey70", size = 0.4)+
          geom_point(data=annotPLOT,
          #geom_jitter(data=annotPLOT,width = 2,height=2,
                      #size=ifelse(annotPLOT$Region=="This Study",5,3),stroke=1.3,
                      size=ifelse(annotPLOT$Region=="This Study",2,2),stroke=1,
                      mapping=aes(y=Latitude,x=Longitude,fill=Population,color=Population,shape=Population))+
          scale_color_manual(values=colorScale)+
          scale_fill_manual(values=fillScale)+
          scale_shape_manual(values=shapeScale)+
          coord_sf(xlim=xrange,ylim=yrange))
  
  #print(Poly)
  #print(Points)
  dev.off()
  pdf(paste("",ifelse(PUTMODERN=="BOTH","WithModern/",
                             ifelse(PUTMODERN=="NONE","OnlyAncient/",
                                    ifelse(PUTMODERN=="ONLY","OnlyModern/",stop()))),"LegendComplex_",ZONE,".pdf",sep=""),width=18,height=14)
  Leg<-c()
  for(region in listRegions[  listRegions %in% annot$Region[ annot$Latitude>yrange[1] & annot$Latitude<yrange[2] &
                                                             annot$Longitude>xrange[1] & annot$Longitude<xrange[2]]]){
    tmp<-annot[ annot$Region==region,]
    tmp<-tmp[order(tmp$Color,tmp$Point),]
    tmp$font=1
    tmp$cex=0.5
    tmp$textColor="black"
    Leg<-rbind(Leg,cbind(Population=region,Color=NA,Point=NA,font=3,cex=0.7,textColor=ifelse(tmp$Color[1]=="grey100","grey90",tmp$Color[1])),
               tmp[,c("Population","Color","Point","font","cex","textColor")])
  }
  Leg$Point<-as.numeric(Leg$Point)
  Leg$cex<-as.numeric(Leg$cex)
  Leg$font<-as.numeric(Leg$font)
  plot(0,0,"n",axes=F,ann=F)
  legend("center",pch=Leg$Point,col=ifelse(Leg$Point<21,Leg$Color,"black"),
         pt.bg=Leg$Color,text.col = Leg$textColor,text.font = Leg$font,legend=Leg$Population,cex=Leg$cex,ncol=ncolLeg,box.col=NA)
  dev.off()
  
  next
  ##############NOW SIMPLIFIED
  
  colorScale<-c()
  fillScale<-c()
  shapeScale<-c()
  
  for(i in unique(annotSimplified$Label[annotSimplified$Latitude>yrange[1] & annotSimplified$Latitude<yrange[2] &
                                        annotSimplified$Longitude>xrange[1] & annotSimplified$Longitude<xrange[2]])){
    colorScale[i]<-ifelse(annotSimplified$Point[ annotSimplified$Label==i]<21,annotSimplified$Color[ annotSimplified$Label==i],"black")
    fillScale[i]<-annotSimplified$Color[ annotSimplified$Label==i]
    shapeScale[i]<-annotSimplified$Point[ annotSimplified$Label==i]
  }
  pdf(paste("",ifelse(PUTMODERN=="BOTH","WithModern/",
                      ifelse(PUTMODERN=="NONE","OnlyAncient/",
                             ifelse(PUTMODERN=="ONLY","OnlyModern/",stop()))),"MapSimplified_",ZONE,".pdf",sep=""))
  print(base+
          geom_point(data=annotSimplifiedPLOT,
                     size=ifelse(annotSimplifiedPLOT$Region=="This Study",5,3),stroke=1.3,,
                     mapping=aes(y=Latitude,x=Longitude,fill=Label,color=Label,shape=Label))+
          scale_color_manual(values=colorScale)+
          scale_fill_manual(values=fillScale)+
          scale_shape_manual(values=shapeScale))
  dev.off()
  
  pdf(paste("",ifelse(PUTMODERN=="BOTH","WithModern/",
                      ifelse(PUTMODERN=="NONE","OnlyAncient/",
                             ifelse(PUTMODERN=="ONLY","OnlyModern/",stop()))),"LegendSimplified_",ZONE,".pdf",sep=""))
  Leg<-c()
  for(region in listRegions[  listRegions %in% annotSimplified$Region[ annotSimplified$Latitude>yrange[1] & annotSimplified$Latitude<yrange[2] &
                                                             annotSimplified$Longitude>xrange[1] & annotSimplified$Longitude<xrange[2]]]){
    tmp<-annotSimplified[ annotSimplified$Region==region,]
    tmp<-tmp[order(tmp$Color,tmp$Point),]
    tmp$font=1
    tmp$cex=0.5
    tmp$textColor="black"
    Leg<-rbind(Leg,unique(tmp[,c("Region","Color","Point","font","cex","textColor")]))
  }
  Leg$Point<-as.numeric(Leg$Point)
  Leg$cex<-as.numeric(Leg$cex)
  Leg$font<-as.numeric(Leg$font)
  plot(0,0,"n",axes=F,ann=F)
  legend("center",pch=Leg$Point,col=ifelse(Leg$Point<21,Leg$Color,"black"),
         pt.bg=Leg$Color,text.col = Leg$textColor,text.font = Leg$font,legend=Leg$Region,cex=Leg$cex,box.col=NA)
  dev.off()
  
  ##############NOW OVERSIMPLIFIED
  colorScale<-c()
  fillScale<-c()
  shapeScale<-c()
  
  for(i in unique(annotOverSimplified$Label[annotOverSimplified$LatitudeMEAN>yrange[1] & annotOverSimplified$LatitudeMEAN<yrange[2] &
                                        annotOverSimplified$LongitudeMEAN>xrange[1] & annotOverSimplified$LongitudeMEAN<xrange[2]])){
    colorScale[i]<-ifelse(annotOverSimplified$Point[ annotOverSimplified$Label==i]<21,annotOverSimplified$Color[ annotOverSimplified$Label==i],"black")
    fillScale[i]<-annotOverSimplified$Color[ annotOverSimplified$Label==i]
    shapeScale[i]<-annotOverSimplified$Point[ annotOverSimplified$Label==i]
  }
  pdf(paste("",ifelse(PUTMODERN=="BOTH","WithModern/",
                      ifelse(PUTMODERN=="NONE","OnlyAncient/",
                             ifelse(PUTMODERN=="ONLY","OnlyModern/",stop()))),"MapOverSimplified_",ZONE,".pdf",sep=""))
  print(base+
          geom_point(data=annotOverSimplifiedPLOT,
                     size=ifelse(annotOverSimplifiedPLOT$Region=="This Study",5,3),stroke=1.3,
                     mapping=aes(y=LatitudeMEAN,x=LongitudeMEAN,fill=Label,color=Label,shape=Label))+
          scale_color_manual(values=colorScale)+
          scale_fill_manual(values=fillScale)+
          scale_shape_manual(values=shapeScale))
  dev.off()
  
  pdf(paste("",ifelse(PUTMODERN=="BOTH","WithModern/",
                      ifelse(PUTMODERN=="NONE","OnlyAncient/",
                             ifelse(PUTMODERN=="ONLY","OnlyModern/",stop()))),"LegendOverSimplified_",ZONE,".pdf",sep=""))
  Leg<-c()
  for(region in listRegions[  listRegions %in% annotOverSimplified$Region[ annotOverSimplified$LatitudeMEAN>yrange[1] & annotOverSimplified$LatitudeMEAN<yrange[2] &
                                                                       annotOverSimplified$LongitudeMEAN>xrange[1] & annotOverSimplified$LongitudeMEAN<xrange[2]]]){
    tmp<-annotOverSimplified[ annotOverSimplified$Region==region,]
    tmp<-tmp[order(tmp$Color,tmp$Point),]
    tmp$font=1
    tmp$cex=0.5
    tmp$textColor="black"
    Leg<-rbind(Leg,unique(tmp[,c("Region","Color","Point","font","cex","textColor")]))
  }
  Leg$Point<-as.numeric(Leg$Point)
  Leg$cex<-as.numeric(Leg$cex)
  Leg$font<-as.numeric(Leg$font)
  plot(0,0,"n",axes=F,ann=F)
  legend("center",pch=Leg$Point,col=ifelse(Leg$Point<21,Leg$Color,"black"),
         pt.bg=Leg$Color,text.col = Leg$textColor,text.font = Leg$font,legend=Leg$Region,cex=Leg$cex,box.col = NA)
  dev.off()
}


REG<-Leg[ !Leg$Population %in% annotPLOT$Population & ! is.na(Leg$Point),]
REG<-REG[ !grepl("LFplc",REG$Population),]
REG<-REG[ !grepl("MFout",REG$Population),]
REG$Population[ grepl("Potrero",REG$Population)]<-"Potrero las Colonias"
REG$Population[ grepl("Tumulo II",REG$Population)]<-"Túmulo II"
REG$Population[ grepl("Barrio",REG$Population)]<-"Barrio Ramos I"
REG$Population[ grepl("Barrancas",REG$Population)]<-"Barrancas B6"

REG<-rbind(REG,cbind("Population"=c("Túmulo I","Túmulo III","A. Avispas","Usina Sur 2"),
           "Color"=c("darkorange4","darkorange4","goldenrod1","goldenrod1"),
           "Point"=c(24,25,22,23),
           "font"=1,
           "cex"=0.5,
           "textColor"="black"))

REG$Point=as.numeric(REG$Point)
REG$Population[ REG$Population=="Chile_LosRieles_12000BP.Capt"]<-"Los Rieles  ~12000BP (Central Chile EH)"
REG$Population[ REG$Population=="Chile_LosRieles_5100BP.Capt"]<-"Los Rieles ~5100BP (Central Chile MH)"
REG$Population[ REG$Population=="Chile_Conchali_700BP.Capt"]<-"Conchalí ~700BP (Central Chile LH)"
REG$Population[ REG$Population=="Argentina_Aconcagua_Inca_500BP.SG"]<-"Aconcagua ~500BP (Qhapaq Hucha Inka LH)"
REG$Population[ REG$Population=="Argentina_Toro_475BP.SG"]<-"El Toro Mountain ~475BP (Qhapaq Hucha Inka LH)"
REG$Population[ REG$Population=="Chile_Plomo_490BP.SG"]<-"El Plomo Mountain ~490BP (Qhapaq Hucha Inka LH)"
REG$Population[ REG$Population=="Argentina_Callingasta_1500BP.SG"]<-"Calingasta Valley ~1500BP (Cuyo LH)"
pdf("LegZONE.pdf")
plot(0,0,"n",axes=F,ann=F)
legend("center",pch=REG$Point,col=ifelse(REG$Point<21,REG$Color,"black"),pt.bg=REG$Color,legend = REG$Population)
dev.off()
