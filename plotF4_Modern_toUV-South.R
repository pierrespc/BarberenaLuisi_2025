#!/bin/Rscript
require(stringr)
require(ggplot2)
require(ggbeeswarm)

setwd("~/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F//")
params=c("../F4_Modern_toSouthRefandUspallata/Lab_with_Compendium_GEHmodern_GenoModern.1240K//TH30000/F4.Modern_toSouthRefandUspallata.finalSet.Lab_with_Compendium_GEHmodern_GenoModern.1240K.TH30000.OUT",
         "../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv")
if(length(params)!=2){
	stop("<F4file> <annotref> ")
}
f4file<-params[1]
annotfile<-params[2]

outfile<-paste("F4_Modern_toSouthRefUspallata/F4.Modern_toSouthRefandUspallata.finalSet.Lab_with_Compendium_GEHmodern_GenoModern.1240K.TH30000.pdf")

annot<-read.table(annotfile,stringsAsFactors = F,header=T,sep="\t",comment.char="@")
annot<-annot[,c("Population","Region","Point","Color")]
f4<-read.table(f4file,stringsAsFactors = F,header=T,sep="\t")
listPop2<-unique(f4$Pop2_X)
f4<-merge(f4,annot,by.x="Pop2_X",by.y="Population")
if(sum(! listPop2%in% f4$Pop2_X)>0){
	print(listPop2[! listPop2 %in% f4$Pop2_X])
	stop()
}
names(f4)[names(f4) %in% names(annot)]<-paste(names(f4)[names(f4) %in% names(annot)],".Pop2_X",sep="")
listPop4<-unique(f4$Pop4_Z)
f4<-merge(f4,annot,by.x="Pop4_Z",by.y="Population")
if(sum(! listPop4%in% f4$Pop4_Z)>0){
  print(listPop4[! listPop4 %in% f4$Pop4_Z])
  stop()
}
names(f4)[names(f4) %in% names(annot)]<-paste(names(f4)[names(f4) %in% names(annot)],".Pop4_Z",sep="")

###remove no Southern Cone modernPop
f4<-f4[ ! f4$Region.Pop2_X %in% c("Amazonia_Modern","Amazonia_UnmaskedModern",
                                  "CentralAndes_Modern","CentralAndes_UnmaskedModern",
                                  "Peru_Modern",
                                  "SouthernNorthAmerica_Modern","SouthernNorthAmerica_UnmaskedModern"),]
###remove no South American  ancientPop
f4<-f4[ ! f4$Region.Pop4_Z %in% c("SouthernNorthAmerica_LH","SouthernNorthAmerica_EH"),]
medianRegionModern<-c()

for(reg in unique(f4$Region.Pop2_X)){
  medianRegionModern[reg]<-median(f4$Z[f4$Region.Pop2_X==reg])
}
medianRegionModern<-sort(medianRegionModern)
medianPopModern<-c()
for(pop in unique(f4$Pop2_X)){
  medianPopModern[pop]<-median(f4$Z[f4$Pop2_X==pop])
}

f4$Compa<-NA
x=0
xPrev=1
listLabels<-c()
listVline<-c()
for(reg in names(medianRegionModern)){
  tmp<-sort(medianPopModern[ names(medianPopModern) %in% f4$Pop2_X[ f4$Region.Pop2_X==reg] ])
  for(pop in names(tmp)){
    x=x+1
    f4$Compa[f4$Pop2_X==pop]<-paste(ifelse(x<10,paste("0",x,sep=""),x),": f4(Mbuti,",pop,";Uspallata, Ancient Group)",sep="")
  }
  listVline<-c(listVline,x+0.5)
  listLabels[str_replace(str_remove(reg,"_Modern"), "TierraDelFuego","Tierra del\nFuego")]<-mean(c(x,xPrev))
  xPrev=x+1
  
}

colorScale=c()
fillScale=c()
pointScale=c()
for(i in c(1:nrow(f4))){
  colorScale[f4$Pop4_Z[i]]=ifelse(f4$Point.Pop4_Z[i]<21,f4$Color.Pop4_Z[i],"black")
  fillScale[f4$Pop4_Z[i]]=f4$Color.Pop4_Z[i]
  pointScale[f4$Pop4_Z[i]]=f4$Point.Pop4_Z[i]
}


breaks=sort(c(3,-3,0,seq(round(min(f4$Z)),round(max(f4$Z)),5)))
pdf(outfile,width=10)
ggplot(f4, aes(x = Compa, y = Z)) +
  coord_cartesian(ylim=c(min(f4$Z),max(f4$Z)))+
  #geom_point(position = position_jitter(seed = 1, width = 0.2),shape=f4$Point,fill = f4$Color,color=ifelse(f4$Point<21,f4$Color,"black"),size=2,stroke=1.3) +
  geom_quasirandom(method = "pseudorandom",size=1.5,stroke=0.5,aes(fill=Pop4_Z,color=Pop4_Z,shape=Pop4_Z))+
  geom_violin(alpha = 0.5,fill=NA,color="black") +
  scale_fill_manual(values = fillScale)+
  scale_color_manual(values = colorScale)+
  scale_shape_manual(values = pointScale)+
  geom_hline(yintercept = c(-3,3),linetype = 2)+ 
  geom_hline(yintercept = 0,linetype = 1)+ 
  geom_vline(xintercept = listVline,linetype=2)+
  annotate(geom="text",x=listLabels,y=rep(max(f4$Z),length(listLabels)),label=names(listLabels),size=1.4)+
  theme_classic()+
  theme(legend.position="none",
        axis.text = element_text(size = 10),
        axis.text.x=element_text(size=6,angle=90))+
  scale_y_continuous(breaks=breaks)+
  labs(x="",y="Z")
dev.off()  
  

write.table(f4,str_replace(outfile,"pdf","tsv"),col.names=T,row.names=F,sep="\t",quote=F)
