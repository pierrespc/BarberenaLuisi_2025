#!/bin/Rscript
require(stringr)
require(ggplot2)
require(ggbeeswarm)


setwd("/Users/pierrespc/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F///")


th=50000
type="SNA"
outTAB<-c()
for(set in c("1240K","1240K.TVs","SG","SG.TVs")){
  
  
  params=c(paste("../F4_ref_toPairwiseUspallata/Lab_with_Compendium_GEHmodern.",set,"/TH",th,"/F4.Ref_toPairwiseUspallata.finalSet.Lab_with_Compendium_GEHmodern.",set,".TH",th,".OUT",sep=""),
           "../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv",
           "../../Uspallata_Annotation.tsv")
  if(length(params)!=3){
    stop("<f4> <annotref> <annotstudy>")
  }
  f4file<-params[1]
  annotfile<-params[2]
  annotStudy=params[3]
  outfile<-paste("F4_ref_toPairwiseUspallata/F4.Ref_toPairwiseUspallata.finalSet.Lab_with_Compendium_GEHmodern.",set,".TH",th,".pdf",sep="")
  
  annot<-read.table(annotfile,stringsAsFactors = F,header=T,sep="\t",comment.char="@")
  annot<-annot[,c("Population","Region","Point","Color")]
  study<-read.table(annotStudy,stringsAsFactors = F,header=T,sep="\t",comment.char="@")
  study<-unique(study[,c("MainRegion","MainRegion","Point","Color")])
  names(study)<-names(annot)
  study<-study[ ! duplicated(study$Population),]
  study$Region="Uspallata"
  study<-unique(study)
  
  #annot<-rbind(annot,study)
  annot<-rbind(annot,study)
  
  
  f4<-read.table(f4file,stringsAsFactors = F,header=T,sep="\t")
  
  f4$Pop3_Y[f4$Pop3_Y=="Local_Early"]<-"LH-HG"
  f4$Pop3_Y[f4$Pop3_Y=="Local_Late"]<-"LH-LF"
  f4$Pop3_Y[f4$Pop3_Y=="Local_PLC"]<-"LH-LFplc"
  f4$Pop3_Y[f4$Pop3_Y=="Migrant"]<-"LH-MF"
  f4$Pop3_Y[f4$Pop3_Y=="Migrant_Outlier"]<-"LH-MFout"
  
  f4$Pop4_Z[f4$Pop4_Z=="Local_Early"]<-"LH-HG"
  f4$Pop4_Z[f4$Pop4_Z=="Local_Late"]<-"LH-LF"
  f4$Pop4_Z[f4$Pop4_Z=="Local_PLC"]<-"LH-LFplc"
  f4$Pop4_Z[f4$Pop4_Z=="Migrant"]<-"LH-MF"
  f4$Pop4_Z[f4$Pop4_Z=="Migrant_Outlier"]<-"LH-MFout"
  
  f4$Compa<-paste("F4(Mbuti,ancient SA group;\n",f4$Pop3_Y,",",f4$Pop4_Z,")",sep="")
  print(nrow(f4))
  f4<-merge(f4,annot,by.x="Pop2_X",by.y="Population")
  print(nrow(f4))
  
  
  if(grepl("SG",set)){
    f4<-f4[ ! str_ends(f4$Pop2_X,".Capt"),]
  }
  
  pdf(outfile)
  if(type=="SouthAmerica"){
    f4<-f4[  (grepl("Patagonia",f4$Region) | grepl("CentralAndes",f4$Region) | 
                grepl("Brazil",f4$Region) | grepl("Pampa",f4$Region) |
                grepl("Uruguay",f4$Region) | grepl("QhapaqHucha_Inka_LH",f4$Region) |
                grepl("CentralChile",f4$Region) | grepl("Cuyo",f4$Region) | grepl("Uspallata",f4$Region)) ,]
  }
  colorScale=c()
  fillScale=c()
  pointScale=c()
  sizeScale=c()
  for(i in c(1:nrow(f4))){
    colorScale[f4$Pop2_X[i]]=ifelse(f4$Point[i]<21,f4$Color[i],"black")
    fillScale[f4$Pop2_X[i]]=f4$Color[i]
    pointScale[f4$Pop2_X[i]]=f4$Point[i]
    sizeScale[ f4$Pop2_X[i]]=ifelse(f4$Region[i]=="Uspallata",4,2)
  }
  
  breaks=seq(round(min(c(-3,f4$Z))),
             round(max(c(3,f4$Z))),1)
  
  
  print(ggplot(f4, aes(x = Compa, y = Z)) +
          geom_violin(alpha = 0.5,fill=NA,color="black",scale = "width") +
          coord_cartesian(ylim=c(min(c(-3,f4$Z)),max(c(3,f4$Z))))+
          #geom_point(position = position_jitter(seed = 1, width = 0.2),shape=f4$Point,fill = f4$Color,color=ifelse(f4$Point<21,f4$Color,"black"),size=2,stroke=1.3) +
          geom_quasirandom(method = "pseudorandom",stroke=1,aes(fill=Pop2_X,color=Pop2_X,shape=Pop2_X,size=Pop2_X))+
          scale_fill_manual(values = fillScale)+
          scale_color_manual(values = colorScale)+
          scale_shape_manual(values = pointScale)+
          scale_size_manual(values = sizeScale)+
          geom_hline(yintercept = c(-3,3),linetype = 2)+
          geom_hline(yintercept=0,linetype = 1)+
          theme_classic()+
          theme(legend.position="none",
                axis.text = element_text(size = 13),
                axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))+
          scale_y_continuous(breaks=breaks)+
          labs(x="",y="Z"))
  
  dev.off()  
  
  
  names(f4)[c(5:10)]<-paste(names(f4)[c(5:10)],"_",set,sep="")
  if(length(outTAB)==0){
    outTAB<-f4
  }else{
    outTAB<-merge(outTAB,f4,by=names(f4)[-c(5:10)],all=T)
  }
}

outTAB$Compa<-str_replace_all(outTAB$Compa,"\n","")
write.table(outTAB,str_remove(str_replace(outfile,".pdf",".tsv"),set),col.names=T,row.names=F,sep="\t",quote=F)                

pdf("F4_ref_toPairwiseUspallata/Legend.pdf")
plot(0,0,"n",axes=F,ann=F)
study$abr<-ifelse(study$Population=="Migrant","LH-MF",
                  ifelse(study$Population=="Migrant_Outlier","LH-MFout",
                         ifelse(study$Population=="Local_Early","LH-HG",
                                ifelse(study$Population=="Local_Late","LH-LF",
                                       ifelse(study$Population=="Local_PLC","LH-LFplc","WTF")))))
  
legend("center",pch=study$Point,pt.bg=study$Color,col="black",
       legend=study$abr)
dev.off()
