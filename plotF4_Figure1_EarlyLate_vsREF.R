#!/bin/Rscript
require(stringr)
require(ggplot2)
require(ggbeeswarm)


setwd("/Users/pierrespc/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F/")


th=50000

for(type in c("SNA")){
  outTAB<-c()
  for(set in c("1240K","1240K.TVs","SG","SG.TVs")){
    
    
    params=c(paste("../F4_Late_toEarlyRefs/Lab_with_Compendium_GEHmodern.",set,"/TH",th,"/F4.Late_toEarlyRefs.finalSet.Lab_with_Compendium_GEHmodern.",set,".TH",th,".OUT",sep=""),
             paste("../F4_ref_toLateEarly/Lab_with_Compendium_GEHmodern.",set,"/TH",th,"/F4.Late_toEarlyRefs.finalSet.Lab_with_Compendium_GEHmodern.",set,".TH",th,".OUT",sep=""),
             "../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv",
             "../../Uspallata_Annotation.tsv",
             paste("F4_EarlyLate_vsRefs/F4_EarlyLate_vs",type,"Refs.",set,".pdf",sep=""))
    if(length(params)!=5){
      stop("<refPOPY> <refPOPX> <annotref> <annotstudy> <pdfout>")
    }
    refpopY<-params[1]
    refpopX<-params[2]
    annotfile<-params[3]
    annotStudy=params[4]
    outfile<-params[5]
    
    annot<-read.table(annotfile,stringsAsFactors = F,header=T,sep="\t",comment.char="@")
    annot<-annot[,c("Population","Region","Point","Color")]
    study<-read.table(annotStudy,stringsAsFactors = F,header=T,sep="\t",comment.char="@")
    study<-unique(study[,c("MainRegion","MainRegion","Point","Color")])
    study$Point[study$Color=="darkorange4"]<-21
    names(study)<-names(annot)
    study<-unique(study)
    
    #annot<-rbind(annot,study)
    annot<-rbind(annot,study)
    
    
    f4popy<-read.table(refpopY,stringsAsFactors = F,header=T,sep="\t")
    print(nrow(f4popy))
    f4popy<-merge(f4popy,annot,by.x="Pop3_Y",by.y="Population")
    print(nrow(f4popy))
    f4popy$Pop2_X[f4popy$Pop2_X=="Early"]<-"LH-HG"
    f4popy$Pop4_Z[f4popy$Pop4_Z=="Early"]<-"LH-HG"
    f4popy$Pop2_X[f4popy$Pop2_X=="Late"]<-"LH-F"
    f4popy$Pop4_Z[f4popy$Pop4_Z=="Late"]<-"LH-F"
    f4popy$Compa<-paste("F4(Mbuti,",f4popy$Pop2_X,";\nancient SA group;\n",f4popy$Pop4_Z,")",sep="")
    f4popy$REF=f4popy$Pop3_Y
    
    
    f4popx<-read.table(refpopX,stringsAsFactors = F,header=T,sep="\t")
    print(nrow(f4popx))
    f4popx<-merge(f4popx,annot,by.x="Pop2_X",by.y="Population")
    print(nrow(f4popx))
    f4popx$Pop3_Y[f4popx$Pop3_Y=="Early"]<-"LH-HG"
    f4popx$Pop4_Z[f4popx$Pop4_Z=="Early"]<-"LH-HG"
    f4popx$Pop3_Y[f4popx$Pop3_Y=="Late"]<-"LH-F"
    f4popx$Pop4_Z[f4popx$Pop4_Z=="Late"]<-"LH-F"
    f4popx$Compa<-paste("F4(Mbuti,ancient ",type,";\n",f4popx$Pop3_Y,",",f4popx$Pop4_Z,")",sep="")
    f4popx$REF=f4popx$Pop2_X
    
    
    f4<-rbind(f4popx[,c("REF","Compa","Color","Point","Region","Dstat","stderr","Z","BABA","ABBA","NSNPs")],f4popy[,c("REF","Compa","Color","Point","Region","Dstat","stderr","Z","BABA","ABBA","NSNPs")])
    
    if(grepl("SG",set)){
      f4<-f4[ ! str_ends(f4$REF,".Capt"),]
    }
    
    pdf(outfile)
    if(type=="SouthAmerica"){
      f4<-f4[  (grepl("Patagonia",f4$Region) | grepl("CentralAndes",f4$Region) | 
                  grepl("Brazil",f4$Region) | grepl("Pampa",f4$Region) |
                  grepl("Uruguay",f4$Region) | grepl("QhapaqHucha_Inka_LH",f4$Region) |
                  grepl("CentralChile",f4$Region) | grepl("Cuyo",f4$Region)) ,]
    }
    colorScale=c()
    fillScale=c()
    pointScale=c()
    for(i in c(1:nrow(f4))){
      colorScale[f4$REF[i]]=ifelse(f4$Point[i]<21,f4$Color[i],"black")
      fillScale[f4$REF[i]]=f4$Color[i]
      pointScale[f4$REF[i]]=f4$Point[i]
    }
    
    breaks=seq(round(min(c(-3,f4$Z))),
               round(max(c(3,f4$Z))),1)
    
    
    print(ggplot(f4, aes(x = Compa, y = Z)) +
            geom_violin(alpha = 0.5,fill=NA,color="black",scale = "width") +
            coord_cartesian(ylim=c(min(c(-3,f4$Z)),max(c(3,f4$Z))))+
            #geom_point(position = position_jitter(seed = 1, width = 0.2),shape=f4$Point,fill = f4$Color,color=ifelse(f4$Point<21,f4$Color,"black"),size=2,stroke=1.3) +
            geom_quasirandom(method = "pseudorandom",size=2,stroke=1,aes(fill=REF,color=REF,shape=REF))+
            scale_fill_manual(values = fillScale)+
            scale_color_manual(values = colorScale)+
            scale_shape_manual(values = pointScale)+
            geom_hline(yintercept = c(-3,3),linetype = 2)+
            geom_hline(yintercept=0,linetype = 1)+
            theme_classic()+
            theme(legend.position="none",
                  axis.text = element_text(size = 13))+
            scale_y_continuous(breaks=breaks)+
            labs(x="",y="Z"))
    
    dev.off()  
    
    names(f4)[1]<-paste(type,"Group",sep="")
    names(f4)[length(f4)-c(5:0)]<-paste(names(f4)[length(f4)-c(5:0)],"_",set,sep="")
    if(length(outTAB)==0){
      outTAB<-f4
    }else{
      outTAB<-merge(outTAB,f4,by=names(outTAB)[c(1:5)],all=T)
    }
  }
  
  outTAB$Compa<-str_replace_all(outTAB$Compa,"\n","")
  write.table(outTAB,paste("F4_EarlyLate_vsRefs/F4_EarlyLate_vs",type,"Refs.tsv",sep=""),col.names=T,row.names=F,sep="\t",quote=F)                
  
}
