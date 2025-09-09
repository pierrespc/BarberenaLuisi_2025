#!/bin/Rscript
require(stringr)
require(ggplot2)
require(ggbeeswarm)

setwd("~/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F//")
annotfile="../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv"
annotStudy="../../Uspallata_Annotation.tsv"
popUSP="Uspallata"
for(setsnp in c("1240K","1240K.TVs","SG","SG.TVs")[1]){
  for(TEMPORALITY in c("","_noLH")){
    print(setsnp)
    f4file=paste("../F4_NorthRef_toUspallataSouthRef/Lab_with_Compendium_GEHmodern.",setsnp,"/TH50000/F4.NorthRef_toUspallataSouthRef.finalSet.Lab_with_Compendium_GEHmodern.",setsnp,".TH50000.OUT",sep="")
    
    
    
    annot<-read.table(annotfile,stringsAsFactors = F,header=T,sep="\t",comment.char="@")
    annot<-annot[,c("Population","Region","Point","Color")]
    
    study<-data.frame(cbind("Population"=c("Uspallata"),
                            "Region"="Uspallata",
                            "Point"=21,
                            "Color"="goldenrod1"),stringsAsFactors=F)
    
    
    annot<-annot[ !annot$Population %in% c("Brazil_Alcobaca_4800-800BP.SG","Peru_SoroMikayaPatjxa_6800BP.SG"),]
    #annot<-rbind(annot,study)
    annot<-rbind(annot,study)
    f4<-read.table(f4file,stringsAsFactors = F,header=T,sep="\t")
    listPop4<-unique(f4$Pop4_Z)
    f4<-merge(f4,annot,by.x="Pop4_Z",by.y="Population")
    if(sum(! listPop4 %in% f4$Pop4_Z)>0){
      warning(listPop4[! listPop4 %in% f4$Pop4_Z])
      
    }
    
    listPop2<-unique(f4$Pop2_X)
    f4<-merge(f4,annot,by.x="Pop2_X",by.y="Population")
    if(sum(! listPop2 %in% f4$Pop2_X)>0){
      warning(listPop2[! listPop2 %in% f4$Pop2_X])
      
    }
    
    names(f4)<-str_replace(names(f4),".x",".South")
    names(f4)<-str_replace(names(f4),".y",".North")
    
    f4$Point.North<-as.numeric(f4$Point.North)
    f4$Point.South<-as.numeric(f4$Point.South)
    
    
    dim(f4)
    if(grepl("SG",f4file)){
      f4<-f4[ ! (str_ends(f4$Pop2_X,".Capt") | str_ends(f4$Pop4_Z,".Capt")),]
    }
    dim(f4)
    
    f4_ALL<-f4
    
    
    
    f4<-f4_ALL[ f4_ALL$Pop3_Y==popUSP,]
    
    outfile<-paste("F4_NorthRef_toUspallataSouthRef/F4.NorthRef_toUspallataSouthRef.finalSet.Lab_with_Compendium_GEHmodern.",setsnp,".TH50000.",popUSP,".summarized",TEMPORALITY,".pdf",sep="")
    
    if(TEMPORALITY!=""){
      f4<-f4[ ! (str_ends(f4$Region.South,"_LH") | str_ends(f4$Region.South,"_Unknown")),]
    }
    f4$Compa<-paste("F4(Mbuti,",f4$Region.North,";\n",popUSP,",Ancient South American group)",sep="")
    
    colorScale=c()
    fillScale=c()
    pointScale=c()
    for(i in c(1:nrow(f4))){
      colorScale[f4$Pop4_Z[i]]=ifelse(f4$Point.South[i]<21,f4$Color.South[i],"black")
      fillScale[f4$Pop4_Z[i]]=f4$Color.South[i]
      pointScale[f4$Pop4_Z[i]]=f4$Point.South[i]
    }
    
    #breaks=sort(c(3,seq(0,round(max(f4$Z)),5)))
    pdf(outfile,width=10)
    print(ggplot(f4, aes(x = Compa, y = Z)) +
            geom_hline(yintercept = 3,linetype = 2)+ 
            geom_hline(yintercept = -3,linetype = 2)+ 
            geom_hline(yintercept = 0,linetype = 1)+
            coord_cartesian(ylim=c(min(f4$Z),max(f4$Z)))+
            #geom_point(position = position_jitter(seed = 1, width = 0.2),shape=f4$Point,fill = f4$Color,color=ifelse(f4$Point<21,f4$Color,"black"),size=2,stroke=1.3) +
            geom_quasirandom(method = "pseudorandom",size=2,stroke=1,aes(fill=Pop4_Z,color=Pop4_Z,shape=Pop4_Z))+
            geom_violin(alpha = 0.5,fill=NA,color="black") +
            scale_fill_manual(values = fillScale)+
            scale_color_manual(values = colorScale)+
            scale_shape_manual(values = pointScale)+
            theme_classic()+
            theme(legend.position="none",
                  axis.text = element_text(size = 9),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            #scale_y_continuous(breaks=breaks)+
            labs(x="",y="Z"))
    dev.off()  
    
    if(TEMPORALITY==""){
      names(f4)[c(5:10)]<-paste(names(f4)[c(5:10)],".",setsnp,sep="")
      if(setsnp=="1240K"){
        outall<-f4
      }else{
        outall<-merge(outall,f4,by=names(f4)[-c(5:10)],all=T)
      }
      print(outall)
      
    }
  }
}
write.table(outall[,names(outall)!="Compa"],str_replace(outfile,"_noLH.pdf",".tsv"),col.names=T,row.names=F,sep="\t",quote=F)