#!/bin/Rscript
require(ggplot2)
require(stringr)
require(ggbeeswarm)

setwd("~/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F/")
orderRegions<-read.table("NeHapROH_vsCondHE/listRegionOrdered.txt",stringsAsFactors = F,header=F)$V1
focusplot="points"
THcond="50000"


listPOP<-c()

for(setCond in c("SG","1240K")[2]){
  
  
  annot<-read.table("../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv",stringsAsFactors = F,header=T,sep="\t",comment.char = "@",quote="\"")
  annotUsp<-read.table("../../Uspallata_Annotation.tsv",stringsAsFactors = F,header=T,comment.char = "@",sep="\t")
  annotUsp$MainRegion<-ifelse(annotUsp$MainRegion=="Migrant","LH-MF",
                              ifelse(annotUsp$MainRegion=="Local_Late","LH-LF",
                                     ifelse(annotUsp$MainRegion=="Local_Early","LH-HG",
                                            ifelse(annotUsp$MainRegion=="Local_PLC","LH-LFplc",
                                                   ifelse(annotUsp$MainRegion=="Migrant_Outlier","LH-MFout","WTF")))))
  annotUsp<-unique(annotUsp[,c("MainRegion","MainRegion","Color","Point")])
  names(annotUsp)[c(1:2)]<-c("Region","Population")
  annotUsp<-annotUsp[ ! (annotUsp$Region=="LH-LF" & annotUsp$Point %in% c(22,24)),]
  
  annot<-annot[   grepl("Brazil",annot$Region) | 
                    grepl("California",annot$Region) |
                    grepl("CentralAmerica",annot$Region) | 
                    grepl("CentralAndes",annot$Region) | 
                    grepl("CentralChile",annot$Region) | 
                    grepl("Cuyo_LH" ,annot$Region) | 
                    grepl("Pampa",annot$Region) | 
                    grepl("Uruguay_LH",annot$Region) | 
                    grepl("SouthPatagonia",annot$Region),]
  
  annot<-rbind(annot[,c("Region","Population","Color","Point")],annotUsp)
  
  if(setCond=="SG"){
    annot<-annot[ ! str_ends(annot$Population,".Capt"),]
    
  }
  annotByReg<-c()
  for(i in unique(annot$Region)){
    annotByReg<-rbind(annotByReg,annot[annot$Region==i,][1,])
  }
  condHE_perReg<-read.table(paste("../ConditionalHeterozygosity_RefByRegion/Lab_with_Compendium.",setCond,"/CondHet.",setCond,".TH",THcond,".ConditionalHeterozygosity.tsv",sep=""),
                            stringsAsFactors = F,header=T)
  condHE_perPop<-read.table(paste("../ConditionalHeterozygosity/Lab_with_Compendium.",setCond,"/CondHet.",setCond,".TH",THcond,".ConditionalHeterozygosity.tsv",sep=""),
                            stringsAsFactors = F,header=T)
  
  condHE_perPop$Population<-ifelse(condHE_perPop$Population=="Migrant","LH-MF",
                                   ifelse(condHE_perPop$Population=="Local_Late","LH-LF",
                                          ifelse(condHE_perPop$Population=="Local_Early","LH-HG",
                                                 ifelse(condHE_perPop$Population=="Local_PLC","LH-LFplc",
                                                        ifelse(condHE_perPop$Population=="Migrant_Outlier","LH-MFout",condHE_perPop$Population)))))
  
  condHE_perReg$Population<-ifelse(condHE_perReg$Population=="Migrant","LH-MF",
                                   ifelse(condHE_perReg$Population=="Local_Late","LH-LF",
                                          ifelse(condHE_perReg$Population=="Local_Early","LH-HG",
                                                 ifelse(condHE_perReg$Population=="Local_PLC","LH-LFplc",
                                                        ifelse(condHE_perReg$Population=="Migrant_Outlier","LH-MFout",condHE_perReg$Population)))))
  NEroh_perReg<-read.csv("../hapROH/NeEstimate_TH400000/RefByRegion_Ne_estimates.tsv",stringsAsFactors = F,header=T,sep="\t")
  NEroh_perReg<-NEroh_perReg[! is.na(NEroh_perReg$std.err),]
  
  NEroh_perPop<-read.csv("../hapROH/NeEstimate_TH400000/RefByPop_Ne_estimates.tsv",stringsAsFactors = F,header=T,sep="\t")
  NEroh_perPop<-NEroh_perPop[! is.na(NEroh_perPop$std.err),]
  
  
  NEroh_Usp<-read.csv("../hapROH/NeEstimate_TH400000/Ne_estimates.tsv",stringsAsFactors = F,header=T,sep="\t")
  NEroh_perReg<-rbind(NEroh_perReg,NEroh_Usp)
  NEroh_perPop<-rbind(NEroh_perPop,NEroh_Usp)
  
  NEroh_perReg$Population<-ifelse(NEroh_perReg$Population=="Migrant","LH-MF",
                                  ifelse(NEroh_perReg$Population=="Local_Late","LH-LF",
                                         ifelse(NEroh_perReg$Population=="Local_Early","LH-HG",
                                                ifelse(NEroh_perReg$Population=="Local_PLC","LH-LFplc",
                                                       ifelse(NEroh_perReg$Population=="Migrant_Outlier","LH-MFout",NEroh_perReg$Population)))))
  
  NEroh_perPop$Population<-ifelse(NEroh_perPop$Population=="Migrant","LH-MF",
                                  ifelse(NEroh_perPop$Population=="Local_Late","LH-LF",
                                         ifelse(NEroh_perPop$Population=="Local_Early","LH-HG",
                                                ifelse(NEroh_perPop$Population=="Local_PLC","LH-LFplc",
                                                       ifelse(NEroh_perPop$Population=="Migrant_Outlier","LH-MFout",NEroh_perPop$Population)))))
  
  
  
  pdf(paste("NeHapROH_vsCondHE/NeHapROH_vsCondHE.",setCond,".TH",THcond,".pdf",sep=""))
  
  for(refs in c("Population","Region")[1]){
    print(refs)
    if(refs=="Population"){
      out<-merge(condHE_perPop,NEroh_perPop,by="Population")
      out_REF<-merge(out,annot,by="Population")
      print("regions removed")
      print(sort(unique(out$Region[ !out$Population %in% out_REF$Population])))
      
    }else{
      out<-merge(condHE_perReg,NEroh_perReg,by="Population")
      out_REF<-merge(out,annotByReg,by.x="Population",by.y="Region")
      print("regions removed")
      print(sort(unique(out$Population[ !out$Population %in% out_REF$Population])))
      
    }
    
    
    print(nrow(out_REF))
    ##get Uspallata entries
    Migrants=out[ grepl("LH-MF",out$Population) ,]
    Local_Early=out[ grepl("LH-HG",out$Population) ,]
    Local_Late=out[ grepl("LH-LF",out$Population),]
    
    ###remove very low confidence Ne estimations
    out_REF<-out_REF[ out_REF$std.err <200,]  
    print(nrow(out_REF))
    out_REF<-out_REF[! out_REF$Population %in% c("LH-MF","LH-HG","LF-LF"),]
    
    
    ###get range for plot
    rangeX=c(min(out_REF$D-2*out_REF$SE),
             max(out_REF$D+2*out_REF$SE))
    rangeY=c(min(out_REF$Ne-2*out_REF$std.err),
             max(out_REF$Ne+2*out_REF$std.err))
    
    plot(out_REF$D,out_REF$Ne,"n",xlim=rangeX,ylim=rangeY,
         xlab="Conditionnal Heterozygosity",
         ylab="Ne from hapROH",
         main=paste("Genomic insights into effective population size\nPublished individuals grouped by",refs))
    
    ###add Migrants 
    if(focusplot =="points"){
      segments(x0=Migrants$D,x1=Migrants$D,
               y0=Migrants$Ne-2*Migrants$std.err,
               y1=Migrants$Ne+2*Migrants$std.err,
               col="goldenrod1",lwd=2)
      segments(y0=Migrants$Ne,y1=Migrants$Ne,
               x0=Migrants$D-2*Migrants$SE,
               x1=Migrants$D+2*Migrants$SE,
               col="goldenrod1",lwd=2)
      points(x=Migrants$D,
             y=Migrants$Ne,pch=21, col="black", bg= "goldenrod1",cex=3)
      
    }else{
      rect(xleft=Migrants$D-2*Migrants$SE,
           xright=Migrants$D+2*Migrants$SE,
           ybottom=Migrants$Ne-2*Migrants$std.err,
           ytop=Migrants$Ne+2*Migrants$std.err,
           border="black",col="goldenrod1")
    }
    
    ###add Local_Early 
    if(focusplot =="points"){
      segments(x0=Local_Early$D,x1=Local_Early$D,
               y0=Local_Early$Ne-2*Local_Early$std.err,
               y1=Local_Early$Ne+2*Local_Early$std.err,
               col="darkorange4",lwd=2)
      segments(y0=Local_Early$Ne,y1=Local_Early$Ne,
               x0=Local_Early$D-2*Local_Early$SE,
               x1=Local_Early$D+2*Local_Early$SE,
               col="darkorange4",lwd=2)
      points(x=Local_Early$D,
             y=Local_Early$Ne,pch=23, col="black", bg= "darkorange4",cex=1.4)
    }else{
      rect(xleft=Local_Early$D-2*Local_Early$SE,
           xright=Local_Early$D+2*Local_Early$SE,
           ybottom=Local_Early$Ne-2*Local_Early$std.err,
           ytop=Local_Early$Ne+2*Local_Early$std.err,
           border="black",col="darkorange4")
    }
    ###add Local_Late 
    if(focusplot =="points"){
      segments(x0=Local_Late$D,x1=Local_Late$D,
               y0=Local_Late$Ne-2*Local_Late$std.err,
               y1=Local_Late$Ne+2*Local_Late$std.err,
               col="darkorange4",lwd=2)
      segments(y0=Local_Late$Ne,y1=Local_Late$Ne,
               x0=Local_Late$D-2*Local_Late$SE,
               x1=Local_Late$D+2*Local_Late$SE,
               col="darkorange4",lwd=2)
      
      points(x=Local_Late$D,
             y=Local_Late$Ne,pch=21, col="black", bg= "darkorange4",cex=3)
    }else{
      rect(xleft=Local_Late$D-2*Local_Late$SE,
           xright=Local_Late$D+2*Local_Late$SE,
           ybottom=Local_Late$Ne-2*Local_Late$std.err,
           ytop=Local_Late$Ne+2*Local_Late$std.err,
           border="black",col="darkorange4")
    }
    ###IC for NE
    segments(x0=out_REF$D,x1=out_REF$D,
             y0=out_REF$Ne-2*out_REF$std.err,
             y1=out_REF$Ne+2*out_REF$std.err,
             col=out_REF$Color)
    ###IC for Het
    segments(y0=out_REF$Ne,y1=out_REF$Ne,
             x0=out_REF$D-2*out_REF$SE,
             x1=out_REF$D+2*out_REF$SE,
             col=out_REF$Color)
    
    points(out_REF$D,out_REF$Ne,
           pch=out_REF$Point,
           col=ifelse(out_REF$Point<21,out_REF$Color,"black"),
           bg=out_REF$Color)
    
    
    if(refs == "Population"){
      leg<-unique(out_REF[,c("Population","Color","Point","Region")])
      leg<-leg[ order(leg$Region),]
      leg$Population=paste(leg$Region,":",leg$Population)
      print(leg[order(leg$Population),c("Population")])
    }else{
      leg<-unique(out_REF[,c("Population","Color","Point")])
      leg<-leg[ order(leg$Population),]
    }
    plot(0,0,"n",axes=F,ann=F)
    legend("top",pch=leg$Point,col=ifelse(leg$Point<21,leg$Color,"black"),pt.bg=leg$Color,legend = leg$Population,cex=0.5)
    
    #legend("bottom",pch=c(21,21,23),pt.bg=c("goldenrod1","darkorange4","darkorange4"), col="black",
    #                                          legend = c("Migrants","Local_Late","Local_Early"))
    legend("bottom",pch=c(21,21,23),pt.bg=c("goldenrod1","darkorange4"), col="black",
           legend = c("LH-MF","LH-LF"))
    
    
    listPOP<-unique(c(listPOP,condHE_perPop$Population,NEroh_perPop$Population))
  }
  
  
  dev.off()
  
  
  
  
  #######NOW ONLY He
  write.table(NEroh_perPop$Population[ NEroh_perPop$Population %in% annot$Group],paste("NeHapROH_vsCondHE/Only.",setCond,".TH",THcond,".ListPops.tsv",sep=""),col.names=F,row.names=F,quote=F)
  
  #pdf(paste("NeHapROH_vsCondHE/OnlyHE.",setCond,".TH",THcond,".pdf",sep=""),width=10)
  svg(paste("NeHapROH_vsCondHE/OnlyHE.",setCond,".TH",THcond,".svg",sep=""),width=10)
  
  for(refs in c("Population","Region")[1]){
    if(refs=="Population"){
      out_REF<-merge(condHE_perPop,annot,by="Population")
      print("populations removed")
      print(sort(unique(condHE_perPop$Population[ !condHE_perPop$Population %in% out_REF$Population])))
      
    }else{
      out_REF<-merge(condHE_perReg,annotByReg[,!names(annotByReg)=="Population"],by.x="Population",by.y="Region")
      out_REF$Region=out_REF$Population
      print("regions removed")
      print(sort(unique(condHE_perReg$Population[ !condHE_perReg$Population %in% out_REF$Population])))
      
    }
    
    nrow(out_REF)
    colorScale=c()
    fillScale=c()
    pointScale=c()
    sizeScale=c()
    for(i in c(1:nrow(out_REF))){
      colorScale[out_REF$Population[i]]=ifelse(out_REF$Point[i]<21,out_REF$Color[i],"black")
      fillScale[out_REF$Population[i]]=out_REF$Color[i]
      pointScale[out_REF$Population[i]]=out_REF$Point[i]
      sizeScale[out_REF$Population[i]]=ifelse(out_REF$Population[i] %in% annotUsp$Population,5,2)
    }
    #out_REF<-out_REF[order(out_REF$D),]
    print(out_REF$Region[ ! out_REF$Region %in% orderRegions])
    outTMP<-c()
    for(reg in orderRegions[orderRegions %in% out_REF$Region]){
      tmp<-out_REF[ out_REF$Region==reg,]
      tmp$mean=mean(tmp$D)
      outTMP<-rbind(outTMP,tmp[order(tmp$D),])
    }
    breaks=c()
    breaksColor=c()
    out_REF<-outTMP[order(outTMP$mean),]
    for(reg in orderRegions[orderRegions %in% out_REF$Region]){
      breaks[reg]<-mean(which(out_REF$Region==reg))
      breaksColor[reg]=out_REF$Color[out_REF$Region==reg][1]
    }
    
    out_REF$x=c(1:nrow(out_REF))
    
    out_REF$lower=out_REF$D-2*out_REF$SE
    out_REF$upper=out_REF$D+2*out_REF$SE
    print(ggplot(data = out_REF,aes(x=x,y=D)) +
            #geom_point(position = position_jitter(seed = 1, width = 0.2),shape=f4$Point,fill = f4$Color,color=ifelse(f4$Point<21,f4$Color,"black"),size=2,stroke=1.3) +
            geom_errorbar(aes(ymin = lower, ymax = upper,color=Population), width = 0.01)+
            
            geom_point(aes(fill=Population,color=Population,shape=Population,size=Population),
                       stroke=1)+
            scale_fill_manual(values = fillScale)+
            scale_color_manual(values = colorScale)+
            scale_shape_manual(values = pointScale)+
            scale_size_manual(values = sizeScale)+
            theme_classic()+
            theme(legend.position="none",
                  axis.text = element_text(size = 10),
                  axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5,color=breaksColor))+
            labs(x="",y="Conditional Heterozygosity")+
            scale_x_continuous(breaks=breaks,name = names(breaks)))
    
    write.table(out_REF$Population,paste("NeHapROH_vsCondHE/Only.",setCond,".TH",THcond,".ListPops.tsv",sep=""),col.names=F,row.names=F,quote=F)
    write.table(out_REF,paste("NeHapROH_vsCondHE/Only.",setCond,".TH",THcond,".tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  }
  dev.off()
  
}


#######NOW ONLY NE

#pdf(paste("NeHapROH_vsCondHE/OnlyNefromHapROH.pdf",sep=""),width=10)
svg(paste("NeHapROH_vsCondHE/OnlyNefromHapROH.svg",sep=""),width=10)

for(refs in c("Population","Region")[1]){
  if(refs=="Population"){
    out_REF<-merge(NEroh_perPop,annot,by="Population")
    print("populations removed")
    print(sort(unique(condHE_perPop$Population[ !condHE_perPop$Population %in% out_REF$Population])))
    
  }else{
    out_REF<-merge(NEroh_perReg,annotByReg[,!names(annotByReg)=="Population"],by.x="Population",by.y="Region")
    out_REF$Region=out_REF$Population
    print("regions removed")
    print(sort(unique(condHE_perReg$Population[ !condHE_perReg$Population %in% out_REF$Population])))
    
  }
  out_REF<-out_REF[ out_REF$std.err <200,]  
  
  nrow(out_REF)
  colorScale=c()
  fillScale=c()
  pointScale=c()
  sizeScale=c()
  for(i in c(1:nrow(out_REF))){
    colorScale[out_REF$Population[i]]=ifelse(out_REF$Point[i]<21,out_REF$Color[i],"black")
    fillScale[out_REF$Population[i]]=out_REF$Color[i]
    pointScale[out_REF$Population[i]]=out_REF$Point[i]
    sizeScale[out_REF$Population[i]]=ifelse(out_REF$Population[i] %in% annotUsp$Population,5,2)
  }
  #out_REF<-out_REF[order(out_REF$D),]
  print(out_REF$Region[ ! out_REF$Region %in% orderRegions])
  outTMP<-c()
  for(reg in orderRegions[orderRegions %in% out_REF$Region]){
    tmp<-out_REF[ out_REF$Region==reg,]
    tmp$mean=mean(tmp$Ne)
    outTMP<-rbind(outTMP,tmp[order(tmp$Ne),])
  }
  breaks=c()
  breaksColor=c()
  out_REF<-outTMP[order(outTMP$mean),]
  for(reg in orderRegions[orderRegions %in% out_REF$Region]){
    breaks[reg]<-mean(which(out_REF$Region==reg))
    breaksColor[reg]=out_REF$Color[out_REF$Region==reg][1]
  }
  
  out_REF$x=c(1:nrow(out_REF))
  
  out_REF$lower=out_REF$Ne-2*out_REF$std.err
  out_REF$upper=out_REF$Ne+2*out_REF$std.err
  print(ggplot(data = out_REF,aes(x=x,y=Ne)) +
          #geom_point(position = position_jitter(seed = 1, width = 0.2),shape=f4$Point,fill = f4$Color,color=ifelse(f4$Point<21,f4$Color,"black"),size=2,stroke=1.3) +
          geom_errorbar(aes(ymin = lower, ymax = upper,color=Population), width = 0.01)+
          
          geom_point(aes(fill=Population,color=Population,shape=Population,size=Population),
                     stroke=1)+
          scale_fill_manual(values = fillScale)+
          scale_color_manual(values = colorScale)+
          scale_shape_manual(values = pointScale)+
          scale_size_manual(values = sizeScale)+
          theme_classic()+
          theme(legend.position="none",
                axis.text = element_text(size = 10),
                axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5,color=breaksColor))+
          labs(x="",y="Ne")+
          scale_x_continuous(breaks=breaks,name = names(breaks)))
  write.table(out_REF$Population,"NeHapROH_vsCondHE/OnlyNefromHapROH.ListPops.tsv",col.names=F,row.names=F,quote=F)
  write.table(out_REF,"NeHapROH_vsCondHE/OnlyNefromHapROH.tsv",col.names=T,row.names=F,quote=F,sep="\t")
}
dev.off()




annot<-annot[ annot$Population %in% listPOP,]

listREG<-read.table("NeHapROH_vsCondHE/listRegionOrdered.txt",stringsAsFactors = F,header=F)$V1

annot$cex=0.5
forLeg<-c()
for(reg in listREG[ listREG %in% annot$Region & ! (grepl("LH-HG",listREG) | grepl("LH-LF",listREG) | grepl("LH-MF",listREG))]){
  forLeg<-rbind(forLeg,cbind("Population"=reg,"Point"=NA,"Color"=NA,"cex"=0.8))
  for(pop in annot$Population[ annot$Region==reg]){
    forLeg<-rbind(forLeg,annot[ annot$Population==pop,c("Population","Point","Color","cex")])
  }
}

forLeg<-rbind(forLeg,cbind("Population"="Uspallata","Point"=NA,"Color"=NA,"cex"=1.5))
forLeg<-rbind(forLeg,annot[ (grepl("LH-HG",listREG) | grepl("LH-LF",listREG) | grepl("LH-MF",listREG)),c("Population","Point","Color","cex")])
forLeg$Population<-ifelse(forLeg$Population=="Migrant","LH-MF",
                          ifelse(forLeg$Population=="Migrant_outlier","LH-MFout",
                                 ifelse(forLeg$Population=="Local_Early","LH-HG",
                                        ifelse(forLeg$Population=="Local_Late","LH-LF",
                                               ifelse(forLeg$Population=="Local_PLC","LH-LFplc",forLeg$Population)))))

forLeg$Point<-as.numeric(forLeg$Point)
forLeg$cex<-as.numeric(forLeg$cex)
pdf("NeHapROH_vsCondHE/Legend.pdf",height=15,width=15)
par(mar=c(0,0,0,0))
plot(0,0,"n",ann=F,axes=F)

legend("center",col=ifelse(forLeg$Point<21,forLeg$Color,"black"),pch=forLeg$Point,pt.bg=forLeg$Color,text.col = forLeg$Color,cex=forLeg$cex,legend = forLeg$Population,ncol=4)

dev.off()

write.table(forLeg,"NeHapROH_vsCondHE/Legend.tsv",col.names=T,row.names = F,sep="\t",quote=F)
