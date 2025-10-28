
require(ggplot2)
require(scatterpie)
require(maps)
require(stringr)

setwd("~/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/AdmixtureMasked/WithoutProjection/Lab_with_Compendium_GEHmodern_GenoModern.1240K/BestRUNperK/")
world <- map_data('world')

K=12

if(K==15){
  listColors<-c("blue1","green4","goldenrod1","orange4","orange3","brown3")
}else{
  if(K==12){
    listColors<-c("blue1","palegreen1","goldenrod1","orange4","orange3","darkorange1")
  }
}

annot<-read.table("~/Documents/PostDocPasteur/aDNA/ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv",stringsAsFactors = F,header=T,sep="\t",comment.char = "@",quote="\"")

ind<-read.table(paste("ForAdmixture.TH30000.GENO0.5.MAF0.01.pruned.",K,".AncestryComponentByIndividual.txt",sep=""),stringsAsFactors = F,header=T)
group<-data.frame(matrix(NA,length(unique(ind$Population)),K))
row.names(group)<-unique(ind$Population)
names(group)<-names(ind)[c(1:K)+1]
groupMEAN<-group

for(comp in names(group)){
  median<-c()
  mean<-c()
  for(pop in unique(ind$Population)){
    median[pop]<-median(ind[ ind$Population==pop,comp])
    mean[pop]<-mean(ind[ ind$Population==pop,comp])
  }
  group[,comp]<-median
  groupMEAN[,comp]<-mean
}

pdf(paste("UVcomponentInModernCOSMO.",K,".pdf",sep=""),width=14,height=14)
par(mar=c(4,10,4,2)+0.1)
for(comp in listColors){
  c2<-list("Component"=comp,
           "Other"="grey90")
  
  
  groupMEAN<-group
  
  
  th=-1
  group<-group[ order(group[,comp],decreasing = T),]
  modern<-ind[ grepl("Luisi",ind$Population) | grepl("Homburger",ind$Population),]
  #modern$Color<-as.numeric(as.factor(modern$Color))
  modernSIGN<-unique(row.names(group[ group[,comp] >th, ]))
  modernSIGN<-modernSIGN[ modernSIGN %in% modern$Population]
  forPlot<-list()
  scatters<-c()
  
  iter=0
  for(pop in modernSIGN ){
    iter=iter+1
    forPlot[[paste(ifelse(iter<10,"0",""),iter,": ",pop,sep="")]]<-modern[modern$Population==pop,comp]
    
    scatters<-rbind(scatters,cbind(
      Label=paste(ifelse(iter<10,"0",""),iter,sep=""),
      Population=pop,
      x=annot$Longitude[annot$Population==pop],
      
      y=annot$Latitude[annot$Population==pop],
      N=sum(ind$Population==pop),
      radius=0.8,
      Component=group[pop,comp],
      Other=1-group[pop,comp]))
  }
  
  
  scatters<-data.frame(scatters,stringsAsFactors = F)
 
  for( col in c(3:length(scatters))){
    scatters[,col]<-as.numeric(scatters[,col])
  }
  
  
  scatters$y[scatters$Population=="AMBA-Luisi"]=scatters$y[scatters$Population=="AMBA-Luisi"]-1.5
  scatters$y[scatters$Population=="Cordoba-Luisi"]=scatters$y[scatters$Population=="Cordoba-Luisi"]+1.5
  scatters$x[scatters$Population=="SantaFe-Homburger"]=scatters$x[scatters$Population=="SantaFe-Homburger"]-1
  scatters$y[scatters$Population=="Hulliche-Reich"]=scatters$y[scatters$Population=="Hulliche-Reich"]+1.5
  scatters$x[scatters$Population=="Yaghan-Reich"]=scatters$x[scatters$Population=="Yaghan-Reich"]-2
  scatters$y[scatters$Population=="Trelew-Luisi"]=scatters$y[scatters$Population=="Trelew-Luisi"]-0.9
  scatters$x[scatters$Population=="Diaguita-Reich"]=scatters$x[scatters$Population=="Diaguita-Reich"]-0.9
  scatters$y[scatters$Population=="Salta-Luisi"]=scatters$y[scatters$Population=="Salta-Luisi"]-0.9
  scatters$x[scatters$Population=="Wichi-Reich"]=scatters$x[scatters$Population=="Wichi-Reich"]-0.9
  scatters$x[scatters$Population=="Maya2-Reich"]=scatters$x[scatters$Population=="Maya2-Reich"]-1.5
  scatters$x[scatters$Population=="Mayan_o2"]=scatters$x[scatters$Population=="Mayan_o2"]+1.5
  scatters$y[scatters$Population=="Mayan"]=scatters$y[scatters$Population=="Mayan"]+1.5
  scatters$x[scatters$Population=="Pima"]=scatters$x[scatters$Population=="Pima"]+1.5
  scatters$x[scatters$Population=="Surui"]=scatters$x[scatters$Population=="Surui"]+1.5
  scatters$x[scatters$Population=="Karitiana"]=scatters$x[scatters$Population=="Karitiana"]+1.5
  scatters$x[scatters$Population=="Mixtec"]=scatters$x[scatters$Population=="Mixtec"]-1.5
  
  scatters$y[scatters$Population=="Aymara"]=scatters$y[scatters$Population=="Aymara"]-1.5
  scatters$y[scatters$Population=="Piapoco"]=scatters$y[scatters$Population=="Piapoco"]-1.5
  scatters$y[scatters$Population=="Trujillo"]=scatters$y[scatters$Population=="Trujillo"]+1.5
  
  
  scatters[ scatters$x %in% scatters$x[duplicated(paste(scatters$x,scatters$y))],]
  
  
  scatters$xtext=scatters$x-1.1
  
  
  
  #svg("PIE_UVcomponentInModern_SouthernCone.15.svg",height=8,width=10)
  print(ggplot(world, aes(long, lat)) +
          geom_map(map=world, aes(map_id=region), fill=NA, color="grey70")+
          geom_scatterpie(aes(x=x, y=y,r=radius),
                          data=scatters, cols=names(c2), alpha=1,color="grey70",size=0.1)+
          scale_fill_manual(values = c2)+
          coord_cartesian(xlim=c(-78,-42),ylim=c(-57,-19.9))+
          #labs(title = "Genetic ancestry proportions (Admixture K= 8)\n") +
          theme_classic())
  
}
dev.off()




pdf(paste("UVcomponentInModernINDIGENOUS.",K,".pdf",sep=""),width=14,height=14)
par(mar=c(4,10,4,2)+0.1)
for(comp in listColors){
  c2<-list("Component"=comp,
           "Other"="grey90")
  
  
  groupMEAN<-group
  
  
  th=-1
  group<-group[ order(group[,comp],decreasing = T),]
  modern<-ind[ ! (grepl("Local",ind$Region) | grepl("Migrant",ind$Region) | grepl("_LH",ind$Region) | grepl("_Unknown",ind$Region) | grepl("_MH",ind$Region) | grepl("_EH",ind$Region) | grepl("Luisi",ind$Population) | grepl("Homburger",ind$Population)),]
  
  #modern$Color<-as.numeric(as.factor(modern$Color))
  modernSIGN<-unique(row.names(group[ group[,comp] >th, ]))
  modernSIGN<-modernSIGN[ modernSIGN %in% modern$Population]
  forPlot<-list()
  scatters<-c()
  
  iter=0
  for(pop in modernSIGN ){
    iter=iter+1
    forPlot[[paste(ifelse(iter<10,"0",""),iter,": ",pop,sep="")]]<-modern[modern$Population==pop,comp]
    
    scatters<-rbind(scatters,cbind(
      Label=paste(ifelse(iter<10,"0",""),iter,sep=""),
      Population=pop,
      x=annot$Longitude[annot$Population==pop],
      
      y=annot$Latitude[annot$Population==pop],
      N=sum(ind$Population==pop),
      radius=0.8,
      Component=group[pop,comp],
      Other=1-group[pop,comp]))
  }
  
  
  scatters<-data.frame(scatters,stringsAsFactors = F)
  
  for( col in c(3:length(scatters))){
    scatters[,col]<-as.numeric(scatters[,col])
  }
  
  
  scatters$y[scatters$Population=="AMBA-Luisi"]=scatters$y[scatters$Population=="AMBA-Luisi"]-1.5
  scatters$y[scatters$Population=="Cordoba-Luisi"]=scatters$y[scatters$Population=="Cordoba-Luisi"]+1.5
  scatters$x[scatters$Population=="SantaFe-Homburger"]=scatters$x[scatters$Population=="SantaFe-Homburger"]-1
  scatters$y[scatters$Population=="Hulliche-Reich"]=scatters$y[scatters$Population=="Hulliche-Reich"]+1.5
  scatters$x[scatters$Population=="Yaghan-Reich"]=scatters$x[scatters$Population=="Yaghan-Reich"]-2
  scatters$y[scatters$Population=="Trelew-Luisi"]=scatters$y[scatters$Population=="Trelew-Luisi"]-0.9
  scatters$x[scatters$Population=="Diaguita-Reich"]=scatters$x[scatters$Population=="Diaguita-Reich"]-0.9
  scatters$y[scatters$Population=="Salta-Luisi"]=scatters$y[scatters$Population=="Salta-Luisi"]-0.9
  scatters$x[scatters$Population=="Wichi-Reich"]=scatters$x[scatters$Population=="Wichi-Reich"]-0.9
  scatters$x[scatters$Population=="Maya2-Reich"]=scatters$x[scatters$Population=="Maya2-Reich"]-1.5
  scatters$x[scatters$Population=="Mayan_o2"]=scatters$x[scatters$Population=="Mayan_o2"]+1.5
  scatters$y[scatters$Population=="Mayan"]=scatters$y[scatters$Population=="Mayan"]+1.5
  scatters$x[scatters$Population=="Pima"]=scatters$x[scatters$Population=="Pima"]+1.5
  scatters$x[scatters$Population=="Surui"]=scatters$x[scatters$Population=="Surui"]+1.5
  scatters$x[scatters$Population=="Karitiana"]=scatters$x[scatters$Population=="Karitiana"]+1.5
  scatters$x[scatters$Population=="Mixtec"]=scatters$x[scatters$Population=="Mixtec"]-1.5
  
  scatters$y[scatters$Population=="Aymara"]=scatters$y[scatters$Population=="Aymara"]-1.5
  scatters$y[scatters$Population=="Piapoco"]=scatters$y[scatters$Population=="Piapoco"]-1.5
  scatters$y[scatters$Population=="Trujillo"]=scatters$y[scatters$Population=="Trujillo"]+1.5
  
  
  scatters[ scatters$x %in% scatters$x[duplicated(paste(scatters$x,scatters$y))],]
  
  
  scatters$xtext=scatters$x-1.1
  
  
  
  #svg("PIE_UVcomponentInModern_SouthernCone.15.svg",height=8,width=10)
  print(ggplot(world, aes(long, lat)) +
          geom_map(map=world, aes(map_id=region), fill=NA, color="grey70")+
          geom_scatterpie(aes(x=x, y=y,r=radius),
                          data=scatters, cols=names(c2), alpha=1,color="grey70",size=0.1)+
          scale_fill_manual(values = c2)+
          coord_cartesian(xlim=range(scatters$x),ylim=range(scatters$y))+
          #labs(title = "Genetic ancestry proportions (Admixture K= 8)\n") +
          theme_classic())
  
}
dev.off()



