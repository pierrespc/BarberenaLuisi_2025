
library(rnaturalearth)
library(sf)
library(ggplot2)
library(scatterpie)
library(ggspatial)
library(ggrepel)
require(stringr)

setwd("~/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F/PieCharts/")

annot<-read.table("~/Documents/PostDocPasteur/aDNA/ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv",stringsAsFactors = F,header=T,sep="\t",comment.char = "@",quote="\"")


c2<-list("CSC"="goldenrod1",
         "Other"="grey95")
             
#pdf("UVcomponentInModern.pdf",width=12)
world <- ne_countries(scale = "medium", returnclass = "sf")

K=15

pdf("UVcomponentInModern.15.pdf",width=7,height=14)
par(mar=c(4,10,4,2)+0.1)
ind<-read.table(paste("../../AdmixtureMasked/WithoutProjection/Lab_with_Compendium_GEHmodern_GenoModern.1240K/BestRUNperK/ForAdmixture.TH30000.GENO0.5.MAF0.01.pruned.",K,".AncestryComponentByIndividual.txt",sep=""),stringsAsFactors = F,header=T)
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

th=-1
group<-group[ order(group$goldenrod1,decreasing = T),]
modern<-ind[ grepl("Modern",ind$Region),]
#modern$Color<-as.numeric(as.factor(modern$Color))
modernSIGN<-unique(row.names(group[ group$goldenrod1 >th, ]))
modernSIGN<-modernSIGN[ modernSIGN %in% modern$Population]
forPlot<-list()
scatters<-c()
col=c()
iter=0
for(pop in modernSIGN ){
  iter=iter+1
  forPlot[[paste(ifelse(iter<10,"0",""),iter,": ",pop,sep="")]]<-modern$goldenrod1[modern$Population==pop]
  col[paste(ifelse(iter<10,"0",""),iter,": ",pop,sep="")]=ifelse(unique(modern$Region[ modern$Population==pop]) %in% c("SanJuan_Modern","SanLuis_Modern","LaRioja_Modern","Mendoza_Modern"),"goldenrod1",
                                                                 ifelse(unique(modern$Region[ modern$Population==pop]) %in% c("Catamarca_Modern","Jujuy_Modern","Salta_Modern","Tucuman_Modern"),"goldenrod3",
                                                                        ifelse(unique(modern$Region[ modern$Population==pop]) %in% c("Cordoba_Modern"),"goldenrod2",
                                                                               ifelse(unique(modern$Region[ modern$Population==pop]) %in% c("AMBA_Modern","BsAs_Modern","LaPampa_Modern"),"goldenrod4",
                                                                                      ifelse(unique(modern$Region[ modern$Population==pop]) %in% c("CentralChile_Modern","Chile_Modern"),"brown1",
                                                                                             ifelse(unique(modern$Region[ modern$Population==pop]) %in% c("Bariloche_Modern","Chubut_Modern","SantaFe_Modern","TierraDelFuego_Modern","SouthPatagonia_Modern"),"brown3",
                                                                                                    "grey70"))))))
  
  scatters<-rbind(scatters,cbind(
    Label=paste(ifelse(iter<10,"0",""),iter,sep=""),
    Population=pop,
    x=annot$Longitude[annot$Population==pop],
    
    y=annot$Latitude[annot$Population==pop],
    N=sum(ind$Population==pop),
    radius=0.8,
    CSC=group[pop,"goldenrod1"],
    Other=1-group[pop,"goldenrod1"]))
}

boxplot(forPlot,las=2,cex.axis=0.7,main=paste("CSC proportion in modern individuals\nK = ",K,sep=""),
        col="goldenrod1",
        #border=col,
        horizontal = T,at = c(length(forPlot):1))

#legend("bottomright",pch=15,col=c("goldenrod1","goldenrod2","goldenrod3","goldenrod4",
#                               "brown1","brown3",
#                               "lightgrey"),
#       legend = c("Central Western Argentina", "Central Argentina","North Western Argentina","Central Eastern Argentina",
#                  "Central Chile","Patagonia","Others"))


dev.off()

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




pdf("PIE_UVcomponentInModern.15.pdf",height=14,width=14)
print(ggplot() +
  geom_sf(data = world, fill = NA, color = "grey20", size = 0.2)+
  scale_fill_manual(values = c2)+
  geom_scatterpie(aes(x=x, y=y,r=radius),
                  data=scatters, cols=names(c2), alpha=1,color="grey70",size=0.1)+
  geom_text(aes(x=xtext,y=y,label=Label),data=scatters,size=2.2,color="grey50")+
  coord_sf(xlim=range(scatters$x),ylim=range(scatters$y))+
  #labs(title = "Genetic ancestry proportions (Admixture K= 8)\n") +
  theme_classic())
dev.off()


pdf("PIE_UVcomponentInModern_SouthernCone.15.pdf",height=8,width=10)
#svg("PIE_UVcomponentInModern_SouthernCone.15.svg",height=8,width=10)
print(ggplot() +
        geom_sf(data = world, fill = NA, color = "grey20", size = 0.2)+
        scale_fill_manual(values = c2)+
        geom_scatterpie(aes(x=x, y=y,r=radius),
                        data=scatters, cols=names(c2), alpha=1,color="grey70",size=0.1)+
        #geom_text(aes(x=xtext,y=y,label=Label),data=scatters,size=2.2,color="grey50")+
        coord_sf(xlim=c(-78,-42),ylim=c(-57,-19.9))+
        #labs(title = "Genetic ancestry proportions (Admixture K= 8)\n") +
        theme_classic())

dev.off()

