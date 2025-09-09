#!/bin/bash
require(stringr)
require(jpeg)

setwd("~/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F/")



dico<-data.frame(rbind(cbind("Site"="Barrancas B6",
                  "Abbr"="BB6",
                  "nAdna" = 6,
                  "nMob"=14,
                  "nDiet"=20,
                  "Color"="darkorange4",
                  "Point"=23),
            cbind("Site"= "Túmulo II",
                  "Abbr"="T-II",
                  "nAdna"=3,
                  "nMob"=9,
                  "nDiet"=15,
                  "Color"="darkorange4",
                  "Point"=22),
            cbind("Site"="Barrio Ramos I",
                  "Abbr"="BR-I",
                  "nAdna" = 2,
                  "nMob"=3,
                  "nDiet"=6,
                  "Color"="darkorange4",
                  "Point"=21),
            cbind("Site"="Potrero Las Colonias",
                  "Abbr"="PLC",
                  "nAdna"=35,
                  "nMob"=33,
                  "nDiet"=52,
                  "Color"="goldenrod1",
                  "Point"=21),
            cbind("Site"="Túmulo I",
                  "Abbr"="T-I",
                  "nAdna" = 0,
                  "nMob"=14,
                  "nDiet"=7,
                  "Color"="darkorange4",
                  "Point"=24),
            cbind("Site"="Túmulo III",
                  "Abbr"="T-III",
                  "nAdna"=0,
                  "nMob"=14,
                  "nDiet"=3,
                  "Color"="darkorange4",
                  "Point"=25),
            cbind("Site"="Usina Sur 2",
                  "Abbr"="US2",
                  "nAdna" = 0,
                  "nMob"=2,
                  "nDiet"=2,
                  "Color"="goldenrod1",
                  "Point"=21),
            cbind("Site"="Aguas de las Avispas",
                  "Abbr"="A.Av",
                  "nAdna" = 0,
                  "nMob"=1,
                  "nDiet"=1,
                  "Color"="goldenrod1",
                  "Point"=22)))


dates<-read.table("Maps/DatsFig1.tsv",stringsAsFactors = F,header=T,sep="\t")

MAIZ<-readJPEG("Maps/Maiz.jpg")
LLAMA<-readJPEG("Maps/Llama.jpg")
dico<-merge(dico,dates,by.x="Abbr",by.y="Site.Event")
dico$Point=as.numeric(dico$Point)
ref<-read.table("../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv",stringsAsFactors = F,header=T,sep="\t",comment.char = "@",quote="\"")

ref<-ref[ grepl("Conchal",ref$Population) | grepl("Rieles",ref$Population) | grepl("Aconcagua",ref$Population) | ref$Study=="delaFuente_2024_GBE" | ref$Study=="Pedersen_2021_MBE",  ]

ref$DateRange<-str_replace(ref$DateRange,"_","")
change<-function(string,split,pos){
  return(strsplit(string,split=split)[[1]][pos])
}
ref$end=sapply(ref$DateRange,change,split="-",pos=1)
ref$start=sapply(ref$DateRange,change,split="-",pos=2)
ref$start[is.na(ref$start)]<-ref$end[is.na(ref$start)]
ref$start<-as.numeric(ref$start)
ref$end<-as.numeric(ref$end)
dico$Start<-as.numeric(dico$Start)
dico$End<-as.numeric(dico$End)

ref$End=ref$end
ref$Start=ref$start
#ref$End[ ref$start>4000]=ref$start[ ref$start>4000]-ref$end[ ref$start>4000]+3000
#ref$Start[ ref$start>4000]=ref$start[ ref$start>4000]+3000

ref$End[ ref$start>10000]=ref$End[ ref$start>10000]-4000
ref$Start[ ref$start>10000]=ref$start[ ref$start>10000]-4000


#svg("Maps/TimeTransect.svg",height=12)
pdf("Maps/TimeTransect.pdf",height=12)
#plot(0,0,"n",xlim=c(-2,nrow(dico)+nrow(ref))+0.5,ylim=c(-3000,max(c(ref$Start,dates$Start))),ann=F,axes=F)
plot(0,0,"n",xlim=c(-2,nrow(dico)+nrow(ref))+0.5,ylim=c(-2700,9000),ann=F,axes=F)

rect(ybottom = dates$Start[dates$Site.Event=="First cultigens"],
     ytop=dates$End[dates$Site.Event=="First cultigens"],
     xleft=-100,xright=100,col="lightseagreen",border="lightseagreen")
text(x=-1,y=dates$Start[dates$Site.Event=="First cultigens"]+300,labels = "First\ncultigens",col="lightseagreen",font=3)
rect(ybottom = dates$Start[dates$Site.Event=="Inka arrival"],
     ytop=dates$End[dates$Site.Event=="Inka arrival"],
     xleft=-100,xright=100,col="lightseagreen",border="lightseagreen",lwd=3)
text(x=-1,y=dates$Start[dates$Site.Event=="Inka arrival"]+300,labels = "Inka\narrival",col="lightseagreen",font=3)
ref<-ref[ order(ref$start),]
dico<-dico[ order(dico$Start),]
axis(3,at=mean(c(-2,nrow(ref))),"Archaeological sites\nwith published ancient DNA data",tick = F)
segments(x0=nrow(ref)+0.5,x1=nrow(ref)+0.5,y0=-2650,y1=50000)
abline(h=0)
axis(3,at=mean(c(1,nrow(dico)))+nrow(ref),"Archaeological sites\nfrom this study",tick=F)


ref$Population<-sapply(ref$Population,change,split="_",2)
ref$Population<-str_remove(ref$Population,".Capt")
ref$Population[ref$Population=="LosRieles"]<-"Los Rieles"
ref$Population[ref$Population=="Conchali"]<-"Conchalí"
ref$Population[ref$Population=="Plomo"]<-"El Plomo"
ref$Population[ref$Population=="Toro"]<-"El Toro"
ref$Population[ref$Population=="Callingasta"]<-"Calingasta"
for(i in c(1:nrow(ref))){
  rect(xleft = i-0.3,xright=i+0.3,ytop=ref$Start[i],ybottom = ref$End[i],col=ref$Color[i],border = ref$Color[i],lwd=3)
  points(x=i,y=-200,pch=ref$Point[i],col=ifelse(ref$Point[i]<21,ref$Color[i],"black"),bg=ref$Color[i])
  text(x=i,y=-1900,labels = ref$N_Start[i])
  text(x=i,y=-2200,labels = "-")
  text(x=i,y=-2500,labels = "-")
  text(x=i,y=-1000,labels = ref$Population[i],col=ref$Color[i],srt=90)
}

x=nrow(ref)
for(i in c(1:nrow(dico))){
  x=x+1
  rect(xleft = x-0.3,xright=x+0.3,ytop=dico$Start[i],ybottom = dico$Start[i],col=dico$Color,border = dico$Color[i],lwd=3)
  points(x=x,y=-200,pch=dico$Point[i],col=ifelse(dico$Point[i]<21,dico$Color[i],"black"),bg=dico$Color[i])
  text(x=x,y=-1900,labels = dico$nAdna[i])
  text(x=x,y=-2200,labels = dico$nDiet[i])
  text(x=x,y=-2500,labels = dico$nMob[i])
  text(x=x,y=-1000,labels = dico$Abbr[i],col=dico$Color[i],srt=90)
  if(dico$Abbr[i]  %in% c("BB6","T-I")){
    rasterImage(LLAMA,x-0.4,-1650,x+0.4,-1300,col=dico$Color[i])
  }else{
    rasterImage(MAIZ,x-0.4,-1650,x+0.4,-1300,col=dico$Color[i])
  }
}
axis(2,at=c(seq(0,6000,1000),seq(7000,9000,1000)),labels = c(seq(0,6000,1000),seq(11000,13000,1000)),las=2)
#axis(2,at=c(0,-2250),labels = c("",""))

text(rep(-1,3),y=c(-1900,-2200,-2500),labels = c("aDNA","Diet","Mobility"))

segments(x0 =-3,x1=nrow(dico)+nrow(ref)+0.5,y0=-1750,y1=-1750)
segments(x0 =-3,x1=nrow(dico)+nrow(ref)+0.5,y0=-2050,y1=-2050)
segments(x0 =-3,x1=nrow(dico)+nrow(ref)+0.5,y0=-2350,y1=-2350)
segments(x0 =-3,x1=nrow(dico)+nrow(ref)+0.5,y0=-2650,y1=-2650)
title(ylab="Years Before Present")
dev.off()