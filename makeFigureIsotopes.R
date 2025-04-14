
require("ggplot2")
require(ggbeeswarm)
require(stringr)

setwd("/Users/pierrespc/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F/Figure2/")

####figure 3A

data<-read.table("Fig3a.txt",stringsAsFactors = F,header=T,sep="\t")

data$Color<-ifelse(data$Site %in% c("PLC","A. Avispas","Usina Sur 2") , "goldenrod1",
                   ifelse(data$Site %in% c("Túmulo I","Túmulo II","Túmulo III","BRI","BB6"),"darkorange4","WTF"))
                   
data$Point<-ifelse(data$Site=="PLC",21,
                   ifelse(data$Site=="A. Avispas",22,
                          ifelse(data$Site=="Usina Sur 2",23,
                                 ifelse(data$Site=="BB6",23,
                                        ifelse(data$Site=="BRI",21,
                                               ifelse(data$Site=="Túmulo II",22,
                                                      ifelse(data$Site=="Túmulo I",24,
                                                             ifelse(data$Site=="Túmulo III",25,"WTF"))))))))

data$Point<-as.numeric(data$Point)
#pdf("Fig3a.pdf",height=5,width=14)
svg("Fig3a.svg",height=5,width=14)
plot(data$Years.BP,data$X13Capat.,col="black",bg=data$Color,pch=data$Point,cex=1.5,xlim=range(data$Years.BP)[c(2,1)],
     xlab="Median Years BP",
     ylab="C13apa",ylim=c(-13,-2))
mod2 <- lm(X13Capat. ~ Years.BP + I(Years.BP^2), data=data)
mod3 <- lm(X13Capat. ~ Years.BP + I(Years.BP^2)+I(Years.BP^3), data=data)
LOESS<-loess(X13Capat. ~ Years.BP, data=data,span=0.5,degree=1)


values<-seq(range(data$Years.BP)[1],range(data$Years.BP)[2],10)
values<-data.frame(matrix(values,length(values),1))
names(values)<-"Years.BP"

pred2<-mod2$coefficients[1]+mod2$coefficients[2]*values+mod2$coefficients[3]*values^2
pred3<-mod3$coefficients[1]+mod3$coefficients[2]*values+mod3$coefficients[3]*values^2+mod3$coefficients[4]*values^3
#lines(pred2$x,pred2$fit, lty='solid', col='darkred', lwd=3)


predLOESS <- predict(LOESS, values$Years.BP, se=TRUE)
predMOD3<- predict.lm(mod3, values,se.fit=TRUE)
points(values$Years.BP, predLOESS$fit,lty='solid', col='darkorange', lwd=2,type='l')
points( values$Years.BP,predLOESS$fit-1.96*predLOESS$se.fit,lty='dashed', col='darkorange', lwd=1,type='l')
points(values$Years.BP,predLOESS$fit+1.96*predLOESS$se.fit, lty='dashed', col='darkorange', lwd=1,type='l')


#points(values$Years.BP,predMOD3$fit, lty='solid', col='darkgreen', lwd=3,type='l')
#points(values$Years.BP,predMOD3$fit-1.96*predMOD3$se.fit, lty='dashed', col='lightgreen', lwd=1,type='l')
#points(values$Years.BP,predMOD3$fit+1.96*predMOD3$se.fit, lty='dashed', col='lightgreen', lwd=1,type='l')

data<-data[ order(data$Years.BP),]
leg<-unique(data[,c("Site","Color","Point")])
leg$Site[leg$Site=="BRI"]<-"Barrio Ramos I"
leg$Site[leg$Site=="BB6"]<-"Barrancas B6"
leg$Site[leg$Site=="PLC"]<-"Potrero Las Colonias"
legend("topleft",pch=leg$Point,pt.bg=leg$Color,col="black",legend = leg$Site,cex=0.8)

#text(x=1800,y=-9.5,labels=paste("R2 = ",round(as.numeric(summary(mod2)["adj.r.squared"]),digits=4)),cex=0.6)
dev.off()



####figure 3B

data<-read.table("Fig3b.txt",stringsAsFactors = F,header=T,sep="\t")
names(data)[2]<-"Site"
names(data)[3]<-"Country"
data$Color<-ifelse(data$Site %in% c("Potrero Las Colonias") , "goldenrod1",
                   ifelse(data$Site %in% c("Tumulo I","Tumulo II","Tumulo III","Barrancas 6","Barrio Ramos I"),"darkorange4",
                          ifelse(data$Country%in% c("Peru","peru"),"orange1",
                                 ifelse(grepl("Chile",data$Country) | grepl("Atacama",data$Country),"brown3",
                                        ifelse(data$Country%in% c("Bolivia"),"sienna2",
                                               ifelse(grepl("Argentina",data$Country) | grepl("Mendoza",data$Country),"coral4","WTF"))))))




data$Point<-ifelse(data$Site %in% c("Potrero Las Colonias","Barrio Ramos I"),21,
                   ifelse(data$Site=="Barrancas 6",23,
                          ifelse(data$Site=="Tumulo II",22,
                                 ifelse(data$Site=="Tumulo I",24,
                                        ifelse(data$Site=="Tumulo III",25,"WTF")))))

for(otcol in c("orange1","brown3","sienna2","coral4")){
  data$Point[ data$Color==otcol]<-c(1:sum(data$Color==otcol))
}

data$Point<-as.numeric(data$Point)

data$es.13Capa<-data$ds.13Capa/sqrt(data$N)
data$es.13Ccol<-data$ds.13Ccol/sqrt(data$N)
data<-data[ order(data$Point),]
svg("Fig2b.svg")
plot(data$media.13Ccol,data$media.13Capa,"n",
     xlab="13Ccol",ylab="13Capa",
     xlim=range(c(data$media.13Ccol-1.96*data$es.13Ccol,data$media.13Ccol+1.96*data$es.13Ccol)),
     ylim=range(c(data$media.13Capa-1.96*data$es.13Capa,data$media.13Capa+1.96*data$es.13Capa)))
segments(x0 = data$media.13Ccol-1.96*data$es.13Ccol,
         x1 = data$media.13Ccol+1.96*data$es.13Ccol,
         y0 = data$media.13Capa,
         y1 = data$media.13Capa,col=data$Color)
segments(x0 = data$media.13Ccol,
         x1 = data$media.13Ccol,
         y0 = data$media.13Capa-1.96*data$es.13Capa,
         y1 = data$media.13Capa+1.96*data$es.13Capa,col=data$Color)
points(data$media.13Ccol,data$media.13Capa,pch=data$Point,col=ifelse(data$Point<21,data$Color,"black"),bg=data$Color,cex=ifelse(data$Point<21,1.2,1.8))


leg1<-data[ data$Point>20,c("Site","Color","Point")]
leg1<-leg1[ order(leg1$Color,leg1$Point),]
leg1$cex=0.65
leg1<-rbind(cbind("Site"="This study","Color"="black","Point"=NA,"cex"=0.8),leg1)
legTMP<-data[ data$Point<21,c("Site","Color","Point","Country")]
legTMP$cex=0.65
leg2<-c()
for(country in c("Peru","Bolivia","Atacama","Central Chile","North West Argentina","South Mendoza")){
  tmp<-legTMP[ legTMP$Country==country,]
  tmp<-tmp[order(tmp$Point),]
  leg2<-rbind(leg2,cbind("Site"=country,"Color"=unique(tmp$Color),"Point"=NA,"Country"=NA,"cex"=0.75),tmp)
}

leg2$Site[ leg2$Site=="SARNC Inland Late Intermediate"]="SARNC LIP"
#legend("topleft",pch=leg2$Point,col=leg2$Color,legend=paste(leg2$Site," (",leg2$X.3,")",sep=""),cex=0.7)
legend("bottomright",pch=as.numeric(leg2$Point),text.col = leg2$Color,col=leg2$Color,legend=leg2$Site,cex=as.numeric(leg2$cex),ncol=2)
legend("topleft",pch=as.numeric(leg1$Point),col="black",pt.bg=leg1$Color,legend=leg1$Site,cex=as.numeric(leg1$cex))
dev.off()






####Fig 4a

table<-read.csv("Fig4a.txt",sep="\t")
data<-matrix(NA,0,2)
i=0
for(n in names(table)){
  i=i+1
  data<-rbind(data,cbind(na.omit(table[,n]),paste(i,": ",str_replace_all(n,"\\."," "),"\nn = ",length(na.omit(table[,n])),sep="")))
}


data<-data.frame(data)
names(data)<-c("value","site")
data$site<-str_replace(data$site,"A  Avis","A. Avis")
data$site<-str_replace(data$site,"PLCpre ","PLC pre-")
data$value<-as.numeric(data$value)

fillScale=c("plum3",rep("goldenrod1",4),
             rep("darkorange4",4))

colorScale<-rep("black",9)

pointScale<-c(21,21,21,22,23,24,22,25,21)
                                        
names(fillScale)<-c("1: Bioavailable strontium\nn = 14","8: PLC migration\nn = 31","9: PLC pre-migration\nn = 19","6: A. Avispas\nn = 1","7: Usina Sur 2\nn = 2","2: Túmulo I\nn = 14","3: Túmulo II\nn = 14","4: Túmulo III\nn = 14","5: Barrio Ramos I\nn = 6")
names(colorScale)<-names(fillScale)
names(pointScale)<-names(fillScale)

svg("Fig4a.svg",width=17,heigh=6)
ggplot(data, aes(x = site, y = value,color=site,fill=site)) +
        #geom_violin(alpha = 0.5,fill=NA,color="black",scale = "width") +
        geom_boxplot(outliers = F)+
        coord_cartesian(xlim=c(1,(i+1)))+
        #geom_point(position = position_jitter(seed = 1, width = 0.2),shape=f4$Point,fill = f4$Color,color=ifelse(f4$Point<21,f4$Color,"black"),size=2,stroke=1.3) +
        geom_quasirandom(method = "pseudorandom",size=4,stroke=1,aes(fill=site,color=site,shape=site))+
        scale_fill_manual(values = fillScale)+
        scale_color_manual(values = colorScale)+
        scale_shape_manual(values = pointScale)+
        annotate("rect",
           xmin = -Inf, xmax = Inf, 
           ymin = min(data$value[ grepl("Bioa",data$site)]), ymax = max(data$value[ grepl("Bioa",data$site)]),  fill = "darkgrey", alpha=.3)+
        annotate("text",x=i+1,y=mean(c(min(data$value[ grepl("Bioa",data$site)]), max(data$value[ grepl("Bioa",data$site)]))),label="Range of local\nbioavailable Sr")+
        annotate("segment", x = i+1, y = mean(c(min(data$value[ grepl("Bioa",data$site)]), max(data$value[ grepl("Bioa",data$site)])))-0.0002, xend = i+1, yend = min(data$value[ grepl("Bioa",data$site)]),
           arrow = arrow(type = "closed", length = unit(0.02, "npc")),size=2)+
        annotate("segment", x = i+1, y = mean(c(min(data$value[ grepl("Bioa",data$site)]), max(data$value[ grepl("Bioa",data$site)])))+0.0002, xend = i+1, yend = max(data$value[ grepl("Bioa",data$site)]),
           arrow = arrow(type = "closed", length = unit(0.02, "npc")),size=2)+
        annotate("text",x=i+1,y=mean(c(min(data$value[ grepl("Bioa",data$site)]), min(data$value))),label="Non-local\nrange")+
        annotate("segment", x = i+1, y = mean(c(min(data$value[ grepl("Bioa",data$site)]), min(data$value)))+0.0002, xend = i+1, yend = min(data$value[ grepl("Bioa",data$site)]),
           arrow = arrow(type = "closed", length = unit(0.02, "npc")),size=2)+
        annotate("segment", x = i+1, y = mean(c(min(data$value[ grepl("Bioa",data$site)]), min(data$value)))-0.0002, xend = i+1, yend = min(data$value),
           arrow = arrow(type = "closed", length = unit(0.02, "npc")),size=2)+
        theme_classic()+
        theme(legend.position="none",axis.text = element_text(size = 13))+
        xlab("")+ylab("87Sr/86Sr")
dev.off()        
