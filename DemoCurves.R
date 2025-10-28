#!/bin/Rscript

setwd("/Users/pierrespc/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F")

######


BSP<-read.table("../Maternal/BSP_DataTable_N24",stringsAsFactors = F,header=T,skip=1)
BSP$TIME=BSP$Time/25
hapNE<-read.table("../HapNe/Lab_with_Compendium.1240K/TH50000/Migrant/HapNe/hapne.csv",stringsAsFactors = F,header=T,sep=",")
hapNE$Time<-(hapNE$TIME-1)*25


svg("DemoCurves/DemoCurves.svg",width=14,height=6)
par(mar=c(5,5,4,1)+0.1)
#plot(0,0,xlim=c(min(c(max(BSP$TIME),max(hapNE$TIME))),0),ylim=range(c(100,hapNE$Q0.025,hapNE$Q0.975,BSP$Lower,BSP$Upper)),log="y",axes=F,ann=F)
plot(0,0,xlim=c(90,0),ylim=range(c(100,hapNE$Q0.025,hapNE$Q0.975,BSP$Lower,BSP$Upper)),log="y",axes=F,ann=F)
title(#main="Effective Population Size Changes Across Time",
      xlab="Generations",ylab="Ne",cex.lab=2)



polygonX=c(BSP$TIME,BSP$TIME[c(nrow(BSP):1)])
polygonY=c(BSP$Lower,BSP$Upper[c(nrow(BSP):1)])
polygon(x=polygonX,y=polygonY,col = "palegreen1",border=NA)
points(BSP$TIME,BSP$Median,type="l",col="green4",lwd=2)



polygonX=c(hapNE$TIME,hapNE$TIME[c(nrow(hapNE):1)])
polygonY=c(hapNE$Q0.025,hapNE$Q0.975[c(nrow(hapNE):1)])
polygon(x=polygonX,y=polygonY,col = "lightblue1",border=NA)
points(hapNE$TIME,hapNE$Q0.5,type="l",col="blue3",lwd=2)
abline(h=c(BSP$Median[1],hapNE$Q0.5[1]),col=c("green4","blue3") ,lty=2,lwd=1.5)

#abline(h=10^seq(1,5,1) ,lty=2,lwd=1)
#legend("bottomleft",
#       lty=1,
#       lwd=c(2,6,2,6),
#       col=c("green4","palegreen1","blue3","lightblue1"),
#       legend=c("Median from Bayesian Skyline Plot (Mitogenomes)",
#                "95% Confidence Interval",
#                "Median from hapNe (Autosomal)",
#                "95% Confidence Interval"),
#       bg="white")

axis(1,cex.axis=2)
#axis(2,at=seq(1000,50000,by=las=2)
axis(2, at=c(seq(200,900,100),
             seq(2000,9000,1000),
             seq(20000,90000,10000),
             seq(200000,900000,100000)), tcl= -0.2,labels=NA)
axis(2, at=c(100,1000,10000,100000), tcl= -0.5,labels=c("1e2","1e3","1e4","1e5"),las=2,cex.axis=2)
box()

text(x=BSP$TIME[1],BSP$Median[1],labels = paste(round(BSP$Median[1]),"\n[",round(BSP$Lower[1]),"-",round(BSP$Upper[1]),"]",sep=""),cex=2,adj = c(1,1),col="green4")
text(x=hapNE$TIME[1],hapNE$Q0.5[1],labels = paste(round(hapNE$Q0.5[1]),"\n[",round(hapNE$Q0.025[1]),"-",round(hapNE$Q0.975[1]),"]",sep=""),cex=2,adj = c(1,1),col="blue3")



dev.off()




svg("DemoCurves/DemoCurves_Veritcal.svg",height=10)
par(mar=c(5,7,0,0))
#plot(0,0,xlim=c(min(c(max(BSP$TIME),max(hapNE$TIME))),0),ylim=range(c(100,hapNE$Q0.025,hapNE$Q0.975,BSP$Lower,BSP$Upper)),log="y",axes=F,ann=F)
plot(0,0,ylim=c(70,-2),xlim=range(c(100,hapNE$Q0.025,hapNE$Q0.975,BSP$Lower,BSP$Upper)),log="x",axes=F,ann=F)
#title(main="Effective Population Size Changes Across Time",
#      ylab="Generations (Years Before Present)",
#      xlab="Ne")



polygonY=c(BSP$TIME,BSP$TIME[c(nrow(BSP):1)])
polygonX=c(BSP$Lower,BSP$Upper[c(nrow(BSP):1)])
polygon(x=polygonX,y=polygonY,col = "palegreen1",border=NA)
points(y=BSP$TIME,x=BSP$Median,type="l",col="green4",lwd=2)



polygonY=c(hapNE$TIME,hapNE$TIME[c(nrow(hapNE):1)])
polygonX=c(hapNE$Q0.025,hapNE$Q0.975[c(nrow(hapNE):1)])
polygon(x=polygonX,y=polygonY,col = "lightblue1",border=NA)
points(y=hapNE$TIME,x=hapNE$Q0.5,type="l",col="blue3",lwd=2)
#abline(v=c(BSP$Median[1],hapNE$Q0.5[1]),col=c("green4","blue3") ,lty=2,lwd=1)
segments(x0=BSP$Median[1], x1=BSP$Median[1], y0=-1,y1=90,col="green4",lty=2,lwd=1)
segments(x0=hapNE$Q0.5[1], x1=hapNE$Q0.5[1], y0=-1,y1=90,col="blue3",lty=2,lwd=1)

#abline(h=10^seq(1,5,1) ,lty=2,lwd=1)
#legend("bottom",
#       lty=1,
#      lwd=c(2,6,2,6),
#       col=c("green4","palegreen1","blue3","lightblue1"),
#       legend=c("Median from Bayesian Skyline Plot (Mitogenomes)",
#                "95% Confidence Interval",
#                "Median from hapNe (Autosomal)",
#                "95% Confidence Interval"),
#       bg="white")

axis(2,at=seq(0,60,10),labels = paste(seq(0,60,10),"\n(",seq(0,60,10)*25+650,")",sep=""),cex.axis=2,las=2)
#axis(2,at=seq(1000,50000,by=las=2)
axis(1, at=c(seq(200,900,100),
             seq(2000,9000,1000),
             seq(20000,90000,10000),
             seq(200000,900000,100000)), tcl= -0.2,labels=NA)
axis(1, at=c(100,1000,10000,100000), tcl= -0.5,labels=c("1e2","1e3","1e4","1e5"),las=1,cex.axis=2)
box()

#text(y=-3,x=BSP$Median[1],labels = paste(round(BSP$Median[1]),"\n[",round(BSP$Lower[1]),"-",round(BSP$Upper[1]),"]",sep=""),adj=c(0.5,1),cex=3,col="green4")
#text(y=-3,x=hapNE$Q0.5[1],labels = paste(round(hapNE$Q0.5[1]),"\n[",round(hapNE$Q0.025[1]),"-",round(hapNE$Q0.975[1]),"]",sep=""),cex=3,adj=c(0.5,1),col="blue3")

dev.off()

       