#!/bin/Rscript

setwd("/Users/pierrespc/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F")

######


BSP<-read.table("../Maternal/BSP_DataTable_N24",stringsAsFactors = F,header=T,skip=1)
BSP$TIME=BSP$Time/25
hapNE<-read.table("../HapNe/Lab_with_Compendium.1240K/TH50000/Migrant/HapNe/hapne.csv",stringsAsFactors = F,header=T,sep=",")
hapNE$Time<-(hapNE$TIME-1)*25


pdf("DemoCurves/DemoCurves.pdf",width=10)
plot(0,0,xlim=c(0,min(c(max(BSP$TIME),max(hapNE$TIME)))),ylim=range(c(100,hapNE$Q0.025,hapNE$Q0.975,BSP$Lower,BSP$Upper)),log="y",axes=F,ann=F)
title(main="Effective Population Size Changes Across Time",
      xlab="Generations",ylab="Ne")
axis(1)
#axis(2,at=seq(1000,50000,by=las=2)
axis(2, at=c(seq(200,900,100),
             seq(2000,9000,1000),
             seq(20000,90000,10000),
             seq(200000,900000,100000)), tcl= -0.2,labels=NA)
axis(2, at=c(100,1000,10000,100000), tcl= -0.5,labels=c("1e2","1e3","1e4","1e5"),las=2)
box()


polygonX=c(BSP$TIME,BSP$TIME[c(nrow(BSP):1)])
polygonY=c(BSP$Lower,BSP$Upper[c(nrow(BSP):1)])
polygon(x=polygonX,y=polygonY,col = "palegreen1",border=NA)
points(BSP$TIME,BSP$Median,type="l",col="green4",lwd=2)

text(x=BSP$TIME[1],BSP$Median[1],labels = paste(round(BSP$Median[1]),"\n[",round(BSP$Lower[1]),"-",round(BSP$Upper[1]),"]",sep=""),cex=0.5,vjust=1)

polygonX=c(hapNE$TIME,hapNE$TIME[c(nrow(hapNE):1)])
polygonY=c(hapNE$Q0.025,hapNE$Q0.975[c(nrow(hapNE):1)])
polygon(x=polygonX,y=polygonY,col = "lightblue1",border=NA)
points(hapNE$TIME,hapNE$Q0.5,type="l",col="blue3",lwd=2)

text(x=hapNE$TIME[1],hapNE$Q0.5[1],labels = paste(round(hapNE$Q0.5[1]),"\n[",round(hapNE$Q0.025[1]),"-",round(hapNE$Q0.975[1]),"]",sep=""),cex=0.5,vjust=1)

legend("bottomright",
       lty=1,
       lwd=c(2,6,2,6),
       col=c("green4","palegreen1","blue3","lightblue1"),
       legend=c("Median from Bayesian Skyline Plot (Mitogenomes)",
                "95% Confidence Interval",
                "Median from hapNe (Autosomal)",
                "95% Confidence Interval"))
dev.off()
       