#!/bin/bash
require(scales)
require(stringr)
setIND="Lab_with_Compendium_GEHmodern"
params<-commandArgs(trailingOnly=T)
setSNP=params[1]
dataset=paste(setIND,".",setSNP,sep="")
th="50000"
folder="/pasteur/helix/projects/Hotpaleo/pierre/Projects/2024-09-01_Uspallata_noHighCov//Analyses/qpWaves_Admixtools2_Uspallata/"
setwd(folder)


outfolder=paste(dataset,"/TH",th,sep="")
info<-read.table("../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv",stringsAsFactors = F,header=T,sep="\t")
listRightPops<-read.table("RightPops_Sets.tsv",stringsAsFactors=F,header=T)

listSETS<-names(listRightPops[-c(1:6)])
print(nrow(listRightPops))
info<-rbind(info[,c("Population","Region","Point","Color")],cbind("Population"=c("Late","Early"),"Region"="Uspallata","Point"=21,"Color"="goldenrod1"))
info$Point=as.numeric(info$Point)
listRightPops<-merge(listRightPops,info[,c("Population","Point","Color")],by.x="Pop",by.y="Population")

if(grepl("SG",setSNP)){
	listRightPops<-listRightPops[ ! str_ends(listRightPops$Pop,".Capt"),]
}

listOrder<-read.table("listRegionOrdered.txt")
tmp<-unique(listRightPops$Region[! listRightPops$Region %in% listOrder$V1])
if(length(tmp)>1){
  print("not in order")
  print(tmp)
  stop()
}

listRightPopsBU<-listRightPops
listRightPops<-c()
for(reg in listOrder$V1[c(length(listOrder$V1):1)]){
  tmp<-listRightPopsBU[listRightPopsBU$Region ==reg,]
  listRightPops<-rbind(listRightPops,tmp[order(tmp$Pop),])
}
print(nrow(listRightPops))
for(target in c("Uspallata")){
	numSet=0
	for(rightset in listSETS){
		print(c(target,rightset))
		tmp<-read.table(paste(outfolder,"/",target,"/qpWave_Outgroups",rightset,"/doublets/doublets.Outgroups",rightset,".",dataset,".TH",th,".tsv",sep=""),header=T,stringsAsFactors=F)
		names(tmp)[-1]<-paste(names(tmp)[-1],"_",rightset,sep="")
		if(rightset==listSETS[1]){
			out<-tmp
		}else{
			out<-merge(out,tmp,by="LeftPops")
		}
		print(nrow(out))
	}
	out<-merge(out,info[,c("Population","Region","Point","Color")],by.x="LeftPops",by.y="Population")
	print(nrow(out))
	tmp<-unique(out$Region[! out$Region %in% listOrder$V1])
	if(length(tmp)>1){
		print("not in order")
		print(tmp)
		stop()
	}
	outBU<-out
	out<-c()
	for(reg in listOrder$V1[c(length(listOrder$V1):1)]){
		tmp<-outBU[outBU$Region ==reg,]
		out<-rbind(out,tmp[order(tmp$LeftPop),])
	}
	write.table(out,paste("Summary_Doublets_",target,"_",length(listSETS),"Sets_",setSNP,".TH",th,".tsv",sep=""),sep="\t",quote=F,col.names=T,row.names=F)
	pdf(paste("Summary_Doublets_",target,"_",length(listSETS),"Sets_",setSNP,".TH",th,".pdf",sep=""),heigh=15,width=8)
	par(mar=c(1,20,4,1))
	plot(0,0,"n",xlim=c(0,(length(listSETS)+1)),ylim=c(0,(nrow(listRightPops)+nrow(out)+10)),
	     axes=F,ann=F)
	title(main=paste("qpWave results testing Rank 0 rejection with doublets\n(South American References & ",target,")\nand using ",length(listSETS)," sets for right populations\n",setSNP," all SNPs TH=",th,sep=""),cex.main=0.9)
	axis(3,at=c(1:length(listSETS)),paste("Set",c(1:length(listSETS))),tick = F,line=-3,cex.axis=0.9)
	y=0
	for(i in c(1:nrow(out))){
	  y=y+1
	  for(set in c(1:length(listSETS))){
	    rect(xleft = set-0.5,xright=set+0.5,
        	 ybottom = y-0.5,ytop=y+0.5,
	         border="black",
        	 col=ifelse(out[i,paste("Rank1_",listSETS[set],sep="")]<0.01,"darkblue",
                    ifelse(out[i,paste("Rank1_",listSETS[set],sep="")]<0.05,"blue",
                           ifelse(out[i,paste("Rank1_",listSETS[set],sep="")]<0.1,"lightblue","white"))))
        
	  }
	  points(x=0.4,y=y,bg=out$Color[i],pch=out$Point[i],col=ifelse(out$Point[i]<21,out$Color[i],"black"))
	  axis(2,at=y,labels = paste(out$Region[i],": ",out$LeftPops[i],sep=""),line=0,las=2,tick=F,cex.axis=0.6)
	}

	y=y+2
	text(x=c(1:length(listSETS)),y=y,labels = paste("Set",c(1:length(listSETS))),cex=0.9)

	y=y+4
	rect(xleft=c(0.5,1.75,3),xright=c(1.5,2.75,4),
	     ybottom=rep(y-0.75,3),ytop=rep(y+0.75,3),
	     col=c("darkblue","blue","lightblue"))
	text(x=c(1,2.25,3.5),y=rep(y,3),labels = paste("P <",c(0.01,0.05,0.1)),col=c("white","white","black"),cex=0.8)

	y=y+4
	for(i in c(1:nrow(listRightPops))){
	  y=y+1
	  for(set in c(1:length(listSETS))){
	    rect(xleft = set-0.5,xright=set+0.5,
	         ybottom = y-0.5,ytop=y+0.5,
	         border="black",
	         col=ifelse(listRightPops[i,listSETS[set]],"black","white"))
    
	  }
	  points(x=0.4,y=y,bg=listRightPops$Color[i],pch=listRightPops$Point[i],col=ifelse(listRightPops$Point[i]<21,listRightPops$Color[i],"black"))
	  axis(2,at=y,labels = paste(listRightPops$Region[i],": ",listRightPops$Pop[i],sep=""),line=0,las=2,tick=F,cex.axis=0.6)
	}
	dev.off()
}

