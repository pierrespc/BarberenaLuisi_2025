#!/bin/Rscript

require(ade4)
require(adegenet)
#require(dartR.base)
#require(dartR)
require(poppr)
library(dartRverse)
require(stringr)

skipOnlyNotRestricted=T

setwd("~/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/FST/")

REF<-read.table("../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv",stringsAsFactors = F,header=T,sep="\t",quote="\"",comment.char="@")
REF<-REF[,c("Region","Population","Point","Color")]
Uspa<-read.table("../../Uspallata_Annotation.tsv",stringsAsFactors = F,header=T,sep="\t",quote="\"",comment.char="@")
Uspa<-unique(cbind("Region"="UV & CV",Uspa[,c("MainRegion","Point","Color")]))
names(Uspa)[2]<-"Population"
Uspa<-Uspa[ ! (Uspa$Point !=21 & Uspa$Population == "Local_Late"),]
REF$Region[REF$Region=="Cuyo_LH"]<-"UV & CV"
REF<-REF[order(REF$Region),]
ColorManifest<-rbind(Uspa,REF[REF$Region=="UV & CV",],REF[REF$Region!="UV & CV",])

#ColorManifest<-ColorManifest[ order(ColorManifest$Region),]
fam<-read.table("finalSet.1240K.TH50000.MAF0.01.plink.fam",stringsAsFactors = F,header=F)
counts<-as.data.frame(table(fam$V1),stringsAsFactors = F)
counts<-merge(counts,ColorManifest,by.x="Var1",by.y="Population")

listNotOne<-counts$Var1[ counts$Freq>1]
listOKCentralAndes<-counts$Var1[  ! (counts$Freq < 3 & grepl("CentralAndes",counts$Region))]
listPopRemove<-counts$Var1[grepl("_100BP",counts$Var1)]
counts<-counts[ counts$Var1 %in% c(listNotOne,listOKCentralAndes) & ! counts$Var1 %in% listPopRemove,]

fam<-merge(fam,ColorManifest,by.x="V1",by.y="Population")

countsRegion<-c()
for(Region in unique(counts$Region)){
  countsRegion<-rbind(countsRegion,cbind(Region,Freq=sum(counts$Freq[ counts$Region==Region]),Npop=sum(counts$Region==Region)))
}
countsRegion<-data.frame(countsRegion)
countsRegion$Freq<-as.numeric(countsRegion$Freq)
countsRegion$Npop<-as.numeric(countsRegion$Npop)

#genlight_object<-read.PLINK("finalSet.1240K.TH50000.MAF0.01.plink.raw",
#                            map.file = "finalSet.1240K.TH50000.MAF0.01.plink.map",
#                            quiet = FALSE,ploidy=2)

genlight_object<-gl.read.PLINK("finalSet.1240K.TH50000.MAF0.01.plink")

row.names(fam)<-fam$V2
fam<-fam[genlight_object@ind.names,]
names(fam)[1]<-"Group"
strata(genlight_object)<-fam[,c("Group","Region")]





amova_result<-list()
amova_test<-list()
amova_perm<-list()
listPops<-list()
amova_result_restricted<-list()
amova_test_restricted<-list()
amova_perm_restricted<-list()


###Among groups in UV
listPops_restricted<-list()
uv<-which(fam$Region=="UV & CV")
listPops[["UV & CV"]]<-unique(fam$Group[uv])
uvRestricted<-which(fam$Region=="UV & CV" & fam$Group %in% counts$Var1[ counts$Freq>1])
listPops_restricted[["UV & CV"]]<-unique(fam$Group[uvRestricted])
if(! skipOnlyNotRestricted){
  print("UV & CV")
  subGenlight<-genlight_object[uv]
  subGenlight<-gl.recalc.metrics(subGenlight)
  subGenlight<-gl.filter.allna(subGenlight)
  subGenlight<-gl.filter.maf(subGenlight,threshold = 0.01)

  amova.result <- poppr.amova(subGenlight, ~Group, missing = "loci",clonecorrect = TRUE)
  amova_result[["UV & CV"]]<-amova.result
  amova.test<-randtest(amova.result,nrepet = 99)
  amova_test[["UV & CV"]]<-amova.test

  amova_perm[["UV & CV"]]<-data.frame(matrix(NA,99,1))
  names(amova_perm[["UV & CV"]])<-"Phi-Group-total"

  GroupStrata<-strata(subGenlight)$Group
  names(amova_perm[["UV & CV"]])<-"Phi-Group-total"
  print(amova.result$statphi)
  for(perm in c(1:99)){
    print(perm)
    strata(subGenlight)$Group<-sample(GroupStrata)
    amova.result <- poppr.amova(subGenlight, ~Group, missing = "loci",clonecorrect = TRUE)
    amova_perm[["UV & CV"]][perm,"Phi-Group-total"]<-amova.result$statphi["Phi-Group-total","Phi"]
  }

  #print(summary(amova_perm[["UV & CV"]][,"Phi-Group-total"]))
}

print("UV & CV restricted")
subGenlight<-genlight_object[uvRestricted]
subGenlight<-gl.recalc.metrics(subGenlight)
subGenlight<-gl.filter.allna(subGenlight)
subGenlight<-gl.filter.maf(subGenlight,threshold = 0.01)
amova.result <- poppr.amova(subGenlight, ~Group, missing = "loci",clonecorrect = TRUE)
amova_result_restricted[["UV & CV"]]<-amova.result
amova.test<-randtest(amova.result,nrepet = 99)
amova_test_restricted[["UV & CV"]]<-amova.test
print(amova.result$statphi)


amova_perm_restricted[["UV & CV"]]<-data.frame(matrix(NA,99,1))
names(amova_perm_restricted[["UV & CV"]])<-"Phi-Group-total"
GroupStrata<-strata(subGenlight)$Group

print(amova.result$statphi)
for(perm in c(1:99)){
  print(c("UV & CV restricted",perm))
  strata(subGenlight)$Group<-sample(GroupStrata)
  amova.result <- poppr.amova(subGenlight, ~Group, missing = "loci",clonecorrect = TRUE)
  amova_perm_restricted[["UV & CV"]][perm,"Phi-Group-total"]<-amova.result$statphi["Phi-Group-total","Phi"]
}


###Among Region and groups in UV and one other region
#for(Region in countsRegion$Region[ grepl("Brazil",countsRegion$Region)| grepl("Patagonia",countsRegion$Region) | 
                                   #grepl("CentralAndes",countsRegion$Region) | grepl("Pampa",countsRegion$Region) |
                                   #grepl("Chile",countsRegion$Region)]){
for(Region in sort(countsRegion$Region[ ! grepl("Beringia",countsRegion$Region) &
                                   ! grepl("Carribean",countsRegion$Region) &
                                   ! grepl("China",countsRegion$Region) &
                                   ! grepl("Indonesia" ,countsRegion$Region) &
                                   ! grepl("Laos" ,countsRegion$Region) &
                                   ! grepl("NorthernNorthAmerica" ,countsRegion$Region) &
                                   ! grepl("Russia" ,countsRegion$Region) &
                                   ! grepl("Saami" ,countsRegion$Region)])){
  
  
  if(! skipOnlyNotRestricted){
    print(paste(Region))
    keep<-which(fam$Region ==Region )
    listPops[[Region]]<-unique(fam$Group[fam$Region ==Region])
    subGenlight<-genlight_object[keep]
    subGenlight<-gl.recalc.metrics(subGenlight)
    subGenlight<-gl.filter.allna(subGenlight)
    subGenlight<-gl.filter.maf(subGenlight,threshold = 0.01)
    amova.result <- poppr.amova(subGenlight, ~Group, missing = "loci",clonecorrect = TRUE)
    amova_result[[Region]]<-amova.result
    amova.test <- randtest(amova.result,nrepet = 99)
    amova_test[[Region]]<-amova.test
    print(amova.result$statphi)
  
    print(paste(Region,"with UV & CV"))
    keep<-c(uv,which(fam$Region ==Region ))
    subGenlight<-genlight_object[keep]
    subGenlight<-gl.recalc.metrics(subGenlight)
    subGenlight<-gl.filter.allna(subGenlight)
    subGenlight<-gl.filter.maf(subGenlight,threshold = 0.01)
    amova.result <- poppr.amova(subGenlight, ~Region/Group, missing = "loci",clonecorrect = TRUE)
    amova_result[[paste(Region,"with UV & CV")]]<-amova.result
    amova.test <- randtest(amova.result,nrepet = 99)
    amova_test[[paste(Region,"with UV & CV")]]<-amova.test
    print(amova.result$statphi)
  }
  
  ####ALL OF THIS IF WE WANT AT LEAST 2 pops in the region (if not no permutation allowed...)
  tmpCounts<-counts[ counts$Freq>1 & counts$Var1 %in% ColorManifest$Population[ColorManifest$Region==Region],]
  #tmpCounts<-counts[ counts$Var1 %in% ColorManifest$Population[ColorManifest$Region==Region],]
  if(nrow(tmpCounts)<3){
    print(c(Region,"less than 3 groups with N>1"))
    print(tmpCounts)
    next
  }
  if(sum(tmpCounts$Freq>2)<2){
    print(c(Region,"Only one group with N>2"))
    print(tmpCounts)
    next
  }
  
  
  listPops_restricted[[Region]]<-tmpCounts$Var1
  print(paste(Region,"only restricted"))
  if(Region %in% names(amova_perm_restricted)){
    print("already done")
    next
  }
  keep<-which(fam$Group %in% tmpCounts$Var1 )
  subGenlight<-genlight_object[keep]
  subGenlight<-gl.recalc.metrics(subGenlight)
  subGenlight<-gl.filter.allna(subGenlight)
  subGenlight<-gl.filter.maf(subGenlight,threshold = 0.01)
  amova.result <- poppr.amova(subGenlight, ~Group, missing = "loci",clonecorrect = TRUE)
  amova_result_restricted[[Region]]<-amova.result
  amova.test <- randtest(amova.result,nrepet = 99)
  amova_test_restricted[[Region]]<-amova.test
  print(amova.result$statphi)
  
  GroupStrata<-strata(subGenlight)$Group
  amova_perm_restricted[[Region]]<-data.frame(matrix(NA,99,1))
  names(amova_perm_restricted[[Region]])<-"Phi-Group-total"
  for(perm in c(1:99)){
    print(c(Region,"restricted only perm",perm))
    strata(subGenlight)$Group<-sample(GroupStrata)
    amova.result <- poppr.amova(subGenlight, ~Group, missing = "loci",clonecorrect = TRUE)
    amova_perm_restricted[[Region]][perm,"Phi-Group-total"]<-amova.result$statphi["Phi-Group-total","Phi"]
  }
  
  
  
  print(paste(Region,"with UV & CV restricted"))
  keep<-c(uvRestricted,which(fam$Group %in% tmpCounts$Var1 ))
  strata(subGenlight)$Group<-droplevels(strata(subGenlight)$Group)
  strata(subGenlight)$Region<-droplevels(strata(subGenlight)$Region)
  subGenlight<-genlight_object[keep]
  subGenlight<-gl.recalc.metrics(subGenlight)
  subGenlight<-gl.filter.allna(subGenlight)
  subGenlight<-gl.filter.maf(subGenlight,threshold = 0.01)
  amova.result <- poppr.amova(subGenlight, ~Region/Group, missing = "loci",clonecorrect = TRUE)
  amova_result_restricted[[paste(Region,"with UV & CV")]]<-amova.result
  amova.test <- randtest(amova.result,nrepet = 99)
  amova_test_restricted[[paste(Region,"with UV & CV")]]<-amova.test
  print(amova.result$statphi)
  GroupStrata<-data.frame("Region"=strata(subGenlight)$Region,"Group"=strata(subGenlight)$Group,"GroupPerm"=strata(subGenlight)$Group)
  RegionStrata<-unique(data.frame("Region"=strata(subGenlight)$Region,"Group"=strata(subGenlight)$Group))
  
  
  amova_perm_restricted[[paste(Region,"with UV & CV")]]<-data.frame(matrix(NA,99,2))
  names(amova_perm_restricted[[paste(Region,"with UV & CV")]])<-c("Phi-Group-Region","Phi-Region-total")
  strata(subGenlight)$GroupPerm<-strata(subGenlight)$Group
  
  for(perm in c(1:99)){
    print(c(Region,"restricted with UV & CV perm group",perm))
    GroupStrata$GroupPerm[ GroupStrata$Region=="UV & CV"]<-sample(GroupStrata$Group[ GroupStrata$Region=="UV & CV"])
    GroupStrata$GroupPerm[ GroupStrata$Region==Region]<-sample(GroupStrata$Group[ GroupStrata$Region==Region])
    ###first permnutation. for groups
    strata(subGenlight)$GroupPerm[ RegionStrata$Region=="UV & CV"]<-GroupStrata$GroupPerm[ RegionStrata$Region=="UV & CV"]
    strata(subGenlight)$GroupPerm[ RegionStrata$Region==Region]<-GroupStrata$GroupPerm[ RegionStrata$Region==Region]
    amova.result <- poppr.amova(subGenlight, ~Region/GroupPerm, missing = "loci",clonecorrect = TRUE)
    amova_perm_restricted[[paste(Region,"with UV & CV")]][perm,"Phi-Group-Region"]<-amova.result$statphi["Phi-GroupPerm-Region","Phi"]
    
    print(c(Region,"restricted with UV & CV perm region",perm))
    ###second permnutation. for regions
    RegionStrata$RegionPerm<-sample(RegionStrata$Region)
    tmpData<-strata(subGenlight)
    tmpData$RegionPerm<-as.character(tmpData$Region)
    for(pop in RegionStrata$Group){
      tmpData$RegionPerm[ tmpData$Group==pop]<-RegionStrata$RegionPerm[RegionStrata$Group==pop]
    }
    strata(subGenlight)<-tmpData
    amova.result <- poppr.amova(subGenlight, ~RegionPerm/Group, missing = "loci",clonecorrect = TRUE)
    amova_perm_restricted[[paste(Region,"with UV & CV")]][perm,"Phi-Region-total"]<-amova.result$statphi["Phi-RegionPerm-total","Phi"]
    
  }
  
}


AmongRegions<-c() #PhiCT
AmongGroupsWithinRegions<-c() #PhiSC
AmongGroupsAlone<-c() #PhiSC/PhiCT

for(Region in c("UV & CV",sort(names(listPops_restricted)[! grepl("UV & CV",names(listPops_restricted))]))){
  # amongroup alone
  AmongGroupsAlone<-rbind(AmongGroupsAlone,
                          cbind(Region,
                                Phi=amova_result_restricted[[Region]]$statphi["Phi-Group-total",1],
                                sd=sd(amova_perm_restricted[[Region]][,"Phi-Group-total"])))
}

maxPoints<--Inf
for(Region in sort(names(listPops_restricted)[ ! grepl("UV & CV",names(listPops_restricted))])){
  # amongroup alone

  regName=paste(Region,"with UV & CV")
  AmongRegions<-rbind(AmongRegions,
                          cbind(Region=Region,
                                Phi=amova_result_restricted[[regName]]$statphi["Phi-Region-total",1],
                                sd=sd(amova_perm_restricted[[regName]][,"Phi-Region-total"])))
  AmongGroupsWithinRegions<-rbind(AmongGroupsWithinRegions,
                      cbind(Region=Region,
                            Phi=amova_result_restricted[[regName]]$statphi["Phi-Group-Region",1],
                            sd=sd(amova_perm_restricted[[regName]][,"Phi-Group-Region"])))
  maxPoints<-max(c(maxPoints,length(listPops_restricted[[Region]])+length(listPops_restricted[["UV & CV"]])))
  
}
AmongGroupsWithinRegions<-as.data.frame(AmongGroupsWithinRegions,stringsAsFactors = F)
AmongGroupsAlone<-as.data.frame(AmongGroupsAlone,stringsAsFactors = F)
AmongRegions<-as.data.frame(AmongRegions,stringsAsFactors = F)

AmongRegions$Phi<-as.numeric(AmongRegions$Phi)
AmongRegions$sd<-as.numeric(AmongRegions$sd)
AmongGroupsAlone$Phi<-as.numeric(AmongGroupsAlone$Phi)
AmongGroupsAlone$sd<-as.numeric(AmongGroupsAlone$sd)
AmongGroupsWithinRegions$Phi<-as.numeric(AmongGroupsWithinRegions$Phi)
AmongGroupsWithinRegions$sd<-as.numeric(AmongGroupsWithinRegions$sd)




xAlone=c(1,seq(1,nrow(AmongGroupsWithinRegions),by=1)*3-1)
xGroupReg=seq(1,nrow(AmongGroupsWithinRegions),by=1)*3
xReg=seq(1,nrow(AmongGroupsWithinRegions),by=1)*3+1

NUMSD=3
VLIM=range(c( AmongGroupsAlone$Phi + rep(c(-1,1)*NUMSD,each=nrow(AmongGroupsAlone))*AmongGroupsAlone$sd,
              AmongGroupsWithinRegions$Phi + rep(c(-1,1)*NUMSD,each=nrow(AmongGroupsWithinRegions))*AmongGroupsWithinRegions$sd,
              AmongRegions$Phi + rep(c(-1,1)*NUMSD,each=nrow(AmongRegions))*AmongRegions$sd) )
VPHI=VLIM[2]
VLIM[2]<-VLIM[2]+(VLIM[2]-VLIM[1])*0.4

yPoints=seq(VLIM[2],VLIM[2]*0.8,length.out=maxPoints)


pdf("AMOVA.pdf",width=12)
plot(0,0,"n",axes=F,ann=F,
     xlim=c(1,(nrow(AmongGroupsAlone)+nrow(AmongRegions)+nrow(AmongGroupsWithinRegions)+1)),
     ylim=VLIM)
axis(2,at=VLIM[2]*0.9,labels = "Groups\nincluded",las=2,tick = F,cex.axis=0.7)

axis(2,at=seq(-0.05,0.15,0.025)[ seq(-0.05,0.15,0.025)<VPHI],cex.axis=0.7)


abline(h=c(0,VPHI*1.05))
title(ylab=paste("Phi +/-",NUMSD," permutation-based SD"))

rect(xleft = 0.5,xright=xAlone[length(xAlone)]+2.5,
     ybottom = AmongGroupsAlone$Phi[AmongGroupsAlone$Region=="UV & CV"]+NUMSD*AmongGroupsAlone$sd[AmongGroupsAlone$Region=="UV & CV"],
     ytop= AmongGroupsAlone$Phi[AmongGroupsAlone$Region=="UV & CV"]-NUMSD*AmongGroupsAlone$sd[AmongGroupsAlone$Region=="UV & CV"],
     col="goldenrod1",border=NA)

segments(x0 = 0.5,x1=xAlone[length(xAlone)]+2.5,
        y0=AmongGroupsAlone$Phi[AmongGroupsAlone$Region=="UV & CV"],lty=1,col="goldenrod3",lwd=2)



points(xAlone,AmongGroupsAlone$Phi,pch=16)
axis(1,at=xAlone,labels=rep("Groups in region",length(xAlone)),las=2,tick=T,cex.axis=0.5)
axis(1,at=xGroupReg,labels=rep("Groups w/in regions",length(xGroupReg)),las=2,tick=T,cex.axis=0.5)
axis(1,at=xReg,labels=rep("w/ UV&CV",length(xReg)),las=2,tick=T,cex.axis=0.5)
axis(3,at=c(1,xGroupReg),labels=str_replace_all(c("UV & CV",AmongGroupsWithinRegions$Region),"_","\n"),tick = F,cex.axis=0.5)
segments(x0=xAlone,x1=xAlone,
         y0=AmongGroupsAlone$Phi-NUMSD*AmongGroupsAlone$sd,
         y1=AmongGroupsAlone$Phi+NUMSD*AmongGroupsAlone$sd)
points(xGroupReg,AmongGroupsWithinRegions$Phi,pch=16)
segments(x0=xGroupReg,x1=xGroupReg,
         y0=AmongGroupsWithinRegions$Phi-NUMSD*AmongGroupsWithinRegions$sd,
         y1=AmongGroupsWithinRegions$Phi+NUMSD*AmongGroupsWithinRegions$sd)

points(xReg,AmongRegions$Phi,pch=16)
segments(x0=xReg,x1=xReg,
         y0=AmongRegions$Phi-NUMSD*AmongRegions$sd,
         y1=AmongRegions$Phi+NUMSD*AmongRegions$sd)

abline(v=c(xAlone,xAlone[length(xAlone)]+3)-0.5,lty=2)



####add points to show which pops were used
points(x=rep(1,length(listPops_restricted[["UV & CV"]])),
       y=yPoints[c(1:length(listPops_restricted[["UV & CV"]]))],
       pch=ColorManifest$Point[ ColorManifest$Population %in%listPops_restricted[["UV & CV"]] ],
       bg=ColorManifest$Color[ ColorManifest$Population %in%listPops_restricted[["UV & CV"]] ],
       col=ifelse(ColorManifest$Point[ ColorManifest$Population %in%listPops_restricted[["UV & CV"]] ]<21,
                  ColorManifest$Color[ ColorManifest$Population %in%listPops_restricted[["UV & CV"]] ],
                  "black"),cex=0.8)
x=1
for(Region in AmongGroupsWithinRegions$Region){
  x=x+1
  print(Region)
  print(listPops_restricted[[Region]])
  points(x=rep(x,length(listPops_restricted[[Region]])+length(listPops_restricted[["UV & CV"]])),
         y=yPoints[c(1:(length(listPops_restricted[[Region]])+length(listPops_restricted[["UV & CV"]])))],
         pch=c(rep(NA,length(listPops_restricted[["UV & CV"]])),ColorManifest$Point[ ColorManifest$Population %in%listPops_restricted[[Region]] ]),
         bg=c(rep(NA,length(listPops_restricted[["UV & CV"]])),ColorManifest$Color[ ColorManifest$Population %in%listPops_restricted[[Region]] ]),
         col=c(rep(NA,length(listPops_restricted[["UV & CV"]])),ifelse(ColorManifest$Point[ ColorManifest$Population %in%listPops_restricted[[Region]] ]<21,
                    ColorManifest$Color[ ColorManifest$Population %in%listPops_restricted[[Region]] ],
                    "black")),
         cex=0.8)
  x=x+1
  points(x=rep(x,length(listPops_restricted[[Region]])+length(listPops_restricted[["UV & CV"]])),
         y=yPoints[c(1:(length(listPops_restricted[[Region]])+length(listPops_restricted[["UV & CV"]])))],
         pch=ColorManifest$Point[ ColorManifest$Population %in% c(listPops_restricted[[Region]],listPops_restricted[["UV & CV"]]) ],
         bg=ColorManifest$Color[ ColorManifest$Population %in% c(listPops_restricted[[Region]],listPops_restricted[["UV & CV"]]) ] ,
         col=ifelse(ColorManifest$Point[ ColorManifest$Population %in% c(listPops_restricted[[Region]],listPops_restricted[["UV & CV"]]) ]<21,
                    ColorManifest$Color[ ColorManifest$Population %in% c(listPops_restricted[[Region]],listPops_restricted[["UV & CV"]]) ],
                    "black"),
         cex=0.8)
  x=x+1
  points(x=rep(x,length(listPops_restricted[[Region]])+length(listPops_restricted[["UV & CV"]])),
         y=yPoints[c(1:(length(listPops_restricted[[Region]])+length(listPops_restricted[["UV & CV"]])))],
         pch=ColorManifest$Point[ ColorManifest$Population %in% c(listPops_restricted[[Region]],listPops_restricted[["UV & CV"]]) ],
         bg=ColorManifest$Color[ ColorManifest$Population %in% c(listPops_restricted[[Region]],listPops_restricted[["UV & CV"]]) ] ,
         col=ifelse(ColorManifest$Point[ ColorManifest$Population %in% c(listPops_restricted[[Region]],listPops_restricted[["UV & CV"]]) ]<21,
                    ColorManifest$Color[ ColorManifest$Population %in% c(listPops_restricted[[Region]],listPops_restricted[["UV & CV"]]) ],
                    "black"),
         cex=0.8)
  # break
}

box()
dev.off()

pdf("corSD-phi.pdf")

plot(c(AmongGroupsAlone$Phi,AmongGroupsWithinRegions$Phi,AmongRegions$Phi),
         c(AmongGroupsAlone$sd,AmongGroupsWithinRegions$sd,AmongRegions$sd),
          pch=16,
         col=c(rep(1,nrow(AmongGroupsAlone)),
             rep(2,nrow(AmongGroupsWithinRegions)),
             rep(3,nrow(AmongRegions))),
     main=paste("Correlation btw Phi & its permutation-based SD\nr2 = ",round(cor.test(c(AmongGroupsAlone$Phi,AmongGroupsWithinRegions$Phi,AmongRegions$Phi),
                                       c(AmongGroupsAlone$sd,AmongGroupsWithinRegions$sd,AmongRegions$sd),
                                       pch=c(rep(15,nrow(AmongGroupsAlone))))$estimate,digits = 3)),
     xlab="Phi",ylab="Permutaion-based SD")
legend("topleft",
       pch=16,
       col=c(1,2,3),
       legend = paste(c("Among Groups in Region",
                        "Among Groups within Regions",
                        "Between Region and UV&CV"),
                      " (r2 = ",
                      c(round(cor.test(AmongGroupsAlone$Phi,AmongGroupsAlone$sd)$estimate,digits=3),
                        round(cor.test(AmongGroupsWithinRegions$Phi,AmongGroupsWithinRegions$sd)$estimate,digits=3),
                        round(cor.test(AmongRegions$Phi,AmongRegions$sd)$estimate,digits=3)),
                        ")",sep=""))

lmAGA<-lm(sd~Phi,data=AmongGroupsAlone)$coefficients
lmAGR<-lm(sd~Phi,data=AmongGroupsWithinRegions)$coefficients
lmBR<-lm(sd~Phi,data=AmongRegions)$coefficients
segments(x0=min(AmongGroupsAlone$Phi),
         x1=max(AmongGroupsAlone$Phi),
         y0=lmAGA[1]+min(AmongGroupsAlone$Phi)*lmAGA[2],
         y1=lmAGA[1]+max(AmongGroupsAlone$Phi)*lmAGA[2],
         col=1)
segments(x0=min(AmongGroupsWithinRegions$Phi),
         x1=max(AmongGroupsWithinRegions$Phi),
         y0=lmAGR[1]+min(AmongGroupsWithinRegions$Phi)*lmAGR[2],
         y1=lmAGR[1]+max(AmongGroupsWithinRegions$Phi)*lmAGR[2],
         col=2)
segments(x0=min(AmongRegions$Phi),
         x1=max(AmongRegions$Phi),
         y0=lmBR[1]+min(AmongRegions$Phi)*lmBR[2],
         y1=lmBR[1]+max(AmongRegions$Phi)*lmBR[2],
         col=3)

dev.off()


#### make table test
tableTest<-c()
tablePhi<-c()

for(i in names(amova_test_restricted)){
  
  
  ### recup composant of variance
  compVar<-amova_result_restricted[[i]]$componentsofcovariance
  compVar<-compVar[c(nrow(compVar):1),]
  compVar<-compVar[-1,]
  tmp<-cbind("Regions included"=i,
             "Type"=str_remove(str_remove(row.names(compVar),"Variations  "),"  "),
             "Obs.Var"=amova_test_restricted[[i]]$obs,
             "Obs.Var2"=compVar[,1],
             "PercentageOfVariance"=compVar[,2],
             "alt.hypothesis"=amova_test_restricted[[i]]$alter,
             "pvalue"=amova_test_restricted[[i]]$pvalue)
  tmp<-data.frame(tmp[ tmp[,2]!="Within samples",],stringsAsFactors = F)
  if(sum(tmp$Obs.Var!=tmp$Obs.Var2)>0){
    stop("pb")
  }
  tableTest<-rbind(tableTest,tmp[,-4])
  
  ###recup phi
  Type<-rownames(amova_result_restricted[[i]]$statphi)
  whichType<-which(amova_result_restricted[[i]]$statphi !=1)
  Type=Type[whichType]
  sd<-c()
  for( tyty in Type){
    sd[tyty]<-sd(amova_perm_restricted[[i]][,tyty])
  }
  tmp<-cbind("Regions included"=i,
             "Type"=str_replace(str_replace(str_replace(Type,"Phi-Group-total","Between Group"),"Phi-Region-total","Between Region"),"Phi-Group-Region","Between Group Within Region"),
             "Phi"=amova_result_restricted[[i]]$statphi[whichType,],
             "sd"=sd)
  tmp<-data.frame(tmp,stringsAsFactors = F)
  tablePhi<-rbind(tablePhi,tmp)
  
  
}


outTable<-merge(tablePhi,tableTest,by=names(tableTest)[c(1,2)],all.y=T)
outTable$Phi[ is.na(outTable$Phi)]<-1
outTable$Phi<-as.numeric(outTable$Phi)
outTable$Phi<-round(outTable$Phi,digits=4)
outTable$sd<-as.numeric(outTable$sd)
outTable$sd<-round(outTable$sd,digits=4)
outTable$Obs.Var<-as.numeric(outTable$Obs.Var)
outTable$Obs.Var<-round(outTable$Obs.Var)
outTable$PercentageOfVariance<-as.numeric(outTable$PercentageOfVariance)
outTable$PercentageOfVariance<-round(outTable$PercentageOfVariance,digits=2)

outTable$Type<-ifelse(outTable$Type=="Between Group","Among population groups",
                      ifelse(outTable$Type=="Between samples Within Group", "Among samples within population group",
                             ifelse(outTable$Type=="Between Group Within Region" ,"Among population groups within major group",
                                    ifelse(outTable$Type=="Between Region","Between major groups","WTF"))))
write.table(outTable,"AMOVA_results_Summarized.tsv",col.names = T,row.names=F,sep="\t",quote=F)
  

