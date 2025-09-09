#!/bin/Rscript

setwd("~/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F/Maps/")
startInd<-list(
  "SG"=read.table("../../DataSets/Lab_with_Compendium_GEHmodern.SG//finalSet.TH50000.ind.txt",stringsAsFactors = F,header=F)[,c(1,3)],
  "SG.TVs"=read.table("../../DataSets/Lab_with_Compendium_GEHmodern.SG.TVs//finalSet.TH50000.ind.txt",stringsAsFactors = F,header=F)[,c(1,3)],
  "1240K"=read.table("../..//DataSets/Lab_with_Compendium_GEHmodern.1240K//finalSet.TH50000.ind.txt",stringsAsFactors = F,header=F)[,c(1,3)],
  "1240K.TVs"=read.table("../../DataSets/Lab_with_Compendium_GEHmodern.1240K.TVs//finalSet.TH50000.ind.txt",stringsAsFactors = F,header=F)[,c(1,3)],
  "f4modern"=read.table("../../DataSets/Lab_with_Compendium_GEHmodern_GenoModern.1240K/finalSet.TH30000.ind.txt",stringsAsFactors = F,header=F)[,c(1,3)],
  "AdmAncient"=read.table("../../DataSets/Lab_with_Compendium.1240K/ForAdmixture.TH30000.MAF0.01.pruned.fam",stringsAsFactors = F,header=F)[,c(2,1)],
  "AdmModern"=read.table("../../DataSets/Lab_with_Compendium_GEHmodern_GenoModern.1240K/ForAdmixture.TH30000.GENO0.5.MAF0.01.pruned.fam",stringsAsFactors = F,header=F)[,c(2,1)])

for(i in names(startInd)){
  names(startInd[[i]])<-c("Individual","Population")
}



ALL<-read.table("../../DataSets_allinds/Lab_with_Compendium_GEHmodern_GenoModern.1240K/finalSet.ind.txt",stringsAsFactors = F,header=F)[,c(1,3)]

names(ALL)<-c("Individual","Population")
ALL<-ALL[ ! ALL$Individual %in% c("I5319","I14922","I14879","I15601","I15676") ,]
annot<-read.table("../../../../ReferenceDataSet/ColorCompendium_ModernAndAncient_NoReichSG.tsv",stringsAsFactors = F,header=T,sep="\t",comment.char = "@",quote="\"")
annot<-unique(annot[,c("Study","Region","Population","Color","Point")])
names(annot)[c(2,3)]<-c("MajorGroup","Group")
annotProj<-read.table("../../../Uspallata_Annotation.tsv",stringsAsFactors = F,header=T,sep="\t",comment.char = "@",quote="\"")

annotProj$MainRegion2<-ifelse(annotProj$MainRegion=="Migrant","LH-MF",
                             ifelse(annotProj$MainRegion=="Migrant_Outlier","LH-MFout",
                                    ifelse(annotProj$MainRegion=="Local_PLC","LH-LFplc",
                                           ifelse(annotProj$MainRegion=="Local_Early","LH-HG",
                                                  ifelse(annotProj$MainRegion=="Local_Late","LH-LF","WTF")))))
annotProj$MainRegion1<-"Uspallata"
annotProj$MainRegion3<-ifelse(annotProj$MainRegion2=="LH-HG","LH-HG","LH-F")

annotProj$Group=paste(annotProj$Site," (",annotProj$MainRegion2,"|",annotProj$MainRegion3,")",sep="")
annotProj$Study<-"PresentStudy"
annotProj$MajorGroup<-"Uspallata"

indREF<-merge(annot,ALL,by.x="Group",by.y="Population")
indREF<-indREF[ order(indREF$Study,indREF$MajorGroup,indREF$Group),]
indPROJ<-merge(ALL,annotProj[,c("Individual","MajorGroup","Group","Study","Color","Point")],by="Individual")

indOUT<-rbind(indPROJ[,names(indREF)],indREF)

###Sex Determination
sex<-read.table("../../../SexDeterminded.tsv",stringsAsFactors = F,header=T)
indOUT$SexDetermination<-indOUT$Individual %in% sex$Sample

###Kinship
kin<-read.table("../../../Kinship/Summarize_Kinship.tsv",stringsAsFactors = F,header=T,sep="\t")
indOUT$Kinship<-indOUT$Individual %in% c(kin$Individual1,kin$Individual2)

###F3
f3<-read.table("../..//Final_Plots_F/F3_IND/Lab_with_Compendium_GEHmodern_MDS.NoNorthernNorthAmerica.NoCarribes.Ancient.tsv",stringsAsFactors = F,header=T,sep="\t")
indOUT$FstatsAncient_SG<-indOUT$Individual %in% f3$Individual[ ! is.na(f3$Dim1_SG)]
indOUT$FstatsAncient_SG.TVs<-indOUT$Individual %in% f3$Individual[ ! is.na(f3$Dim1_SG.TVs)]
indOUT$FstatsAncient_1240K<-indOUT$Individual %in% f3$Individual[ ! is.na(f3$Dim1_1240K)]
indOUT$FstatsAncient_1240K.TVs<-indOUT$Individual %in% f3$Individual[ ! is.na(f3$Dim1_1240K.TVs)]


##admixture ancient
indOUT$AdmixtureAncient<-indOUT$Individual %in% startInd[["AdmAncient"]]$Individual


### f4 modern
f4_modern<-read.table("../../Final_Plots_F/F4_Modern_toSouthRefUspallata/F4.Modern_toSouthRefandUspallata.finalSet.Lab_with_Compendium_GEHmodern_GenoModern.1240K.TH30000.tsv",stringsAsFactors = F,header=T,sep="\t")
indOUT$F4_Modern<-indOUT$Individual %in% startInd[["f4modern"]]$Individual & ( indOUT$MajorGroup =="Uspallata" | indOUT$Group %in% unique(c(f4_modern$Pop4_Z,f4_modern$Pop2_X)))

##admixture Modern
indOUT$AdmixtureWithModern<-indOUT$Individual %in% startInd[["AdmModern"]]$Individual

###AMOVA
AMOVA=read.table("../../FST/indINCLUDED_forAMOVA.txt",stringsAsFactors = F,header=T,sep="\t")[,c(1,2)]
indOUT$AMOVA<-indOUT$Individual %in% AMOVA$id


###condHet 
listPOPs_SG<-read.table("../../Final_Plots_F/NeHapROH_vsCondHE/Only.SG.TH50000.ListPops.tsv",stringsAsFactors = F,header=F)$V1
listPOPs_1240K<-read.table("../../Final_Plots_F/NeHapROH_vsCondHE/Only.1240K.TH50000.ListPops.tsv",stringsAsFactors = F,header=F)$V1
listPOPs_SG<-c(listPOPs_SG,unique(indOUT$Group[ indOUT$Study=="PresentStudy"]))

listPOPs_1240K<-c(listPOPs_1240K,unique(indOUT$Group[ indOUT$Study=="PresentStudy"]))


indOUT$ConditionalHeterozygosity_SG<-indOUT$Individual %in% startInd[["SG"]]$Individual & indOUT$Group %in%listPOPs_SG
indOUT$ConditionalHeterozygosity_1240K<-indOUT$Individual %in% startInd[["1240K"]]$Individual & indOUT$Group %in%listPOPs_1240K

###Ne
listPOPs<-read.table("../../Final_Plots_F/NeHapROH_vsCondHE/OnlyNefromHapROH.ListPops.tsv",stringsAsFactors = F,header=F,sep="\t")$V1
listPOPs<-c(listPOPs,unique(indOUT$Group[ indOUT$Study=="PresentStudy"]))

listIND<-rbind(read.table("../..//hapROH/NeEstimate_TH400000/Ne_estimates.IndIncludedLocal_Late.tsv",stringsAsFactors = F,header=T,sep="\t"),
      read.table("../../hapROH/NeEstimate_TH400000/Ne_estimates.IndIncludedMigrant.tsv",stringsAsFactors = F,header=T,sep="\t"))$iid
listIND<-c(listIND,read.table("../../../../ReferenceDataSet/CompendiumPostKin_Analyses_NoReichSG/hapROH/ListINDincluded.forNeEstimates.tsv",stringsAsFactors = F,header=F)$V2)

indOUT$Ne.hapROH<-indOUT$Individual %in% listIND & indOUT$Group %in%listPOPs


### hapNE
keep<-read.table("../../HapNe/Lab_with_Compendium.1240K/TH50000/Migrant/Migrant.KEEP",stringsAsFactors = F,header=F)
indOUT$HapNe<-indOUT$Individual %in% keep$V2


indOUT$Included<-apply(indOUT[,-c(1:6)],1,sum,USE.NAMES = F)

indOUT<-indOUT[ indOUT$Included>0 | indOUT$MajorGroup=="Uspallata", names(indOUT)!="Included"]

indOUT<-indOUT[,c(2,3,1,c(4:length(indOUT)))]

write.table(indOUT,"AnnotIndividualsIncluded.tsv",col.names=T,row.names=F,sep="\t",quote=F)


###check
indOUT$USED<-apply(indOUT[,-c(1:6)],1,sum)
