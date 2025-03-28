

setwd("~/Documents/PostDocPasteur/aDNA/2024-09-01_Uspallata_noHighCov/Analyses/Final_Plots_F/")



tableIND<-read.table("Maps/AnnotIndividualsIncluded.tsv",stringsAsFactors = F,header=T,sep="\t")
###Table A: f3
a1240<-read.table("../F3_IND/Lab_with_Compendium_GEHmodern.1240K/TH50000/Lab_with_Compendium_GEHmodern.1240K.OUT",stringsAsFactors = F,header=T)
a1240<-a1240[a1240$Source1 %in% tableIND$Individual[ tableIND$FstatsAncient_1240K] & a1240$Source2 %in% tableIND$Individual[ tableIND$FstatsAncient_1240K],]
names(a1240)[-c(1:3)]<-paste(names(a1240)[-c(1:3)],"1240K",sep=".")

a1240.tvs<-read.table("../F3_IND/Lab_with_Compendium_GEHmodern.1240K.TVs//TH50000/Lab_with_Compendium_GEHmodern.1240K.TVs.OUT",stringsAsFactors = F,header=T)
a1240.tvs<-a1240.tvs[a1240.tvs$Source1 %in% tableIND$Individual[ tableIND$FstatsAncient_1240K.TVs] & a1240.tvs$Source2 %in% tableIND$Individual[ tableIND$FstatsAncient_1240K.TVs],]
names(a1240.tvs)[-c(1:3)]<-paste(names(a1240.tvs)[-c(1:3)],"1240K.TVs",sep=".")

aSG<-read.table("../F3_IND/Lab_with_Compendium_GEHmodern.SG/TH50000/Lab_with_Compendium_GEHmodern.SG.OUT",stringsAsFactors = F,header=T)
aSG<-aSG[aSG$Source1 %in% tableIND$Individual[ tableIND$FstatsAncient_SG] & aSG$Source2 %in% tableIND$Individual[ tableIND$FstatsAncient_SG],]
names(aSG)[-c(1:3)]<-paste(names(a1240)[-c(1:3)],"SG",sep=".")

aSG.tvs<-read.table("../F3_IND/Lab_with_Compendium_GEHmodern.SG.TVs//TH50000/Lab_with_Compendium_GEHmodern.SG.TVs.OUT",stringsAsFactors = F,header=T)
aSG.tvs<-aSG.tvs[aSG.tvs$Source1 %in% tableIND$Individual[ tableIND$FstatsAncient_SG.TVs] & aSG.tvs$Source2 %in% tableIND$Individual[ tableIND$FstatsAncient_SG.TVs],]
names(aSG.tvs)[-c(1:3)]<-paste(names(aSG.tvs)[-c(1:3)],"SG.TVs",sep=".")

tableA<-merge(a1240,a1240.tvs,by=names(a1240)[c(1:3)],all=T)
tableA<-merge(tableA,aSG,by=names(a1240)[c(1:3)],all=T)
tableA<-merge(tableA,aSG.tvs,by=names(a1240)[c(1:3)],all=T)

tableA<-merge(tableA,tableIND[,c("Individual","MajorGroup","Group","Color","Point")],by.x="Source2",by.y="Individual")
names(tableA)[length(tableA)-c(0:3)]<-paste(names(tableA)[length(tableA)-c(0:3)],"2",sep="")
tableA<-merge(tableA,tableIND[,c("Individual","MajorGroup","Group","Color","Point")],by.x="Source1",by.y="Individual")
names(tableA)[length(tableA)-c(0:3)]<-paste(names(tableA)[length(tableA)-c(0:3)],"1",sep="")

write.table(tableA,"TableF4/tableA_f3inds.tsv",sep="\t",col.names = T,row.names=F,quote=F)

###TableB: MDS from f3
b1240<-read.table("../F3_IND/Lab_with_Compendium_GEHmodern.1240K/TH50000/Lab_with_Compendium_GEHmodern.1240K_MDS.NoNorthernNorthAmerica.tsv",stringsAsFactors = F,header=T,sep="\t")
b1240<-b1240[b1240$Individual %in% tableIND$Individual[ tableIND$FstatsAncient_1240K],]
names(b1240)[c(2:13)]<-paste("Dim",c(1:12),"_1240K",sep="")

b1240.tvs<-read.table("../F3_IND/Lab_with_Compendium_GEHmodern.1240K.TVs//TH50000/Lab_with_Compendium_GEHmodern.1240K.TVs_MDS.NoNorthernNorthAmerica.tsv",stringsAsFactors = F,header=T,sep="\t")
b1240.tvs<-b1240.tvs[b1240.tvs$Individual %in% tableIND$Individual[ tableIND$FstatsAncient_1240K.TVs],]
names(b1240.tvs)[c(2:13)]<-paste("Dim",c(1:12),"_1240K.TVs",sep="")

bSG<-read.table("../F3_IND/Lab_with_Compendium_GEHmodern.SG/TH50000/Lab_with_Compendium_GEHmodern.SG_MDS.NoNorthernNorthAmerica.tsv",stringsAsFactors = F,header=T,sep="\t")
bSG<-bSG[bSG$Individual %in% tableIND$Individual[ tableIND$FstatsAncient_SG],]
names(bSG)[c(2:13)]<-paste("Dim",c(1:12),"_SG",sep="")

bSG.tvs<-read.table("../F3_IND/Lab_with_Compendium_GEHmodern.SG.tvs/TH50000/Lab_with_Compendium_GEHmodern.SG.TVs_MDS.NoNorthernNorthAmerica.tsv",stringsAsFactors = F,header=T,sep="\t")
bSG.tvs<-bSG[bSG$Individual %in% tableIND$Individual[ tableIND$FstatsAncient_SG.TVs],]
names(bSG.tvs)[c(2:13)]<-paste("Dim",c(1:12),"_SG.TVs",sep="")

tableB<-merge(b1240[,c(1:13)],b1240.tvs[,c(1:13)],by="Individual",all=T)
tableB<-merge(tableB,bSG[,c(1:13)],by="Individual",all=T)
tableB<-merge(tableB,bSG.tvs[,c(1:13)],by="Individual",all=T)

tableB<-merge(tableB,tableIND[,c("Individual","MajorGroup","Group","Color","Point")],by="Individual")
write.table(tableB,"TableF4/tableB_MDS.f3inds.tsv",sep="\t",col.names = T,row.names=F,quote=F)

###table C: vs modern
tableC<-read.table("F4_Modern_toSouthRefUspallata/F4.Modern_toSouthRefandUspallata.finalSet.Lab_with_Compendium_GEHmodern_GenoModern.1240K.TH30000.tsv",stringsAsFactors = F,header=T,sep="\t")
names(tableC)[c(5:10)]<-paste(names(tableC)[c(5:10)],"_1240K",sep="")
names(tableC)<-str_replace(names(tableC),"Region","MajorGroup")
write.table(tableC,"TableF4/tableC_f4.modern.tsv",sep="\t",col.names = T,row.names=F,quote=F)

###table D: vs ref vs pairsUV
tableD<-read.table("F4_ref_toPairwiseUspallata/F4.Ref_toPairwiseUspallata.finalSet.Lab_with_Compendium_GEHmodern..TH50000.tsv",stringsAsFactors = F,header=T,sep="\t")
names(tableD)<-str_replace(names(tableD),"Region","MajorGroup")
write.table(tableD,"TableF4/tableD_f4.refs_vsPairsUV.tsv",sep="\t",col.names = T,row.names=F,quote=F)

###table E: UV vs pairREFs
tableE<-read.table("F4_UspallataOneGroup_toRefs//F4.Uspallata_toPairwiseRefs.finalSet.Lab_with_Compendium_GEHmodern.TH50000.tsv",stringsAsFactors = F,header=T,sep="\t")
names(tableE)<-str_replace(names(tableE),"Region","MajorGroup")
write.table(tableE,"TableF4/tableE_f4.UV_vsPairsRefs.tsv",sep="\t",col.names = T,row.names=F,quote=F)

###table F:  UV-SA vs North
tableF<-read.table("F4_NorthRef_toUspallataSouthRef/F4.NorthRef_toUspallataSouthRef.finalSet.Lab_with_Compendium_GEHmodern.SG.TVs.TH50000.Uspallata.summarized.tsv",stringsAsFactors = F,header=T,sep="\t")
names(tableF)<-str_replace(names(tableF),"Region","MajorGroup")
write.table(tableF,"TableF4/tableF_f4.North_vsUV-South.tsv",sep="\t",col.names = T,row.names=F,quote=F)

###table G:  Refs vs LH-MF,LH-F
tableG<-read.table("F4_EarlyLate_vsRefs//F4_EarlyLate_vsSNARefs.tsv",stringsAsFactors = F,header=T,sep="\t")
names(tableG)<-str_replace(names(tableG),"Region","MajorGroup")
write.table(tableG,"TableF4/tableG_f4.EarlyLate_vsRefs.tsv",sep="\t",col.names = T,row.names=F,quote=F)

