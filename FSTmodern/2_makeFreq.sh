#!/bin/bash

module load plink/1.90b6.16 

rootFILE=$(pwd)/../DataSets/
TH=30000

setIND=Lab_with_Compendium_GEHmodern_GenoModern
setSNP=1240K

prefIN=$rootFILE/$setIND.$setSNP/ForAdmixture.TH$TH.GENO0.5.MAF0.01.pruned


mkdir DiffAF

grep 'Luisi\|Reich\|delaFuente\|Homburger' $prefIN.fam > DiffAF/Modern.KEEP
grep 'AR0'  $prefIN.fam | grep -v LAR | grep -v EGAR >> DiffAF/Uspallata.KEEP

cat DiffAF/{Modern,Uspallata}.KEEP > DiffAF/IndSet.KEEP

plink --bfile $prefIN --keep DiffAF/IndSet.KEEP --make-bed --out DiffAF/forDiffAF

awk '{if($1==$2){$1="Uspallata"};$6=$1;print $1,$2,$3,$4,$5,$6}'  DiffAF/forDiffAF.fam >  DiffAF/forDiffAF.TMP
mv DiffAF/forDiffAF.{TMP,fam}

for i in $(awk '{print $1}' DiffAF/forDiffAF.fam | sort | uniq )
do
	break
	mkdir DiffAF/$i
	awk -v i=$i '{if($1==i)print $0}' DiffAF/forDiffAF.fam > DiffAF/$i/$i.KEEP
	plink --bfile DiffAF/forDiffAF --keep DiffAF/$i/$i.KEEP --keep-allele-order --freq --out DiffAF/$i/$i
done
for i in $(awk '{print $1}' DiffAF/forDiffAF.fam | sort | uniq | grep -v Uspallata)
do
	break
	cat DiffAF/$i/$i.KEEP DiffAF/Uspallata/Uspallata.KEEP > DiffAF/$i/$i.UV.KEEP
	plink --bfile DiffAF/forDiffAF --keep DiffAF/$i/$i.UV.KEEP --make-bed --out DiffAF/$i/$i.UV
	awk '{if($1=="Uspallata"){$6=2}else{$6=1};print "Together",$2,$3,$4,$5,$6}' DiffAF/$i/$i.UV.fam > DiffAF/$i/$i.UV.TMP
	mv DiffAF/$i/$i.UV.{TMP,fam}
	plink --bfile DiffAF/$i/$i.UV -maf 0.01 --fst 'case-control' --out DiffAF/$i/$i.UV
done


setIND=Lab_with_Compendium_GEHmodern_GenoModernUnmasked
prefIN=$rootFILE/$setIND.$setSNP/ForAdmixture.TH$TH.GENO0.5.MAF0.01.pruned

grep 'French\|Mbuti' $prefIN.fam > DiffAF/NotAm.KEEP


plink --bfile $prefIN --keep DiffAF/NotAm.KEEP --make-bed --out DiffAF/NotAm
for i in $(awk '{print $1}' DiffAF/NotAm.fam
do
	mkdir DiffAF/$i
	awk -v i=$i '{if($1==i)print $0}' DiffAF/NotAm.fam > DiffAF/$i/$i.KEEP
	plink --bfile DiffAF/NotAm --keep DiffAF/$i/$i.KEEP --keep-allele-order --freq --out DiffAF/$i/$i
done
