#!/bin/bash

module load plink/1.90b6.16 

rootFILE=$(pwd)/../DataSets/
TH=30000

setIND=Lab_with_Compendium_GEHmodern_GenoModern
setSNP=1240K

prefIN=$rootFILE/$setIND.$setSNP/ForAdmixture.TH$TH.GENO0.5.MAF0.01.pruned


mkdir DiffAFshuffled

grep 'Luisi\|Reich\|delaFuente\|Homburger' $prefIN.fam > DiffAF/Modern.KEEP
grep 'AR0'  $prefIN.fam | grep -v LAR | grep -v EGAR >> DiffAF/Uspallata.KEEP

cat DiffAF/{Modern,Uspallata}.KEEP > DiffAF/IndSet.KEEP

plink --bfile $prefIN --keep DiffAF/IndSet.KEEP --make-bed --out DiffAF/forDiffAF

awk '{if($1==$2){$1="Uspallata"};$6=$1;print $1,$2,$3,$4,$5,$6}'  DiffAF/forDiffAF.fam >  DiffAF/forDiffAF.TMP
mv DiffAF/forDiffAF.{TMP,fam}

for comp in "Aymara-Reich___Cabecar-Reich" "Aymara-Reich___Chono-Reich" "Aymara-Reich___Diaguita-Reich" "Aymara-Reich___Guarani-Reich" "Aymara-Reich___Huilliche-delaFuente"  "Aymara-Reich___Kaweskar-delaFuente"  "Cabecar-Reich___Chono-Reich"  "Cabecar-Reich___Diaguita-Reich" "Cabecar-Reich___Guarani-Reich" "Cabecar-Reich___Huilliche-delaFuente" "Cabecar-Reich___Kaweskar-delaFuente"  "Chono-Reich___Diaguita-Reich" "Chono-Reich___Guarani-Reich"  "Chono-Reich___Huilliche-delaFuente" "Chono-Reich___Kaweskar-delaFuente" "Diaguita-Reich___Guarani-Reich" "Diaguita-Reich___Huilliche-delaFuente"  "Diaguita-Reich___Kaweskar-delaFuente"  "Guarani-Reich___Huilliche-delaFuente" "Guarani-Reich___Kaweskar-delaFuente"  "Huilliche-delaFuente___Kaweskar-delaFuente"
do
	
	mkdir DiffAFshuffled/$comp
	pop1=$(echo $comp | awk -F "___" '{print $1}')
	pop2=$(echo $comp | awk -F "___" '{print $2}')
	awk -v i=$pop1 -v j=$pop2  '{if($1==i || $1==j)print $0}' DiffAF/forDiffAF.fam > DiffAFshuffled/$comp/$comp.KEEP
	plink --bfile DiffAF/forDiffAF --keep DiffAFshuffled/$comp/$comp.KEEP --keep-allele-order --make-bed --out DiffAFshuffled/$comp/$comp
	head DiffAFshuffled/$comp/$comp.fam
	tail DiffAFshuffled/$comp/$comp.fam
	echo "######"
	Rscript shuffleFam.R DiffAFshuffled/$comp/$comp.fam
	head DiffAFshuffled/$comp/$comp.fam
	tail DiffAFshuffled/$comp/$comp.fam
	grep $pop1 DiffAFshuffled/$comp/$comp.fam > DiffAFshuffled/$comp/$pop1.KEEP 
	grep $pop2 DiffAFshuffled/$comp/$comp.fam > DiffAFshuffled/$comp/$pop2.KEEP
	plink --bfile DiffAFshuffled/$comp/$comp --keep  DiffAFshuffled/$comp/$pop1.KEEP  --keep-allele-order --freq --out DiffAFshuffled/$comp/$pop1
	plink --bfile DiffAFshuffled/$comp/$comp --keep  DiffAFshuffled/$comp/$pop2.KEEP  --keep-allele-order --freq --out DiffAFshuffled/$comp/$pop2
done

exit
