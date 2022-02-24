#!/bin/bash
#SBATCH --job-name=%j_apt
#SBATCH --ntasks=40
#SBATCH --nodelist=compute1

##############################################################################################################
################## Defining the path variables ########################################################
##############################################################################################################

#enter the name of the project
proj=Ria_ICC

#Enter the path to the the list of DEGs
list=~/ICC

#Enter the path to the list of Transcription Factors (TFs)
pwms=~/Homo_sapiens_2021_04_20_2_40_am/pwms_all_motifs

#enter the path to the file containing TFs from cis-bp
tf=~/Homo_sapiens_2021_04_20_2_40_am/TF_Information.txt

#Define the path to the folder containing the PWMs 

mkdir ~/$proj/Sequences/background/
mkdir ~/$proj/Sequences/meme_motifs/
mkdir ~/$proj/Results/
mkdir ~/$proj/extra/
mkdir ~/$proj/TFs/

sequence=~/$proj/Sequences
pseudog=~/$proj/extra
result=~/$proj/Results
mscan=~/$proj/TFs
echo "Paths stored"

################################################################################################################
#################### Editing the genelist for any windows-based error/duplications #############################
################################################################################################################
:<<'END'
tr -d '\15\32' < $list > genelist_w1

tr '[:lower:]' '[:upper:]' < genelist_w1 > genelist_w2

sort -u genelist_w2 > $list

rm genelist_w1

rm genelist_w2 

#cp $list/genelist $mscan/genelist
END
###############################################################################################################
################### Retrieving the start and end coordinates of genes #########################################
###############################################################################################################

#The subsequent step requires RSAT installation

#retrieve-seq  -org Arabidopsis_thaliana.TAIR10.42 -feattype gene -type upstream -format fasta -label id,name -from -1000 -to -1 -noorf -i $list -o $sequence/ref-seq

grep -v "WARNING" $sequence/ref-seq > $sequence/fasta
cat $sequence/fasta > $sequence/ref-seq
rm $sequence/fasta 

echo "Upstream regions of the target genes extracted"

cd $pseudog
grep -i ">" $sequence/ref-seq | sed "s/[|>:;]/ /g" | awk '{print $2,$21,$22,$23,$2}' | awk '{if($2<$3) print}' > $pseudog/Gene_Coord
sed -i 's/D/+/g' $pseudog/Gene_Coord
sed -i 's/R/-/g' $pseudog/Gene_Coord
echo "Gene Coordinates of the upstream regions extracted"

echo "Files with upstream coordinates transferred to Pseudogenomes"

echo "Extracting the list of TFs from the DEGs"

rm $mscan/tflist
while read -r line
do
	grep $line $tf | awk '{if($9 != "N") print $4,$6,$7,$9}' >> $mscan/tflist
done < $list

w=($(wc $mscan/tflist))

if [[ $w -gt 0 ]]
then

	echo "List of transcription factors extracted"

###########################################################################################################################
####################### Compilation of all motifs in a file in TRANSFAC format ############################################
###########################################################################################################################


rm $sequence/collected

#cd $mscan

while read -r line
do

	a=($line)
	motif="${a[0]}.txt"
	w=($(wc $pwms/$motif))
	if [[ $w -gt 1 ]]
	then

		echo "AC ${a[1]}" >> $sequence/collected
		echo "XX" >> $sequence/collected
		echo "ID ${a[2]}" >> $sequence/collected
		echo "XX" >> $sequence/collected
		cat $pwms/$motif >> $sequence/collected
		echo "XX" >> $sequence/collected
	fi
	
done < $mscan/tflist

echo "//" >> $sequence/collected

sed -i "s/Pos/PO/g" $sequence/collected

echo "Concatenation of motifs done"

cd $sequence
ls | grep "sequences" | grep -v "list" > $sequence/list_of_sequences
##################################################################################################################################
################### Scanning of sequences by single compiled motif file using FIMO ###############################################
##################################################################################################################################

export PATH=/home/cluster/kshattry/meme/bin:/home/cluster/kshattry/meme/libexec/meme-5.3.2:$PATH
fasta-get-markov -m 1 -dna $sequence/ref-seq $sequence/background/ref-seq
transfac2meme -bg $sequence/background/ref-seq -use_acc $sequence/collected > $sequence/meme_motifs/ref-seq
fimo -o $result/ref-seq -bfile $sequence/background/ref-seq $sequence/meme_motifs/ref-seq $sequence/ref-seq

#rm $mscan/collected

###################################################################################################################################
###################################### Removing duplicated rows from result #######################################################
###################################################################################################################################

fi


