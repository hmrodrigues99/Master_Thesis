#!/bin/bash

export PATH=$PATH:~/meme/bin
export PATH=$PATH:~/data/scripts

usage()
{
cat << EOF

usage $0 <options>

A TF-Target Validator based on TFBS (transcription factor binding sites) being present on the promoter region sequences of their putative targets

OPTION
 -h	show this Help message
 -n	path to file containing the Co-expression network (delimiter=" ")
 -p	path to file containing the promoter sequences of all genes
 -t	string with the TF ID (ex: LOC1234) from Quercus suber
 -v	minimum pvalue to establish a valid match in FIMO (default = 1.0E-4
 -o	output directory name (default = TF_Targets_match)

DETAILS

This script receives as input: 1)The co-expression network (3 mandatory columns: LOC1 | LOC2 | edge_score/irp_score); 2)The table with the gene promoter sequences (2 mandatory columns: Gene | Sequence); 3)A list with TF IDs (ex:LOC123 OR AT123).
It creates a output directory with severall output files, one of which is a table with 2 columns: TF | Validated Target
The TF putative targets are the genes that the TF is linked to by an edge in the network, and the validation conducted in this script is a check using FIMO (MEME Suite) of a motif match between the TFBS (transcription factor binding site) and the promoter sequences of its putative targets.

I conducted my work using this script with the following input files: --network Cork_network.txt --promoter_sequences Promoter_sequences_df.csv --TF [relevant TFs from the previous TF_analysis step] --TFBS_matrix TFBS_TF_[relevant TF LOCs]_fasta.fa
I stored my output tables into ...

EOF
}

printf "\n\n --------------------  HELLO AND WELCOME TO THE TF-TARGET VALIDATOR  --------------------- \n\n"

#get options
while getopts "h:n:p:t:v:o:" OPTION
do
   case $OPTION in
     h) usage; exit 1;;
     n) network=$OPTARG;;
     p) promoters=$OPTARG;;
     t) tf=$OPTARG;;
     v) pvalue=$OPTARG;;
     o) outdir=$OPTARG;;
     ?) usage; exit;;
    esac
done

#Check that all required arguments were passed

if [[ -z $network ]] || [[ -z $promoters ]] || [[ -z $tf ]]
then
	printf "\n-----------------\n  --  ERROR: options -n ; -p and -t are required! -- \n------------------\n\n"
	usage
	exit 1
fi

if [[ -z $pvalue ]];
then
	printf " | Using the FIMO default pvalue threshold of 1.0E-4                        |\n"
	pvalue="1.0E-4"
else
	pvalue = $pvalue
	printf " | Using the user defined pvalue treshold = $pvalue                            |\n"
fi

if [[ -z $outdir ]];
then
	printf " | Using the default output directory = TF_Target_Validator_out             |\n"
	outdir="TF_Target_Validator_out"
else
	outdir=$outdir
	printf " | Outputting the files to the Out-directory = $outdir                         |\n"
fi


#   ---------------------------------       MAIN CODE     --------------------------------   #

#Create an output directory
mkdir -p $outdir

#Get the correspondent AT IDs from the LOC IDs
while read f
do
	grep $f arabidopsis_concise.txt > tf_loc.txt
	sed 's/.*\t.*\tAT/AT/' tf_loc.txt > at_ID.txt
	at_tf=$(<at_ID.txt)
	rm tf_loc.txt
	rm at_ID.txt

	#This script receives the TF ID and returns a file with the correspondent TFBS matrix
	printf "\n | The TFBS of ${f} is being searched within the JASPAR database using your input TF ID query   | \n"i
	#The following subscript stopped working properly
	#get_TFBS.sh $at_tf
	curl -X GET "https://jaspar.genereg.net/api/v1/matrix/?search=$at_tf" -H "accept: application/json" -H "X-CSRFToken: NzxtdFE3eih6q6L4ALnKx9Jf00L50eQO17X3J1cyFlS44E8svSYy7CBl90H1bgyZ" > ./matrixID1.txt

	grep -o 'base_id.*version' ./matrixID1.txt > ./matrixID2.txt
	sed 's/base_id":"//' ./matrixID2.txt > ./matrixID3.txt
	sed 's/","version//' ./matrixID3.txt > ./matrixID.txt

	matrixID=$(<matrixID.txt)
	jaspar_out=$(<matrixID1.txt)

	#Downloading the TFBS matrix by the matrix ID
	wget https://jaspar.genereg.net/api/v1/matrix/$matrixID.meme
	mv ./$matrixID.meme ./$at_tf.meme

	rm matrixID1.txt
	rm ./matrixID2.txt
	rm ./matrixID3.txt
	rm ./matrixID.txt

	#Second Queried Database - PlantTFDB
	#If the TFBS does not exist within the JASPAR database, the script tries another search in PlantTFDB
	if [ $jaspar_out=='{"count":0,"next":null,"previous":null,"results":[]}' ];
	then
        	printf "\n | TF ${f} NOT FOUND ON JASPAR DATABASE, SEARCHING IN PLANTTFDB INSTEAD             | \n"
        	rm ./${at_tf}.meme
        	wget http://planttfdb.gao-lab.org/motif/Ath/${at_tf}.meme
	fi
	
	#Worst case possible where the script doesnt find the TFBS in all of the queried databases
	FILE=./${at_tf}.meme
	if test -f "$FILE"; then
		mv ./$at_tf.meme ./$outdir/$at_tf.meme

		#This script receives all the previous inputs plus the TFBS matrix to prepare a fasta file ready to enter into FIMO
		match_TFBS_promoter.py --network $network --promoter_sequences $promoters --TF $f
		mv ./TFBS_TF_${f}_fasta.fa ./$outdir/TFBS_TF_${f}_fasta.fa
		grep '>' ./$outdir/TFBS_TF_${f}_fasta.fa > ./$outdir/TF_${f}_target_list1.txt
        	sed 's/>//g' ./$outdir/TF_${f}_target_list1.txt > ./$outdir/TF_${f}_target_list2.txt
		sed 's/ .*//g' ./$outdir/TF_${f}_target_list2.txt > ./$outdir/TF_${f}_target_list.txt
		rm ./$outdir/TF_${f}_target_list1.txt
		rm ./$outdir/TF_${f}_target_list2.txt	

		#This line runs the FIMO tool from MEME Suite and performs the motif match between the TFBS and the promoter region sequences of its putative targets
		printf "  | Running FIMO on TF ${f} for motif matching          | \n"
		fimo --oc ./$outdir/fimo_${f} --verbosity 1 --thresh $pvalue ./$outdir/$at_tf.meme ./$outdir/TFBS_TF_${f}_fasta.fa

		printf "\n  |  Printing all matching targets of the TF ${f}                  | \n"
		#cat ./$outdir/fimo_${f}/fimo.gff | awk '{print $1,$6}' > ./$outdir/motif_match_${f}.txt
		cat ./$outdir/fimo_${f}/fimo.gff | awk '{print $1}' > ./$outdir/motif_match_${f}.txt
        	awk -F' ' '!a[$1]++' ./$outdir/motif_match_${f}.txt > ./$outdir/best_motif_match_${f}.txt
	else
		printf "${f}	NO_TFBS\n" > ./$outdir/TF_${f}_target_list.txt
		printf "${f}	NO_TFBS\n" > ./$outdir/TF_${f}_matches.txt
		NO_TFBS=TRUE
	fi

	#Starts building the Output File Table based on the FIMO results. YES - match | NO - no match | NO_TFBS - no TFBS found
	while read target
	do
		if grep -Fxq "$target" ./$outdir/best_motif_match_${f}.txt
		then
			printf "${target}	YES\n" >> ./$outdir/TF_${f}_matches.txt
		else
			printf "${target}	NO\n" >> ./$outdir/TF_${f}_matches.txt
		fi
	done < ./$outdir/TF_${f}_target_list.txt	
done < $tf

