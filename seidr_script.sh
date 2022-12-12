#!/bin/bash

export PATH=$PATH:~/bin/seidr/build
export PATH=$PATH:~/data/scripts

usage()
{
cat << EOF

usage: $0 <options>

A Crowd Network generator for visualization of Co-Expression Networks.

OPTIONS:
 -h	show this Help message
 -f     path to file containing the Raw_Counts.txt table (delimiter=" ")
 -n	number of computer threads for the seidr run
 -a     crowd network aggregating method (default = irp)
 -o     output directory of the generated files and name of the output network.txt file (default = seidr_output)
 -p     p-value treshold (default = 1.28, aprox. 10%), for 5% use 1.64 and 1% use 2.32
 -d     depth of the network generation where higher runtime extracts more information (default = SLOW)
 -t     path to file containing a list of genes (file with 1 column and one gene per row), for the generation of a Targeted Network
 -c     cutoff of the crowd network
 -s     outputs a network_stats.txt and a network_nodestats.txt with information about the network and some node statistics, respectively
 -e     outputs an additional separate file only with the edge correlation scores

DETAILS:

This script receives a Raw_Counts.txt table delimited by spaces and computes a Crowd Network via the Seidr toolkit.

As an input the Raw_Counts.txt table with row samples and gene columns should be provided.

There are 4 available aggregating methods of the crowd network (borda; top1; top2; irp), being the irp method implemented by default.

The computation of the network can be a tedious process, therefore 4 depth settings (FAST; MEDIUM; SLOW; VERY_SLOW) can be specified. The default option is set to SLOW. Increasing the depth will extract more information for the generation of the network but will take more time.
FAST setting - only constructs the network using correlation scores (calculated by 3 correlation algorithms).
MEDIUM setting - addition of some mutual information algorithms.
SLOW setting - incorporates machine learning algorithms with higher runtimes.
VERY_SLOW setting - adds one more algorithm at a great time cost (may not be worth the minimal information gain and leads to crashes with a large Raw_Counts.txt).

With the -t option it is possible to specify a gene list to generate a Targeted Network. This type of guided aproach originates a network that contains the gene/nodes present in the target list, and other genes (not specified in the list) if they happen to be co-expressed with the genes on the list.

###
My work uses a Raw_Counts.txt table generated by the Raw_Counts.R script which includes the gene counts for 74 samples from 6 projects sample data.
The final GCN used was obtained using the following arguments: -p (path to)/Raw_Counts.txt -n 20 -a irp -o DEG_Cork_SLOWP -p 2.32 -d SLOW -t DEG_Cork.txt -c 0.4 -s TRUE -e TRUE
###

EOF
}

printf "\n\n ----------------------------- HELLO AND WELCOME TO THE SEIDR SCRIPT ------------------------------- \n\n"

#Get options
#while getopts "h:f:n:a:o:p:d:t:c:s:m:" OPTION
while getopts "h:f:n:a:o:p:d:t:c:s:e:" OPTION
do
  case $OPTION in
    h) usage; exit 1;;
    f) rawcountsfile=$OPTARG;;
    n) thread=$OPTARG;;
    a) aggregate=$OPTARG;;
    o) outdir=$OPTARG;;
    p) pvalue=$OPTARG;;
    d) depth=$OPTARG;;
    t) target=$OPTARG;;
    c) cutoff=$OPTARG;;
    s) stats=$OPTARG;;
    e) edgecorscores=$OPTARG;;
   #m) mcounts=$OPTARG;;
    ?) usage; exit;;
  esac
done

#Check that all required arguments were passed

if [[ -z $rawcountsfile ]] || [[ -z $thread ]]
then
  printf "\n---------------------\n -- ERROR: options -f and -n are required!-- \n--------------------\n\n"
  usage
  exit 1
fi

if [[ -z $aggregate ]];
then
	printf " | Using the Seidr default irp method for aggregating the Crowd Network                              |\n"
	aggregate="irp"
fi

if [[ -z $outdir ]];
then
	printf " | Using the default outdir = seidr_output                                                           |\n"
	outdir="seidr_output"
else
	outdir=$outdir
	printf " | Outputting the files to the Out-Directory = $outdir                                               |\n"
fi

if [[ -z $pvalue ]];
then
	printf " | Applying the default p-value for prunning noise of 1.28 ( p-value = 0.10 )                       |\n"
	pvalue=$pvalue
else
	pvalue=$pvalue
	printf " | Applying the non-default p-value for prunning noise of $pvalue                                   |\n"
fi

if [[ -z $depth ]];
then
	printf " | Applying the default computational speed ( depth = SLOW )                                        |\n"
	depth="SLOW"
fi

if [[ -z $target ]];
then
	printf " | Generating a non-targeted network                                                                |\n"
	target="FALSE"
fi

if [[ -z $cutoff ]];
then
	printf " | Applying cutoof to the network: FALSE (can be done after the network is generated/visualized)|\n"
       cutoff="None"
else
	cutoff=$cutoff
	printf " | Applying cutoff to the network: TRUE (cutoff = $cutoff )                                        |\n"
fi


if [[ -z $stats ]];
then
	stats="TRUE"
elif [[ $stats == 'TRUE' ]];
then
	stats="TRUE"
elif [[ $stats == 'FALSE' ]];
then
	stats="FALSE"
else
	printf " |           ! Your input for the optional argument -stats/-s was invalid. Valid inputs are:     |\n"
	printf " |             TRUE or FALSE. On argument omission or misspelling, the setting TRUE is selected  |\n"
	printf " |             by default for -stats!                                                             |\n"
	stats="TRUE"
fi
printf " | Calculate network and node statistics: $stats                                                         |\n"


if [[ -z $edgecorscores ]];
then
	edgecorscores="TRUE"
elif [[ $edgecorscores == 'TRUE' ]];
then
	edgecorscores="TRUE"
elif [[ $edgecorscores == "FALSE" ]];
then
	edgecorscores="FALSE"
else
	printf " |          ! Your input for the optional argument -edgecorscores/-e was invalid. Valid inputs  |\n"
	printf " |            TRUE or FALSE. On argument omission or misspelling, the setting TRUE is selected  |\n"
	printf " |            by default for -edgecorscores!                                                     |\n"
	edgecorscores="TRUE"
fi
printf " | Generate a separate file with only the correlation scores: $edgecorscores                            |\n"



# ------------------------------------------------------        MAIN CODE        ----------------------------------------------------- #

# Create an output directory
mkdir -p $outdir

#This script receives a Raw_Counts.txt and outputs an expression.tsv and a genes.txt file to the specified outdir.
#These two files are the base material for all calculus for the Seidr crowd network generation toolkit.
raw_counts_normalization_seidr.R --rawcountsfile $rawcountsfile --output $outdir

printf "\n\n | --------------------- Started computation with $thread threads and a $depth depth ---------------------- |\n\n"

#Checks if the network to be generated is in targeted mode or not
if [ $target == 'FALSE' ]
then

#Running the algorithms

# FAST
printf "Calculating the Pearson Correlation scores 'fast'.\n"
correlation -m pearson -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/pearson_scores.tsv
seidr import -A -r -u -n PEARSON -o $PWD/$outdir/pearson_scores.sf -F lm -i $PWD/$outdir/pearson_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the Spearman Correlation scores 'fast'.\n"
correlation -m spearman -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt
seidr import -A -r -u -n SPEARMAN -o $PWD/$outdir/spearman_scores.sf -F lm -i $PWD/$outdir/spearman_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the PCorrelation scores 'fast'.\n"
pcor -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt
seidr import -A -r -u -n PCOR -o $PWD/$outdir/pcor_scores.sf -F lm -i $PWD/$outdir/pcor_scores.tsv -g $PWD/$outdir/genes.txt

#MEDIUM
if [ $depth == 'MEDIUM' ] || [ $depth == 'SLOW' ] || [ $depth == 'VERY_SLOW' ]
then

printf "Calculating the RAW scores 'medium'.\n"
mi -m RAW -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/mi_scores.tsv
seidr import -r -u -n MI -o $PWD/$outdir/mi_scores.sf -F lm -i $PWD/$outdir/mi_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the CLR scores 'medium'.\n"
mi -m CLR -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -M $PWD/$outdir/mi_scores.tsv -o $PWD/$outdir/clr_scores.tsv
seidr import -r -u -z -n CLR -o $PWD/$outdir/clr_scores.sf -F lm -i $PWD/$outdir/clr_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the ARACNE scores 'medium'.\n"
mi -m ARACNE -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -M $PWD/$outdir/mi_scores.tsv -o $PWD/$outdir/aracne_scores.tsv
seidr import -r -u -z -n ARACNE -o $PWD/$outdir/aracne_scores.sf -F lm -i $PWD/$outdir/aracne_scores.tsv -g $PWD/$outdir/genes.txt
fi

#SLOW
if [ $depth == 'SLOW' ] || [ $depth == 'VERY_SLOW' ]
then

printf "Calculating the Narromi scores 'slow'\n"
narromi -O $thread -m interior-point -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/narromi_scores.tsv
seidr import -r -z -n NARROMI -o $PWD/$outdir/narromi_scores.sf -F m -i $PWD/$outdir/narromi_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the Plsnet scores 'slow'\n"
plsnet -O $thread -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/plsnet_scores.tsv
seidr import -r -z -n PLSNET -o $PWD/$outdir/plsnet_scores.sf -F m -i $PWD/$outdir/plsnet_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the LLR-Ensemble scores 'slow'\n"
llr-ensemble -O $thread -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/llr_scores.tsv
seidr import -r -z -n LLR -o $PWD/$outdir/llr_scores.sf -F m -i $PWD/$outdir/llr_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the SVM-Ensemble scores 'slow'\n"
svm-ensemble -O $thread -k POLY -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/svm_scores.tsv
seidr import -r -z -n SVM -o $PWD/$outdir/svm_scores.sf -F m -i $PWD/$outdir/svm_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the Genie3 scores 'slow'.\n"
genie3 -O $thread -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/genie3_scores.tsv
seidr import -r -z -n GENIE3 -o $PWD/$outdir/genie3_scores.sf -F m -i $PWD/$outdir/genie3_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the Tigress scores 'slow'\n"
tigress -O $thread -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/tigress_scores.tsv
seidr import -r -z -n TIGRESS -o $PWD/$outdir/tigress_scores.sf -F m -i $PWD/$outdir/tigress_scores.tsv -g $PWD/$outdir/genes.txt
fi

#VERY_SLOW
if [ $depth == 'VERY_SLOW' ]
then

printf "Calculating the El-Ensemble scores 'very slow'.\n"
el-ensemble -O $thread -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/elnet_scores.tsv
seidr import -r -z -n ELNET -o $PWD/$outdir/elnet_scores.sf -F m -i $PWD/$outdir/elnet_scores.tsv -g $PWD/$outdir/genes.txt
fi

else

#Running the Algorithms in targeted mode

# FAST
printf "Calculating the Pearson Correlation scores 'fast'.\n"
correlation -t $target -m pearson -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/pearson_scores.tsv
seidr import -A -r -u -n PEARSON -o $PWD/$outdir/pearson_scores.sf -F el -i $PWD/$outdir/pearson_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the Spearman Correlation scores 'fast'.\n"
correlation -t $target -m spearman -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/spearman_scores.tsv
seidr import -A -r -u -n SPEARMAN -o $PWD/$outdir/spearman_scores.sf -F el -i $PWD/$outdir/spearman_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the PCorrelation scores 'fast'.\n"
pcor -t $target -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/pcor_scores.tsv
seidr import -A -r -u -n PCOR -o $PWD/$outdir/pcor_scores.sf -F el -i $PWD/$outdir/pcor_scores.tsv -g $PWD/$outdir/genes.txt

#MEDIUM
if [ $depth == 'MEDIUM' ] || [ $depth == 'SLOW' ] || [ $depth == 'VERY_SLOW' ]
then

printf "Calculating the RAW scores 'medium'.\n"
mi -t $target -m RAW -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -M $PWD/$outdir/mi_full_scores.tsv -o $PWD/$outdir/mi_scores.tsv
seidr import -r -u -n MI -o $PWD/$outdir/mi_scores.sf -F el -i $PWD/$outdir/mi_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the CLR scores 'medium'.\n"
mi -t $target -m CLR -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -M $PWD/$outdir/mi_full_scores.tsv -o $PWD/$outdir/clr_scores.tsv
seidr import -r -u -z -n CLR -o $PWD/$outdir/clr_scores.sf -F el -i $PWD/$outdir/clr_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the ARACNE scores 'medium'.\n"
mi -t $target -m ARACNE -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -M $PWD/$outdir/mi_full_scores.tsv -o $PWD/$outdir/aracne_scores.tsv
seidr import -r -u -z -n ARACNE -o $PWD/$outdir/aracne_scores.sf -F el -i $PWD/$outdir/aracne_scores.tsv -g $PWD/$outdir/genes.txt
fi

#SLOW
if [ $depth == 'SLOW' ] || [ $depth == 'VERY_SLOW' ]
then

printf "Calculating the Narromi scores 'slow'\n"
narromi -t $target -O $thread -m interior-point -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/narromi_scores.tsv
seidr import -r -z -n NARROMI -o $PWD/$outdir/narromi_scores.sf -F el -i $PWD/$outdir/narromi_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the Plsnet scores 'slow'\n"
plsnet -t $target -O $thread -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/plsnet_scores.tsv
seidr import -r -z -n PLSNET -o $PWD/$outdir/plsnet_scores.sf -F el -i $PWD/$outdir/plsnet_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the LLR-Ensemble scores 'slow'\n"
llr-ensemble -t $target -O $thread -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/llr_scores.tsv
seidr import -r -z -n LLR -o $PWD/$outdir/llr_scores.sf -F el -i $PWD/$outdir/llr_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the SVM-Ensemble scores 'slow'\n"
svm-ensemble -t $target -O $thread -k POLY -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/svm_scores.tsv
seidr import -r -z -n SVM -o $PWD/$outdir/svm_scores.sf -F el -i $PWD/$outdir/svm_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the Genie3 scores 'slow'.\n"
genie3 -t $target -O $thread -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/genie3_scores.tsv
seidr import -r -z -n GENIE3 -o $PWD/$outdir/genie3_scores.sf -F el -i $PWD/$outdir/genie3_scores.tsv -g $PWD/$outdir/genes.txt
printf "Calculating the Tigress scores 'slow'\n"
tigress -t $target -O $thread -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/tigress_scores.tsv
seidr import -r -z -n TIGRESS -o $PWD/$outdir/tigress_scores.sf -F el -i $PWD/$outdir/tigress_scores.tsv -g $PWD/$outdir/genes.txt
fi

#VERY_SLOW
if [ $depth == 'VERY_SLOW' ]
then

printf "Calculating the El-Ensemble scores 'very slow'.\n"
el-ensemble -t $target -O $thread -i $PWD/$outdir/expression.tsv -g $PWD/$outdir/genes.txt -o $PWD/$outdir/elnet_scores.tsv
seidr import -r -z -n ELNET -o $PWD/$outdir/elnet_scores.sf -F el -i $PWD/$outdir/elnet_scores.tsv -g $PWD/$outdir/genes.txt
fi

fi

#Aggregating the Crowd Network

#FAST
if [ $depth == 'FAST' ]; then seidr aggregate -o $PWD/$outdir/network.sf -m $aggregate $PWD/$outdir/pcor_scores.sf $PWD/$outdir/pearson_scores.sf $PWD/$outdir/spearman_scores.sf; fi

#MEDIUM
if [ $depth == 'MEDIUM' ]; then seidr aggregate -o $PWD/$outdir/network.sf -m $aggregate $PWD/$outdir/aracne_scores.sf $PWD/$outdir/clr_scores.sf $PWD/$outdir/mi_scores.sf $PWD/$outdir/pcor_scores.sf $PWD/$outdir/pearson_scores.sf $PWD/$outdir/spearman_scores.sf; fi

#SLOW
if [ $depth == 'SLOW' ]; then seidr aggregate -o $PWD/$outdir/network.sf -m $aggregate $PWD/$outdir/aracne_scores.sf $PWD/$outdir/clr_scores.sf $PWD/$outdir/genie3_scores.sf $PWD/$outdir/llr_scores.sf $PWD/$outdir/mi_scores.sf $PWD/$outdir/narromi_scores.sf $PWD/$outdir/pcor_scores.sf $PWD/$outdir/pearson_scores.sf $PWD/$outdir/plsnet_scores.sf $PWD/$outdir/spearman_scores.sf $PWD/$outdir/svm_scores.sf $PWD/$outdir/tigress_scores.sf; fi

#VERY SLOW
if [ $depth == 'VERY_SLOW' ]; then seidr aggregate -o $PWD/$outdir/network.sf -m $aggregate $PWD/$outdir/aracne_scores.sf $PWD/$outdir/clr_scores.sf $PWD/$outdir/elnet_scores.sf $PWD/$outdir/genie3_scores.sf $PWD/$outdir/llr_scores.sf $PWD/$outdir/mi_scores.sf $PWD/$outdir/narromi_scores.sf $PWD/$outdir/pcor_scores.sf $PWD/$outdir/pearson_scores.sf $PWD/$outdir/plsnet_scores.sf $PWD/$outdir/spearman_scores.sf $PWD/$outdir/svm_scores.sf $PWD/$outdir/tigress_scores.sf; fi

#Pruning noise from the network, dropping edges bellow a specific P-value (1.28 or 0.10% by default)
#If the desired P-value is 0.05% the -p should be 1.64, and for 0.01% the -p should be 2.32
seidr backbone -F $pvalue $PWD/$outdir/network.sf

#Network Post-Processing for Visualization in Cytoscape
#seidr view --column-headers $PWD/$outdir/network.bb.sf > cyto1_network.txt
seidr view --column-headers -o cyto1_network.txt $PWD/$outdir/network.bb.sf
sed 's/;/\t/g' cyto1_network.txt > cyto2_network.txt
#To discard the NC_Scores use:
#awk '{print $1,$2,$(NF-3)}' cyto2_network.txt > cyto3_network.txt
#To keep the NC_Scores use:
awk '{print $1,$2,$(NF-3),$(NF-1),$(NF)}' cyto2_network.txt > $PWD/$outdir/${outdir}_network.txt

#Applying a cutoof to the network score with network_cut.py if -cutoof value is specified
if [ $cutoff != 'None' ]; then network_cut.py --network $PWD/$outdir/${outdir}_network.txt --cutoff $cutoof --method ${aggregate}_score --out $PWD/$outdir/${outdir}_network_${cutoof}.txt; fi

#Option to keep the graph and node statistics into a *name_of_outdir*_network_stats.txt and *name_of_outdir*_network_nodestats.txt
if [ $stats == 'TRUE' ]
then
	
#seidr reheader $PWD/$outdir/network.bb.sf
#seidr graphstats -o $PWD/$outdir/${outdir}_network_stats.txt $PWD/$outdir/network.bb.sf
seidr stats --metrics PR,BTW,CLO $PWD/$outdir/network.bb.sf
seidr view --centrality -c -o $PWD/$outdir/${outdir}_network_nodestats.txt $PWD/$outdir/network.bb.sf
fi

#Option to collect the correlation scores in a separate file called *name*_network_corscores.txt.
#This option can be used to import the scores in Cytoscape and color the edges of the network, giving a sense of positive and negative correlation of the network.
if [ $edgecorscores == 'TRUE' ] 
then

if [ $depth == 'FAST' ]
then
        awk '{print $1,$2,$4,$6,$8}' cyto2_network.txt > cor1_scores.txt
	#awk '{print $1,$2,$4,$5,$6 }' cyto2_network.txt > cor1_scores.txt
fi

if [ $depth == 'MEDIUM' ]
then
	awk '{print $1,$2,$7,$8,$9}' cyto2_network.txt > cor1_scores.txt
fi

if [ $depth == "SLOW" ]
then
	awk '{print $1,$2,$16,$18,$22}' cyto2_network.txt > cor1_scores.txt
	#awk '{print $1,$2,$10,$11,$14 }' cyto2_network.txt > cor1_scores.txt
fi

if [ $depth == "VERY_SLOW" ]
then
	awk '{print $1,$2,$11,$12,$14 }' cyto2_network.txt > cor1_scores.txt
fi

sed 's/ /\t/g' cor1_scores.txt > cor2_scores.txt

#Only needed to make sure that the correlations are all positive or all negative for the means
#sed "s/;[[:digit:]][[:digit:]]*\.[[:digit:]]//g" cor1_scores.txt > cor2_scores.txt
#sed "s/;[[:digit:]][[:digit:]]*//g" cor2_scores.txt > cor3_scores.txt
#echo "Source Target PCOR_score PEARSON_score SPEARMAN_score" > cor4_scores.txt
#grep -e "[[:space:]]0.*[[:space:]]0.*[[:space:]]0" cor3_scores.txt >> cor4_scores.txt
#grep -e "-0.*-0.*-0" cor3_scores.txt >> cor4_scores.txt

#This script calculates the means of the correlation scores in cor1_scores and makes an extra key collumn to allow the import in Cytoscape
network_scores.py --cor_scores cor2_scores.txt --out $PWD/$outdir/${outdir}_network_corscores.txt
fi

#Lastly, eliminate intermediate files.
#rm cyto1_network.txt
#rm cyto2_network.txt
#rm cor1_scores.txt
#rm cor2_scores.txt


#This python script applies a cutoof to the network by their correlation scores if specified
#if [ "$cutoof"  != "None" ]; then network_cut.py --network cork_cor_scores_mean.txt --cutoff $cutoof --out cork_cor_scores_mean_cut.txt; else echo "Applying  #no cutoff to the correlation scores"; fi

printf "\n\n | -------------------------- END OF SEIDR_SCRIPT ----------------------------- |\n"
