#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script receives as input a 1)Co-expression Network from Seidr (seidr_script.sh); the 2)Promoter_sequences_df.csv file, which is a table with the genes and their respective promoter regions (tf_target_promoter_sequences.py); and the 3)ID of a TF, or a list of them.
It checks if the TF Binding Sites (retrieved by get_TFBS.sh) appear/match within the promoter regions of their targets in the network.
As output, it writes a table with TF | Target(s) that matched. A separate list will be created for TFs without TFBS."""

import pandas as pd
import argparse

def match_TFBS_promoter(network, promoter_sequences, TF, outname):
    network = pd.read_csv(network, sep=",")
    query_promoter_sequences = pd.read_csv(promoter_sequences, sep=",")
    target_list=[]
    filtered_network = network[network["LOC1"].str.contains(TF) | network["LOC2"].str.contains(TF)]
    for LOC in filtered_network["LOC1"]:
        if not LOC == TF:
            target_list.append(LOC)
    for LOC in filtered_network["LOC2"]:
        if not LOC == TF:
            target_list.append(LOC)
    print(target_list)
    queried_rows = query_promoter_sequences.loc[query_promoter_sequences["Gene"].isin(target_list)]
    fasta_output = open(outname + "_TF_" + TF + "_fasta.fa", "w")
    for row_index in range(0, len(queried_rows)):
        S = queried_rows["Scafold"].values[row_index]
        G = queried_rows["Gene"].values[row_index]
        seq = queried_rows["Sequence"].values[row_index]
        one_entry = ">" + G + " " + S + "\n" + seq + "\n"
        one_entry = one_entry.upper()
        fasta_output.write(one_entry)
    fasta_output.close()

def main():
    parser=argparse.ArgumentParser(description="This script receives as input a 1)Co-expression Network from Seidr (seidr_script.sh); the 2)Promoter_sequences_df.csv file, which is a table with the genes and their respective promoter regions (tf_target_promoter_sequences.py); and the 3)ID of a TF, or a list of them. It checks if the TF Binding Sites (retrieved by get_TFBS.sh) appear/match within the promoter regions of their targets in the network. As output, it writes a table with TF | Target | TFBS | 'match':YES or NO, and a summary table with TF | number of targets | number of targets that matched.")
    parser.add_argument('--network', required=True, metavar="tab", help="co-expression network from seidr")
    parser.add_argument('--promoter_sequences', required=True, metavar="tab", help="pandas dataframe with gene, scaffold and promoter region sequence columns")
    parser.add_argument('--TF', required=True, metavar="STR", help="One single query TF ID or a list of them (separated by enter/one TF ID per line)")
    parser.add_argument('--out', type=str, metavar="STR", default="TFBS")

    args = parser.parse_args()

    match_TFBS_promoter(args.network, args.promoter_sequences, args.TF, args.out)

main()


