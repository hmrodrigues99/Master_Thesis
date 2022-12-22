#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Script used to process the blastp results obtained from the quercus_blastp_arabidopsis.sh script.
Keeps only the highest quality (lower evalue) quercus suber - arabidopsis thaliana homology link.
Run arguments used during the Master's work: blastp_results(quercus_blastp_thaliana_link_top20.txt, quercus_arabidopsis)
The output file (e.g. quercus_link_thaliana.txt) should be imported in Cytoscape as a node table.
This will add the Arabidopsis IDs to the respective Quercus suber genes (using as key column: Quercus_gene)
"""

import pandas as pd
import argparse


def blastp_results(blastp_results, outname):
    df = pd.read_table(blastp_results, sep='\t')
    df.columns = ['Quercus_gene', 'Arabidopsis_gene', 'alignment_len', 'query_len', 'subject_len', 'query_start', 'query_end', 'subject_s', 'subject_e', 'evalue']
    df_unduplicated = df.drop_duplicates(subset=['Quercus_gene'], keep='first')
    df_link = df_unduplicated[['Quercus_gene', 'Arabidopsis_gene']]
    df_unduplicated.to_csv(outname+'.txt', sep="\t")
    df_link.to_csv(outname+'_linked.txt', sep="\t")

def main():
    parser=argparse.ArgumentParser(description='process blastp results')
    parser.add_argument('--blastp_results', required=True, metavar='tab', 
                        help='blastp tabular results - format 6')
    parser.add_argument('--out', metavar='STR', type=str, default="quercus_link_arabidopsis.txt",
                        help="Output file name")
    
    args = parser.parse_args()
    
    blastp_results(args.blastp_results, args.out)

    
main()
