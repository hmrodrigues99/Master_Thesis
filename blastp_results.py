#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
blastp result post-processing
"""

import pandas as pd
import argparse


def blastp_results(blastp_results, outname):
    df = pd.read_table(blastp_results, sep='\t')
    df.columns = ['Quercus_gene', 'Arabidopsis_gene', 'alignment_len', 'query_len', 'subject_len', 'query_start', 'query_end', 'subject_s', 'subject_e', 'evalue']
    df_unduplicated = df.drop_duplicates(subset=['Quercus_gene'], keep='first')
    df_link = df_unduplicated[['Quercus_gene', 'Arabidopsis_gene']]
    #print(df.shape)
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
