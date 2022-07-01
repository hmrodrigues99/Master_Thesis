#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script load an network.txt and applies a cut to the global scores
"""

import pandas as pd
import argparse

def network_cut(network, cutoff, method_score, outname):
    df = pd.read_table(network, sep=" ")
    cut_max = cutoff
    cut_min = cutoff * -1
    df_cut = df[~df[method_score].between(cut_min, cut_max)]
    df_cut.to_csv(outname, sep=" ", index=False)
    
def main():
    parser=argparse.ArgumentParser(description="Apply a cutoff to the irp_score column of a network.txt")
    parser.add_argument('--network', required=True, metavar='tab',
                        help='network where the cutoff should be applied')
    parser.add_argument('--cutoff', required=True, type=float, metavar="float",
                        help='cutoff value')
    parser.add_argument('--method_score', type=str, metavar="STR", default="cor_mean",
                        help='If it is a Crowd Network please specify the method in which it was agregated (example: irp_score)')
    parser.add_argument('--out', type=str, metavar="STR", default="network_cut.txt",
                        help='Output file name')
    
    args = parser.parse_args()
    
    network_cut(args.network, args.cutoff, args.method_score, args.out)
    

main()
