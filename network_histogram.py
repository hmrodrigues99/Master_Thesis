#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script loads a network.txt and creates a histogram with the distribution of the seidr aggregated scores
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse

def network_histogram(network, score, outname):
    df = pd.read_table(network, sep=" ")
    fig, ax= plt.subplots()
    df.hist(column='irp_score' ,ax=ax)
    plt.xlabel('irp score')
    plt.ylabel('number of edges')
    plt.title(outname)
    fig.savefig(outname + '.png')
    
def main():
    parser=argparse.ArgumentParser(description="Make an histogram of the aggregated scores of a network.txt")
    parser.add_argument('--network', required=True, metavar='tab',
                        help='network from which the scores will generate a histogram')
    parser.add_argument('--method_score', type=str, metavar="STR", default="cor_mean",
                        help='If it is a Crowd Network please specify the method in which it was agregated (example: irp_score)')
    parser.add_argument('--out', type=str, metavar="STR", default="network_histogram",
                        help='Output file name')
    
    args = parser.parse_args()
    
    network_histogram(args.network, args.method_score, args.out)
    

main()
