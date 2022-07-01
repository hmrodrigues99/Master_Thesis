# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 11:16:46 2022

@author: Hugo Rodrigues
"""

#Receive list of TFs, check if they are DEGs in Cork tissue, check number of targets in the network and how many of those are either up-expressed or down-expressed
#This script receives as input a TF_list, my Cork network and the lists of DEG identified in the 3 tissues (cork, phloem and xylem) with the respective log2FoldChange values.
#This was used to generate table X to Y, by receiving different TF lists: a custom list of TFs which had validated connections in my network by ConnecTF (table X1) and the TF lists identified by PlantTFDB within relevant network modules (Table Y2, Y4...)

import pandas as pd

#input TF list
#TFs with atleast one validated network tf-target interaction by ConnecTF
TF_list = pd.read_csv("C:/Users/Hugo Rodrigues/Documents/TF_Targets/Cork_network_TFs.txt", sep="\t")
TF_list = TF_list["TFs"]

#TF list of Module 1
TFs_Module1 = pd.read_csv("C:/Users/Hugo Rodrigues/Documents/Thesis Material/Module1.csv")
TFs_Module1 = TFs_Module1.dropna(subset = ["Family"])
TF_list = TFs_Module1["name"]

#TF list of Module 2
TFs_Module2 = pd.read_csv("C:/Users/Hugo Rodrigues/Documents/Thesis Material/Module2.csv")
TFs_Module2 = TFs_Module2.dropna(subset = ["Family"])
TF_list = TFs_Module2["name"]

#TF list of Module 4
TFs_Module4 = pd.read_csv("C:/Users/Hugo Rodrigues/Documents/Thesis Material/Module4.csv")
TFs_Module4 = TFs_Module4.dropna(subset = ["Family"])
TF_list = TFs_Module4["name"]


#input Cork network
network = pd.read_table("C:/Users/Hugo Rodrigues/Documents/Thesis Material/Task_2/Cork_network.txt")
#print(network)

#input DEGs lists
deg_all = pd.read_table("C:/Users/Hugo Rodrigues/Documents/Cytoscape Networks/DEG_All.txt", header=0, names=["DEGs"])
deg_cork = pd.read_table("C:/Users/Hugo Rodrigues/Documents/Mercator/Locs_Cork_Log2Fold.txt")
deg_phloem = pd.read_table("C:/Users/Hugo Rodrigues/Documents/Mercator/Locs_Phloem_Log2Fold.txt")
deg_xylem  = pd.read_table("C:/Users/Hugo Rodrigues/Documents/Mercator/Locs_Xylem_Log2Fold.txt")

#Final Output dataframe collumns
is_deg = []
is_deg_cork = []
nTargets = []
nTargetsDEG = []
nTargetsDEG_p = []
nTargetsDEG_p_per = []
nTargetsDEG_n = []
nTargetsDEG_n_per = []
nTargetsDEG_Cork = []
nTargetsDEG_Cork_p = []
nTargetsDEG_Cork_p_per = []
nTargetsDEG_Cork_n = []
nTargetsDEG_Cork_n_per = []

for TF in TF_list:
    Targets = []
    filtered_network = network[network["LOC1"].str.contains(TF) | network["LOC2"].str.contains(TF)]
    for LOC in filtered_network["LOC1"]:
        if not LOC == TF:
            Targets.append(LOC)
    for LOC in filtered_network["LOC2"]:
        if not LOC == TF:
            Targets.append(LOC)
    ntargets = len(Targets)
    nTargets.append(ntargets)
    #Check if a TF is a deg in drought and if was also identified in Cork tissue
    check_deg = deg_all[deg_all["DEGs"].str.contains(TF)]
    if len(check_deg) != 0:
        is_deg.append("Yes")
        check_deg_cork = deg_cork[deg_cork["Gene"].str.contains(TF)]
        if len(check_deg_cork) != 0:
            is_deg_cork.append("Yes")
        else:
            is_deg_cork.append("No")
    else:
        is_deg.append("No")
        is_deg_cork.append("No")
    
    #Checks if the Target is a deg (if it is present in the Deg_all list), and if it is, appends it to the new Targets_degs list
    # Also check if its a deg target identified in cork (it isnt necessary to be a deg exclusive to cork)
    Targets_degs = []
    Targets_degs_cork = []
    for target in Targets:
        check_deg = deg_all[deg_all["DEGs"].str.contains(target)]
        if len(check_deg) != 0:
            Targets_degs.append(target)
            check_deg = deg_cork[deg_cork["Gene"].str.contains(target)]
            if len(check_deg) != 0:
                Targets_degs_cork.append(target)
    ntargetsDEG = len(Targets_degs)
    nTargetsDEG.append(ntargetsDEG)
    ntargetsDEG_Cork = len(Targets_degs_cork)
    nTargetsDEG_Cork.append(ntargetsDEG_Cork)
    
    #Check number of positive and negative deg targets in cork, xylem and phloem tissues
    positive_deg = []
    negative_deg = []
    positive_deg_cork = []
    negative_deg_cork = []
    for target in Targets_degs:
        target_log2fc_phloem = deg_phloem[deg_phloem["Gene"] == target]
        target_log2fc_xylem = deg_xylem[deg_xylem["Gene"] == target]
        target_log2fc_cork = deg_cork[deg_cork["Gene"] == target]
        if len(target_log2fc_phloem) != 0: 
            target_log2fc_phloem = target_log2fc_phloem.iloc[0,1]
            if target_log2fc_phloem > 0:
                positive_deg.append(target)
            else:
                negative_deg.append(target)
        if len(target_log2fc_xylem) != 0:
            target_log2fc_xylem = target_log2fc_xylem.iloc[0,1]
            if target_log2fc_xylem > 0:
                positive_deg.append(target)
            else:
                negative_deg.append(target)
        if len(target_log2fc_cork) != 0:
            target_log2fc_cork = target_log2fc_cork.iloc[0,1]
            if target_log2fc_cork > 0:
                positive_deg_cork.append(target)
            else:
                negative_deg_cork.append(target)
    #drop list duplicates (both in positive and negative deg lists)
    positive_deg = list(set(positive_deg))
    negative_deg = list(set(negative_deg))
    
    #Get number of all positive and negative deg targets in both lists and get respective %
    for target in positive_deg_cork:
        if target not in positive_deg:
            positive_deg.append(target)
    ntargetsDEG_p = len(positive_deg)
    nTargetsDEG_p.append(ntargetsDEG_p)
    if ntargetsDEG != 0:
        ntargetsDEG_p_per = ntargetsDEG_p/ntargetsDEG*100
    else:
        ntargetsDEG_p_per = "-"
    nTargetsDEG_p_per.append(ntargetsDEG_p_per)
    for target in negative_deg_cork:
        if target not in negative_deg:
            negative_deg.append(target)
    ntargetsDEG_n = len(negative_deg)
    nTargetsDEG_n.append(ntargetsDEG_n)
    if ntargetsDEG != 0:
        ntargetsDEG_n_per = ntargetsDEG_n/ntargetsDEG*100
    else:
        ntargetsDEG_n_per = "-"
    nTargetsDEG_n_per.append(ntargetsDEG_n_per)
    
    #Check number of positive and negative deg_cork targets and respective %
    ntargetsDEG_Cork_p = len(positive_deg_cork)
    nTargetsDEG_Cork_p.append(ntargetsDEG_Cork_p)
    ntargetsDEG_Cork_n = len(negative_deg_cork)
    nTargetsDEG_Cork_n.append(ntargetsDEG_Cork_n)
    if ntargetsDEG_Cork != 0:    
        ntargetsDEG_Cork_p_per = ntargetsDEG_Cork_p/ntargetsDEG_Cork*100
        ntargetsDEG_Cork_n_per = ntargetsDEG_Cork_n/ntargetsDEG_Cork*100
    else:
        ntargetsDEG_Cork_p_per = "-"
        ntargetsDEG_Cork_n_per = "-"
    nTargetsDEG_Cork_p_per.append(ntargetsDEG_Cork_p_per)
    nTargetsDEG_Cork_n_per.append(ntargetsDEG_Cork_n_per)

print(len(TF_list))
print(len(is_deg))
print(len(is_deg_cork))

#Output
data = {"TF": TF_list, "DEG": is_deg, "DEG_Cork": is_deg_cork, "Number of Targets": nTargets, "DEG Targets": nTargetsDEG, "pDEG Targets": nTargetsDEG_p, "pDEG Targets per": nTargetsDEG_p_per, "nDEG Targets": nTargetsDEG_n, "nDEG Targets per": nTargetsDEG_n_per, "DEG_Cork Targets": nTargetsDEG_Cork, "pDEG_Cork Targets": nTargetsDEG_Cork_p, "pDEG_Cork Targets per": nTargetsDEG_Cork_p_per, "nDEG_Cork Targets": nTargetsDEG_Cork_n, "nDEG_Cork Targets per": nTargetsDEG_Cork_n_per}
df = pd.DataFrame(data)
df = df.sort_values(by=["Number of Targets", "DEG Targets", "DEG_Cork Targets"], ascending=False)
df.to_csv("C:/Users/Hugo Rodrigues/Documents/Thesis Material/data.csv", index=False)
