# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 09:31:47 2022

@author: Hugo Rodrigues
"""

from infomap import Infomap

filepath = r"C:\Users\Hugo Rodrigues\Documents\Hugo_Networks\df_IDs_score.txt"
outdir =   r"C:\Users\Hugo Rodrigues\Documents\Hugo_Networks\df_IDs_score"


im = Infomap(silent=True)
im.read_file(filepath)
im.run()

im.write_clu(outdir+".clu")
im.write_flow_tree(outdir+".ftree")


