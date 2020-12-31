__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "pull out necessary 10x cells"

import os
import re
import gzip
import subprocess

def extract_viable_10x(path10x):
    print ("STATUS: pull out viable cells")
    barcodes=[]
    if os.path.isdir(path10x):
        with gzip.open(path10x+'/outs/filtered_feature_bc_matrix/barcodes.tsv.gz', 'rb') as f:
            file_content = f.read().decode("UTF-8").split("\n")
            if file_content[-1] =='': 
                del file_content[-1]
            barcodes=[ line for line in file_content]
        return barcodes
    else:
        AssertionError("Could not find 10x path "+path10x+"/outs/")
