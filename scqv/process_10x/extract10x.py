__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "pull out necessary 10x cells"

import os
import re
import gzip

def extract_viable_10x(path10x):
    print ("STATUS: pull out viable cells")
    barcodes=[]
    if os.path.isdir(path10x):
        with gzip.open(path10x+'/outs/filtered_feature_bc_matrix/barcodes.tsv.gz', 'rb') as f:
            file_content = f.read().decode("UTF-8").split("\n")
            if file_content[-1] =='': 
                del file_content[-1]
            barcodes=[ line for line in file_content]
        print ("STATUS: Number of viable barcodes is "+str(len(barcodes))+"... should match what is in 10x web_summary")
        return barcodes
    else:
        AssertionError("Could not find 10x path "+path10x+"/outs/")
