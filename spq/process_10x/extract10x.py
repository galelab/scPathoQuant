__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "pull out necessary 10x cells"

import os
import re
import sys
import gzip
from spq.map_reads import map_reads as mr

def checkgenomefile(args):
    files_genome = mr._get_genome_file(args)
    if len(files_genome) > 1:
        print("WARNING: two fasta files in genome folder using this one "+files_genome[0])
    with open(files_genome[0], "r") as fin:
        for line in fin:
            line = line.strip()
            if line.startswith(">"):
              if re.search("\s", line) :
                    print("ERROR: space in line "+line.strip()+"\nin file "+files_genome[0]+"\nRemove ALL spaces and rerun code")
                    sys.exit()

def extract_viable_10x(path10x):
    print ("STATUS: pull out viable cells")
    barcodes=[]
    if os.path.isdir(path10x):
        with gzip.open(os.path.join(path10x, 'outs', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz'), 'rb') as f:
            file_content = f.read().decode("UTF-8").split("\n")
            if file_content[-1] =='': 
                del file_content[-1]
            barcodes=[ line for line in file_content]
        print ("STATUS: Number of viable barcodes is "+str(len(barcodes))+"... should match what is in 10x web_summary")
        return barcodes
    else:
        raise AssertionError("Could not find 10x path "+path10x+"/outs/")
