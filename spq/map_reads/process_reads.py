__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "processes unmapped reads"

import re
import os
from sys import platform
import pysam
PATH = os.path.dirname(os.path.abspath(__file__))
PATH = re.sub("map_reads", "", PATH)
if platform == "linux":
    bowtie2path = os.path.join(PATH, "aligntools", "bowtie2-2.4.2-linux-x86_64")
elif platform == "OS":
    bowtie2path = os.path.join(PATH, "aligntools","bowtie2-2.4.2-macos-x86_64")
else:
    ValueError("Program wont run on this operating system "+platform)

def process_unmapped_reads(args, viable_cb):
    try: 
       os.mkdir(args.output_path)
    except:
        pass 
    try: 
        print("STATUS: Generating _tmp/ files")
        os.mkdir(os.path.join(args.output_path, "_tmp"))
        # -- pull out unmapped reads

        pysam.view( "-@", str(args.processors), "-b", "-h", "-f", "4", os.path.join(args.path10x, "outs", "possorted_genome_bam.bam"), "-o", os.path.join(args.output_path, "_tmp","unmapped.bam"), catch_stdout=False)
        pysam.view("-@", str(args.processors), "-h", os.path.join(args.output_path,"_tmp", "unmapped.bam"), "-o", os.path.join(args.output_path, "_tmp", "unmapped.sam"), catch_stdout=False)

        # -- pull out umi and cell barcode information 
        print ("STATUS: extracting cell barcodes and UMIs...")
        with open(os.path.join(args.output_path, "_tmp", "barcode_umi_read_table.csv"), "w") as fout:
            fout.write("cell_barcode,umi,read\n")
            with open(os.path.join(args.output_path,"_tmp","unmapped.sam")) as fin:
                for line in fin:
                    if not line.startswith("@"):
                        line = line.strip()
                        larray = line.split("\t")
                        if larray[-1].startswith("UB") and larray[-4].startswith("CB"):
                            if re.sub("CB:Z:", "", larray[-4]) in viable_cb: 
                                fout.write(re.sub("CB:Z:", "", larray[-4])+","+re.sub("UB:Z:", "", larray[-1])+","+larray[0]+"\n")

        # -- convert to fastq file
        pysam.bam2fq('-n', '-0', os.path.join(args.output_path, "_tmp", "unmapped.fq.gz"), os.path.join(args.output_path,"_tmp","unmapped.bam"),  catch_stdout=False)

    except:
        print("STATUS: _tmp folder already exists so not regenerating data, path "+str(os.path.join(args.output_path, "_tmp")))
        pass