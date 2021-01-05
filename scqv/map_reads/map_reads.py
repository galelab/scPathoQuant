__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "map reads with bowtie (currently only version)"

import re
import os
import glob
import subprocess
import scqv.extra_functions as ef
from sys import platform

PATH = os.path.dirname(os.path.abspath(__file__))
PATH = re.sub("map_reads", "", PATH)
if platform == "linux":
    bowtie2path = os.path.join(PATH, "aligntools/bowtie2-2.4.2-linux-x86_64/")
    samtoolspath = os.path.join(PATH, "extra_tools/samtoolsv1.11_linux/bin/")
elif platform == "OS":
    bowtie2path = os.path.join(PATH, "aligntools/bowtie2-2.4.2-macos-x86_64/")
else:
    ValueError("Program wont run on this operating system "+platform)

def map2viralgenome(args):

    # -- check if libraries are present if not generate 
    files = glob.glob(os.path.join(args.path2genome,"*.bt2"))
    if len(files) == 0:
        files_genome = glob.glob(os.path.join(args.path2genome, "*.fa"))
        if len(files_genome) > 1:
            print ("WARNING:  two fasta files in genome folder.")
        arg=[bowtie2path+"bowtie2-build", files_genome[0], os.path.join(args.path2genome, "genome")]
        ef._run_subprocesses(arg, "STATUS: generating libraries ...", "extracting unmapped")

    else:
        print("STATUS: libraries have already been made for this genome")

    # -- align reads 
    arg=[bowtie2path+"bowtie2", "-p", args.processors, "-q", os.path.join(args.output_path,"_tmp/unmapped.fq"), "-x", os.path.join(args.path2genome,"genome"), "-S", os.path.join(args.output_path, "virus_al.sam")]
    ef._run_subprocesses(arg, "STATUS: Align reads ...", "aligning reads")
    
    # --generate necessary output files
    arg=[samtoolspath+"samtools", "view", "-@", args.processors, "-F", "4", "-Sb", os.path.join(args.output_path,"virus_al.sam"), "-o", os.path.join(args.output_path, "virus_al.bam")]
    ef._run_subprocesses(arg, "STATUS: Extracting mapped reads and converting to bam", "extracting mapped reads and converting to bam")

    arg=[samtoolspath+"samtools", "sort", os.path.join(args.output_path,"virus_al.bam"), "-o", os.path.join(args.output_path,"virus_al_sort.bam")]
    ef._run_subprocesses(arg, "STATUS: Sorting bam", "sorting bam")

    arg=[samtoolspath+"samtools", "index", os.path.join(args.output_path,"virus_al_sort.bam"), os.path.join(args.output_path,"virus_al_sort.bam.bai")]
    ef._run_subprocesses(arg, "STATUS: generating bam index", "generating bam index")


