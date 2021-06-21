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
    bbmap2path = os.path.join(PATH, "aligntools", "bbmapv38.9")
    bowtie2path = os.path.join(PATH, "aligntools","bowtie2-2.4.2-linux-x86_64")
    samtoolspath = os.path.join(PATH, "extra_tools", "samtoolsv1.11_linux","bin")
elif platform == "OS":
    bowtie2path = os.path.join(PATH, "aligntools", "bowtie2-2.4.2-macos-x86_64")
else:
    ValueError("Program wont run on this operating system "+platform)

def map2viralgenome(args):

    # -- check if libraries are present if not generate 
    if args.aligner == "bowtie2":
        files = glob.glob(os.path.join(args.path2genome,"*.bt2"))
        if len(files) == 0:
            files_genome = glob.glob(os.path.join(args.path2genome, "*.fa"))
            if len(files_genome) > 1:
                print ("WARNING:  two fasta files in genome folder.")
            arg=[bowtie2path+"bowtie2-build", files_genome[0], os.path.join(args.path2genome, "genome")]
            ef._run_subprocesses(arg, "STATUS: generating libraries bowtie2 ...", "extracting unmapped")

        else:
            print("STATUS: bowtie2 indexes have already been made for this genome")
    elif args.aligner == "bbmap":
        if os.path.isdir(os.path.join(args.path2genome, "ref")) is False:
            files_genome = glob.glob(os.path.join(args.path2genome, "*.fa"))
            bbmapgenomefile = files_genome[0]
    else:
        raise ValueError(args.aligner+" not an available aligner specify bbmap or bowtie2")

    # -- align reads 
    if args.alignment_type == "local":
        if args.aligner == "bowtie2":
            arg=[os.path.join(bowtie2path, "bowtie2"), "-p", args.processors, "-q", os.path.join(args.output_path,"_tmp", "unmapped.fq"), "--local", "-x", os.path.join(args.path2genome,"genome"), "-S", os.path.join(args.output_path, "virus_al.sam")]
            ef._run_subprocesses(arg, "STATUS: Align reads bowtie2 local...", "aligning reads")
        elif arg.aligner == "bbmap":
            arg=[os.path.join(bbmap2path, "bbmap.sh"), "ref="+bbmapgenomefile, "in="+os.path.join(args.output_path,"_tmp", "unmapped.fq"), "nodisk=t" ,"local=t",
                "covstats="+os.path.join(args.output_path, "constats.txt"), "covhist="+os.path.join(args.output_path,"covhist.txt"),
                "basecov=+"+os.path.join(args.output_path,"basecov.txt"), "bincov="+os.path.join(args.output_path, "bincov.txt"),
                "out="+os.path.join(args.output_path, "virus_al.sam")]
            ef._run_subprocesses(arg, "STATUS: Align reads bbmap local...", "aligning reads")
    else:
        if args.aligner == "bowtie2":
            arg=[os.path.join(bowtie2path, "bowtie2"), "-p", args.processors, "-q", os.path.join(args.output_path,"_tmp", "unmapped.fq"), "-x", os.path.join(args.path2genome,"genome"), "-S", os.path.join(args.output_path, "virus_al.sam")]
            ef._run_subprocesses(arg, "STATUS: Align reads bowtie2 global...", "aligning reads")
        elif arg.aligner == "bbmap":
            arg=[os.path.join(bbmap2path, "bbmap.sh"), "ref="+bbmapgenomefile, "in="+os.path.join(args.output_path,"_tmp", "unmapped.fq"), "nodisk=t" ,"local=f",
                "covstats="+os.path.join(args.output_path, "constats.txt"), "covhist="+os.path.join(args.output_path,"covhist.txt"),
                "basecov=+"+os.path.join(args.output_path,"basecov.txt"), "bincov="+os.path.join(args.output_path, "bincov.txt"),
                "out="+os.path.join(args.output_path, "virus_al.sam")]
            ef._run_subprocesses(arg, "STATUS: Align reads bbmap global...", "aligning reads")

    # --generate necessary output files
    arg=[samtoolspath+"samtools", "view", "-@", args.processors, "-F", "4", "-Sb", os.path.join(args.output_path,"virus_al.sam"), "-o", os.path.join(args.output_path, "virus_al.bam")]
    ef._run_subprocesses(arg, "STATUS: Extracting mapped reads and converting to bam", "extracting mapped reads and converting to bam")

    arg=[samtoolspath+"samtools", "sort", os.path.join(args.output_path,"virus_al.bam"), "-o", os.path.join(args.output_path,"virus_al_sort.bam")]
    ef._run_subprocesses(arg, "STATUS: Sorting bam", "sorting bam")

    arg=[samtoolspath+"samtools", "index", os.path.join(args.output_path,"virus_al_sort.bam"), os.path.join(args.output_path,"virus_al_sort.bam.bai")]
    ef._run_subprocesses(arg, "STATUS: generating bam index", "generating bam index")


