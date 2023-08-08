__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "map reads with bowtie/bbmap"

import re
import os
import glob
import subprocess
import scqv.extra_functions as ef
from sys import platform
import pysam

PATH = os.path.dirname(os.path.abspath(__file__))
PATH = re.sub("map_reads", "", PATH)
if platform == "linux":
    bbmap2path = os.path.join(PATH, "aligntools", "bbmapv38.9")
    bowtie2path = os.path.join(PATH, "aligntools","bowtie2-2.4.2-linux-x86_64")
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
            if len(files_genome)==0:
                files_genome = glob.glob(os.path.join(args.path2genome, "*.fa"))
                if (len(files_genome) > 0):
                    bowtie2genomefile = files_genome[0]
                else: 
                    files_genome = glob.glob(os.path.join(args.path2genome, "*.fasta"))
                    if (len(files_genome) > 0):
                        bowtie2genomefile = files_genome[0]            
                    else:
                        files_genome = glob.glob(os.path.join(args.path2genome, "*.fna"))
                        if (len(files_genome) > 0):
                            bowtie2genomefile = files_genome[0]            
            if len(files_genome) > 1:
                print ("WARNING:  two fasta files in genome folder.")
            arg=[os.path.join(bowtie2path, "bowtie2-build"), bowtie2genomefile, os.path.join(args.path2genome, "genome")]
            ef._run_subprocesses(arg, "STATUS: generating libraries bowtie2 ...", "extracting unmapped")

        else:
            print("STATUS: bowtie2 indexes have already been made for this genome")
    elif args.aligner == "bbmap":
        files_genome = glob.glob(os.path.join(args.path2genome, "*.fa"))
        if (len(files_genome) > 0):
            bbmapgenomefile = files_genome[0]
        else: 
            files_genome = glob.glob(os.path.join(args.path2genome, "*.fasta"))
            if (len(files_genome) > 0):
                bbmapgenomefile = files_genome[0]            
            else:
                files_genome = glob.glob(os.path.join(args.path2genome, "*.fna"))
                if (len(files_genome) > 0):
                    bbmapgenomefile = files_genome[0]            
    else:
        raise ValueError(args.aligner+" not an available aligner specify bbmap or bowtie2, all lowercase")

    # -- align reads 
    if args.aligner == "bowtie2":
        if args.bowtie2_params is None:
            arg=[os.path.join(bowtie2path, "bowtie2"), "-p", args.processors, 
            "-q", os.path.join(args.output_path,"_tmp", "unmapped.fq"), 
            "-x", os.path.join(args.path2genome,"genome"), 
            "-S", os.path.join(args.output_path, "virus_al.sam")]
        else:
            arg=[os.path.join(bowtie2path, "bowtie2"), "-p", args.processors, 
            "-q", os.path.join(args.output_path,"_tmp", "unmapped.fq"), 
            "-x", os.path.join(args.path2genome,"genome"), 
            "-S", os.path.join(args.output_path, "virus_al.sam")] + args.bowtie2_params
            
        ef._run_subprocesses(arg, "STATUS: Align reads bowtie2 global...", "aligning reads")
    elif args.aligner == "bbmap":
        if args.bbmap_params is None:
            arg=[os.path.join(bbmap2path, "bbmap.sh"), "ref="+bbmapgenomefile, "in="+os.path.join(args.output_path,"_tmp", "unmapped.fq"), "nodisk=t" ,"local=f",
                "covstats="+os.path.join(args.output_path, "constats.txt"), "covhist="+os.path.join(args.output_path,"covhist.txt"), "minid=0.95",
                "basecov="+os.path.join(args.output_path,"basecov.txt"), "bincov="+os.path.join(args.output_path, "bincov.txt"),
                "out="+os.path.join(args.output_path, "virus_al.sam")]
        else:
            arg=[os.path.join(bbmap2path, "bbmap.sh"), "ref="+bbmapgenomefile, "in="+os.path.join(args.output_path,"_tmp", "unmapped.fq"), "nodisk=t" ,"local=f",
                "covstats="+os.path.join(args.output_path, "constats.txt"), "covhist="+os.path.join(args.output_path,"covhist.txt"), "minid=0.95",
                "basecov="+os.path.join(args.output_path,"basecov.txt"), "bincov="+os.path.join(args.output_path, "bincov.txt"),
                "out="+os.path.join(args.output_path, "virus_al.sam")] + args.bbmap_params
        ef._run_subprocesses(arg, "STATUS: Align reads bbmap global...", "aligning reads")

    # --generate necessary output files
    print("STATUS: Extracting mapped reads and converting to bam")
    pysam.view("-@", str(args.processors),"-S", "-b", os.path.join(args.output_path,"virus_al.sam"), "-o", os.path.join(args.output_path, "virus_al.bam"), catch_stdout=False)

    print("STATUS: Converting mapped reads to from bam to sam file")
    pysam.view("-@", str(args.processors), "-h", os.path.join(args.output_path, "virus_al.bam"), "-o",  os.path.join(args.output_path, "virus_al_mapped.sam"), catch_stdout=False)

    print("STATUS: Sorting bam")
    pysam.sort("-@", str(args.processors), os.path.join(args.output_path,"virus_al.bam"), "-o",  os.path.join(args.output_path,"virus_al_sort.bam"))

    print( "STATUS: generating bam index")
    pysam.index("-@", str(args.processors), os.path.join(args.output_path,"virus_al_sort.bam"))