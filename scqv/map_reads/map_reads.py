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
    star2path = os.path.join(PATH, "aligntools", "STAR-2.7.10a", "bin", "Linux_x86_64_static")
    # app_path = os.path.join(PATH, 'align', 'STAR')
    # os.environ["PATH"] += os.pathsep + app_path
    bowtie2path = os.path.join(PATH, "aligntools","bowtie2-2.4.2-linux-x86_64")
    samtoolspath = os.path.join(PATH, "extra_tools", "samtoolsv1.11_linux","bin")
elif platform == "OS":
    bowtie2path = os.path.join(PATH, "aligntools", "bowtie2-2.4.2-macos-x86_64")
    star2path = os.path.join(PATH, "aligntools", "STAR-2.7.10a", "bin", "MacOSX_x86_64")
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
            arg=[os.path.join(bowtie2path, "bowtie2-build"), files_genome[0], os.path.join(args.path2genome, "genome")]
            ef._run_subprocesses(arg, "STATUS: generating libraries bowtie2 ...", "extracting unmapped")

        else:
            print("STATUS: bowtie2 indexes have already been made for this genome")
    elif args.aligner == "bbmap":
        files_genome = glob.glob(os.path.join(args.path2genome, "*.fa"))
        bbmapgenomefile = files_genome[0]
    elif args.aligner == "star":
        files_genome = glob.glob(os.path.join(args.path2genome, "*.fa"))
        genomefile = files_genome[0]

        arg=[os.path.join(star2path, "STAR"), "--runThreadN", args.processors, "--runMode", "genomeGenerate", 
            "--genomeDir", os.path.join(args.path2genome, "STAR_indicies"),"--genomeFastaFiles", genomefile]
        ef._run_subprocesses(arg, "STATUS: generating indicies for STAR ...", "generating indicies for STAR")
    else:
        raise ValueError(args.aligner+" not an available aligner specify star, bbmap or bowtie2, all lowercase")

    # -- align reads 
    if args.alignment_type == "local":
        if args.aligner == "bowtie2":
            arg=[os.path.join(bowtie2path, "bowtie2"), "-p", args.processors, "-q", os.path.join(args.output_path,"_tmp", "unmapped.fq"), "--local", "-x", os.path.join(args.path2genome,"genome"), "-S", os.path.join(args.output_path, "virus_al.sam")]
            ef._run_subprocesses(arg, "STATUS: Align reads bowtie2 local...", "aligning reads")
        elif args.aligner == "bbmap":
            arg=[os.path.join(bbmap2path, "bbmap.sh"), "ref="+bbmapgenomefile, "in="+os.path.join(args.output_path,"_tmp", "unmapped.fq"), "nodisk=t" ,"local=t",
                "covstats="+os.path.join(args.output_path, "constats.txt"), "covhist="+os.path.join(args.output_path,"covhist.txt"),
                "basecov="+os.path.join(args.output_path,"basecov.txt"), "bincov="+os.path.join(args.output_path, "bincov.txt"),
                "out="+os.path.join(args.output_path, "virus_al.sam")]
            ef._run_subprocesses(arg, "STATUS: Align reads bbmap local...", "aligning reads")
        elif args.aligner == "star":
            raise ValueError(args.aligner+" local alignment doesn't work with star --alignment_type needs to be set to global")

    elif args.alignment_type == "global":
        if args.aligner == "bowtie2":
            arg=[os.path.join(bowtie2path, "bowtie2"), "-p", args.processors, "-q", os.path.join(args.output_path,"_tmp", "unmapped.fq"), "-x", os.path.join(args.path2genome,"genome"), "-S", os.path.join(args.output_path, "virus_al.sam")]
            ef._run_subprocesses(arg, "STATUS: Align reads bowtie2 global...", "aligning reads")
        elif args.aligner == "bbmap":
            arg=[os.path.join(bbmap2path, "bbmap.sh"), "ref="+bbmapgenomefile, "in="+os.path.join(args.output_path,"_tmp", "unmapped.fq"), "nodisk=t" ,"local=f",
                "covstats="+os.path.join(args.output_path, "constats.txt"), "covhist="+os.path.join(args.output_path,"covhist.txt"),
                "basecov="+os.path.join(args.output_path,"basecov.txt"), "bincov="+os.path.join(args.output_path, "bincov.txt"),
                "out="+os.path.join(args.output_path, "virus_al.sam")]
            ef._run_subprocesses(arg, "STATUS: Align reads bbmap global...", "aligning reads")
        elif args.aligner == "star":
            arg=[os.path.join(star2path, "STAR"), "--runThreadN", args.processors, "--genomeDir",  os.path.join(args.path2genome, "STAR_indicies"),
                "--outSAMtype", "SAM", "--readFilesIn", os.path.join(args.output_path,"_tmp", "unmapped.fq"), "--outFileNamePrefix", os.path.join(args.output_path, "virus_al")]
            ef._run_subprocesses(arg, "STATUS: Align reads star global..", "aligning reads star global")
    else:
        raise ValueError(args.alignment_type+" not an available set to global or local, all lowercase")

    # --generate necessary output files
    arg=[os.path.join(samtoolspath, "samtools"), "view", "-@", args.processors, "-F", "4", "-Sb", os.path.join(args.output_path,"virus_al.sam"), "-o", os.path.join(args.output_path, "virus_al.bam")]
    ef._run_subprocesses(arg, "STATUS: Extracting mapped reads and converting to bam", "extracting mapped reads and converting to bam")

    arg=[os.path.join(samtoolspath, "samtools"), "view", "-@", args.processors, "-h", os.path.join(args.output_path, "virus_al.bam"), "-o", os.path.join(args.output_path, "virus_al_mapped.sam")]
    ef._run_subprocesses(arg, "STATUS: Converting mapped reads to from bam to sam file", "Converting mapped reads to from bam to sam file")

    arg=[os.path.join(samtoolspath, "samtools"), "sort", os.path.join(args.output_path,"virus_al.bam"), "-o", os.path.join(args.output_path,"virus_al_sort.bam")]
    ef._run_subprocesses(arg, "STATUS: Sorting bam", "sorting bam")

    arg=[os.path.join(samtoolspath, "samtools"), "index", os.path.join(args.output_path,"virus_al_sort.bam"), os.path.join(args.output_path,"virus_al_sort.bam.bai")]
    ef._run_subprocesses(arg, "STATUS: generating bam index", "generating bam index")


