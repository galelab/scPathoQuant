__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "map reads with bowtie/bbmap"

import re
import os
import glob
import subprocess
import spq.extra_functions as ef
import pysam

PATH = os.path.dirname(os.path.abspath(__file__))
PATH = re.sub("map_reads", "", PATH)

def _get_genome_file(args):
    files_genome = glob.glob(os.path.join(args.path2genome, "*.fa"))
    if (len(files_genome)== 0):
        files_genome = glob.glob(os.path.join(args.path2genome, "*.fasta"))
        if (len(files_genome) == 0):
            files_genome = glob.glob(os.path.join(args.path2genome, "*.fna"))
            if (len(files_genome) == 0):
                AssertionError("NO GENOME FILE WAS FOUND, NEEDS TO END IN .fa, .fasta, or .fna")
    return(files_genome)

def map2pathogengenome(args):
    # -- check if libraries are present if not generate 
    if args.aligner == "bowtie2":
        files = glob.glob(os.path.join(args.path2genome,"*.bt2"))
        if len(files) == 0:
            files_genome = _get_genome_file(args)
            bowtie2genomefile = files_genome[0]
            if len(files_genome) > 1:
                print("WARNING:  two fasta files in genome folder using this one "+files_genome[0])
            arg=["bowtie2-build", bowtie2genomefile, os.path.join(args.path2genome, "genome")]
            ef._run_subprocesses(arg, "STATUS: generating libraries bowtie2 ...", "extracting unmapped")

        else:
            print("STATUS: bowtie2 indexes have already been made for this genome")
    elif args.aligner == "bbmap":
        files_genome = _get_genome_file(args)
        bbmapgenomefile = files_genome[0]
        if len(files_genome) > 1:
            print("WARNING:  two fasta files in genome folder using this one "+files_genome[0])
    else:
        raise ValueError(args.aligner+" not an available aligner specify bbmap or bowtie2, all lowercase")

    # -- align reads 
    if args.aligner == "bowtie2":
        if args.bowtie2_params is None:
            arg=["bowtie2", "-p", args.processors,
            "-q", os.path.join(args.output_path,"_tmp", "unmapped.fq.gz"), 
            "-x", os.path.join(args.path2genome,"genome"), 
            "-S", os.path.join(args.output_path, "pathogen_al.sam")]
        else:
            arg=["bowtie2", "-p", args.processors,
            "-q", os.path.join(args.output_path,"_tmp", "unmapped.fq.gz"), 
            "-x", os.path.join(args.path2genome,"genome"), 
            "-S", os.path.join(args.output_path, "pathogen_al.sam")] + args.bowtie2_params
            
        ef._run_subprocesses(arg, "STATUS: Align reads bowtie2 global...", "aligning reads")
    elif args.aligner == "bbmap":
        if args.bbmap_params is None:
            arg=["bbmap.sh", "ref="+bbmapgenomefile, "in="+os.path.join(args.output_path,"_tmp", "unmapped.fq.gz"), "nodisk=t" ,"local=f",
                "covstats="+os.path.join(args.output_path, "constats.txt"), "covhist="+os.path.join(args.output_path,"covhist.txt"), "minid=0.95",
                "basecov="+os.path.join(args.output_path,"basecov.txt"), "bincov="+os.path.join(args.output_path, "bincov.txt"),
                "out="+os.path.join(args.output_path, "pathogen_al.sam")]
        else:
            arg=["bbmap.sh", "ref="+bbmapgenomefile, "in="+os.path.join(args.output_path,"_tmp", "unmapped.fq.gz"), "nodisk=t" ,"local=f",
                "covstats="+os.path.join(args.output_path, "constats.txt"), "covhist="+os.path.join(args.output_path,"covhist.txt"), "minid=0.95",
                "basecov="+os.path.join(args.output_path,"basecov.txt"), "bincov="+os.path.join(args.output_path, "bincov.txt"),
                "out="+os.path.join(args.output_path, "pathogen_al.sam")] + args.bbmap_params
        ef._run_subprocesses(arg, "STATUS: Align reads bbmap global...", "aligning reads")

    # --generate necessary output files
    print("STATUS: Extracting mapped reads and converting to bam")
    pysam.view("-@", str(args.processors),"-S", "-b", os.path.join(args.output_path,"pathogen_al.sam"), 
               "-o", os.path.join(args.output_path, "pathogen_al.bam"), catch_stdout=False)

    print("STATUS: Converting mapped reads to from bam to sam file")
    pysam.view("-@", str(args.processors),  "-b", "-F", "4", os.path.join(args.output_path, "pathogen_al.bam"), 
                "-o",  os.path.join(args.output_path, "pathogen_al_mapped.sam"), catch_stdout=False)

    print("STATUS: Sorting bam")
    pysam.sort("-@", str(args.processors), "-o", os.path.join(args.output_path,"pathogen_al_mapped_sort.bam"), 
            os.path.join(args.output_path,"pathogen_al_mapped.sam"))

    print( "STATUS: generating bam index")
    pysam.index( os.path.join(args.output_path,"pathogen_al_mapped_sort.bam"))
    
    print("STATUS: generating fastq file for all mapped reads")
    pysam.bam2fq('-n', '-0', os.path.join(args.output_path, "pathogen_al_mapped.fq.gz"), os.path.join(args.output_path,"pathogen_al_mapped_sort.bam"),  catch_stdout=False)
