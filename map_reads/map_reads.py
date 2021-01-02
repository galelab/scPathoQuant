__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "map reads with bowtie (currently only version)"

import os
import glob
import subprocess
import extra_functions as ef


def map2viralgenome(args, bowtie2path, samtoolspath):

    # -- check if libraries are present if not generate 
    files = glob.glob(args.path2genome+"/*.bt2")
    if len(files) == 0:
        files_genome = glob.glob(args.path2genome+".fa")
        if len(files_genome) > 1:
            print ("WARNING:  two fasta files in genome folder.")
        arg=[bowtie2path+"bowtie2-build", files_genome[0], args.path2genome+"genome"]
        ef._run_subprocesses(arg, "STATUS: generating libraries ...", "extracting unmapped")

    else:
        print("STATUS: libraries have already been made for this genome")

    # -- align reads 
    arg=[bowtie2path+"bowtie2", "-p", args.processors, "-q", os.path.join(args.output_path,"_tmp/unmmapped.fq"), "-x", os.path.join(args.path2genome,"genome"), "-S", os.path.join(args.output_path, "virus_aligned.sam")]
    ef._run_subprocesses(arg, "STATUS: Align reads ...", "aligning reads", verbose=True)
    
    # --generate necessary output files
    arg=[samtoolspath+"samtools", "view", "-@", args.processors, "-F", "4", "-Sb", os.path.join(args.output_path,"virus_aligned.sam"), "-o", os.path.join(args.output_path, "virus_aligned.bam")]
    ef._run_subprocesses(arg, "STATUS: Extracting mapped reads and converting to bam", "extracting mapped reads and converting to bam")

    arg=[samtoolspath+"samtools", "sort", os.path.join(args.output_path,"virus_aligned.bam"), "-o", os.path.join(args.output_path,"viral_aligned_sort.bam")]
    ef._run_subprocesses(arg, "STATUS: Sorting bam", "sorting bam")

    arg=[samtoolspath+"samtools", "index", os.path.join(args.output_path,"virus_aligned_sort.bam"), os.path.join(args.output_path,"viral_aligned_sort.bam.bai")]
    ef._run_subprocesses(arg, "STATUS: generating bam index", "generating bam index")


