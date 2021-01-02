__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "Quantifies single cell viral reads"

# ----Libraries
import os
import argparse
from sys import platform
from map_reads import process_reads as pr
from map_reads import map_reads as mr
from process_10x import extract10x as ex
from quantify import viral_copies as vc
from visualization import viz 
PATH = os.path.dirname(os.path.abspath(__file__))
samtoolspath = PATH+"/extra_tools/samtoolsv1.11_linux/bin/"
if platform == "linux":
    bowtie2path = PATH+"/aligntools/bowtie2-2.4.2-linux-x86_64/"
elif platform == "OS":
    bowtie2path = PATH+"/aligntools/bowtie2-2.4.2-macos-x86_64/"
else:
    ValueError("Program wont run on this operating system "+platform)

args = parse_arguments()

def parse_arguments():
    parser = argparse.ArgumentParser(prog="scViralQuant", description="Quantify sc viral mapping reads")
    parser.add_argument("-10x", "--path10x", type=str, required=True)
    parser.add_argument("-path2bam", "--path2possorted_genome_bam", required=True)
    parser.add_argument("-op", "--output_path", default="output_scviralquant/")
    parser.add_argument("-p", "--processors", default=1)
    parser.add_argument("-p2genome", "--path2genome", type=str, required=True)

    return parser.parse_args()

args = parse_arguments()
viable_cb=ex.extract_viable_10x(args.path10x)
pr.process_unmapped_reads(args, samtoolspath, viable_cb)
mr.map2viralgenome(args, bowtie2path, samtoolspath)
dfvc, gene_name = vc.htseq_run(args)
viz.generate_viral_copy_plots(args, dfvc, gene_name)
