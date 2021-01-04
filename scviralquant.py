__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "Quantifies single cell viral reads"

# ----Libraries
import os
import argparse
from sys import platform
from scqv.map_reads import process_reads as pr
from scqv.map_reads import map_reads as mr
from scqv.process_10x import extract10x as ex
from scqv.quantify import viral_copies as vc
from scqv.quantify import viral_genes as vg
from scqv.quantify import integrate
from scqv.visualization import viz 
PATH = os.path.dirname(os.path.abspath(__file__))
if platform == "linux":
    bowtie2path = PATH+"/aligntools/bowtie2-2.4.2-linux-x86_64/"
    samtoolspath = PATH+"/extra_tools/samtoolsv1.11_linux/bin/"
elif platform == "OS":
    bowtie2path = PATH+"/aligntools/bowtie2-2.4.2-macos-x86_64/"
else:
    ValueError("Program wont run on this operating system "+platform)

def parse_arguments():
    parser = argparse.ArgumentParser(prog="scViralQuant", description="Quantify sc viral mapping reads")
    parser.add_argument("-10x", "--path10x", type=str, required=True)
    parser.add_argument("-op", "--output_path", default="output_scviralquant/")
    parser.add_argument("-p", "--processors", default=1)
    parser.add_argument("-p2genome", "--path2genome", type=str, required=True)

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    viable_cb=ex.extract_viable_10x(args.path10x)
    pr.process_unmapped_reads(args, samtoolspath, viable_cb)
    mr.map2viralgenome(args, bowtie2path, samtoolspath)
    dfvc, gene_name = vc.htseq_run(args)
    viz.generate_viral_copy_plots(args, gene_name, dfvc)
    integrate.integrate_data_2_matrix(args.path10x, dfvc, gene_name)
    dfgenes = vg.htseq_run(args)
    if dfgenes is not False:
        viz.generate_viral_gene_plots(args, dfgenes)
        integrate.integrate_viralgenes_data_2_matrix(args.path10x, dfgenes)
