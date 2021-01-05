__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "Quantifies single cell viral reads"

# ----Libraries
import os
import argparse
import shutil
from scqv.map_reads import process_reads as pr
from scqv.map_reads import map_reads as mr
from scqv.process_10x import extract10x as ex
from scqv.quantify import viral_copies as vc
from scqv.quantify import viral_genes as vg
from scqv.quantify import integrate
from scqv.visualization import viz 


def parse_arguments():
    parser = argparse.ArgumentParser(prog="scViralQuant", description="Quantify sc viral mapping reads")
    parser.add_argument("-10x", "--path10x", type=str, required=True)
    parser.add_argument("-op", "--output_path", default="output_scviralquant/")
    parser.add_argument("-p", "--processors", default=1)
    parser.add_argument("-p2genome", "--path2genome", type=str, required=True)
    parser.add_argument("--tmp_removal", type=bool, required=False, default=True)

    return parser.parse_args()
if __name__ == '__main__':
    args = parse_arguments()
    viable_cb = ex.extract_viable_10x(args.path10x)
    pr.process_unmapped_reads(args, viable_cb)
    mr.map2viralgenome(args)
    dfvc, gene_name = vc.htseq_run(args)
    viz.generate_viral_copy_plots(args, gene_name, dfvc)
    integrate.integrate_data_2_matrix(args.path10x, dfvc, gene_name)
    dfgenes = vg.htseq_run(args)
    if dfgenes is not False:
        viz.generate_viral_gene_plots(args, dfgenes)
        integrate.integrate_viralgenes_data_2_matrix(args.path10x, dfgenes)
    if args.tmp_removal:
        shutil.rmtree(os.path.join(args.output_path,"_tmp"))