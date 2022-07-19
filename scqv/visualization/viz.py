__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "visualize counts"

import os
import pandas as pd
import seaborn as sns

def generate_viral_copy_plots(args, gene_name, dfumi):
    print ("STATUS: Generating violin plots for viral copies")
    if dfumi.shape[0] > 0:
        sns.set_theme(style="whitegrid")
        dfumi[gene_name]=[gene_name]*dfumi.shape[0]
        ax = sns.violinplot(x=gene_name, y="umi",data=dfumi)
        ax.set(xlabel='', ylabel="UMI counts")
        g = ax.get_figure()
        g.subplots_adjust(left=0.2)
        g.set_figheight(4)
        g.set_figwidth(4)
        g.savefig(os.path.join(args.output_path, "viral_copy.png"), dpi=500)
        g.savefig(os.path.join(args.output_path, "viral_copy.svg"), dpi=500)

    else:
        print ("STATUS: No cells have viral reads not generating violin plot")

def generate_viral_gene_plots(args, dfumi):
    print ("STATUS: Generating violin plots for viral copies")
    if dfumi.shape[0] > 0:
        sns.set_theme(style="whitegrid")
        # dfumi[gene_name]=[gene_name]*dfumi.shape[0]
        ax = sns.violinplot(x="gene", y="umi",data=dfumi)
        xlabels = ax.get_xticklabels()
        ax.set_xticklabels(labels=xlabels, rotation=90)
        ax.set(xlabel='genes', ylabel="UMI counts")
        g = ax.get_figure()
        g.subplots_adjust(left=0.2,bottom=0.5)
        g.set_figheight(4)
        g.set_figwidth(6)
        g.savefig(os.path.join(args.output_path, "viral_genes.png"), dpi=500)
        g.savefig(os.path.join(args.output_path, "viral_genes.svg"), dpi=500)

    else:
        print ("STATUS: No cells have viral reads not generating violin plot for viral genes")
