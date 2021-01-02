__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "visualize counts"

import pandas as pd
import seaborn as sns

def generate_viral_copy_plots(args, gene_name, dfumi):
    print ("STATUS: Generating violin plots for viral copies")
    sns.set_theme(style="whitegrid")
    dfumi[gene_name]=["HIV"]*dfumi.shape[0]
    ax = sns.violinplot(x="HIV", y="umi",data=dfumi)
    g = ax.get_figure()
    g.savefig(args.output_path+"viral_copy.png", dpi=500)

