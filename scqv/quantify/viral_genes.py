__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "quantify viral genes"

import os
import re
import glob
import pandas as pd
import extra_functions as ef

def quantify_reads(output_path,filename):
    
    print ("STATUS: quantifying viral copies ...")
    ## -- get reads mapping to virus
    reads_mapping = dict()
    with open(filename, 'r') as fin:
        for line in fin:
            larray = line.strip().split("\t")
            if larray[-1] != "XF:Z:__no_feature" and larray[-1] !="XF:Z:__too_low_aQual":
                if larray[-1].startswith("XF:Z:__ambiguous"):
                    pass
                else:
                    reads_mapping[larray[0]]=re.sub("XF:Z:", "", larray[-1])

    ## -- get umi and cell barcode information 
    tmpdf = pd.DataFrame.from_dict(reads_mapping,orient='index')   
    df = pd.read_csv(os.path.join(output_path,"_tmp", "barcode_umi_read_table.csv"))
    df1 = df[df['read'].isin(reads_mapping.keys())]
    df1 = df1.set_index('read')
    df1 = pd.concat([df1, tmpdf],  axis=1, join="inner")
    df1 = df1.reset_index() 
    df1 = df1.rename(columns={"index":"read", 0: "gene"})
    df2 = df1.groupby(['cell_barcode', 'umi', 'gene']).agg({'gene': lambda x: ','.join(x)})
    df2.columns = ['combined_genes']
    df2 = df2.reset_index()
    drop_reads = []
    for i in range(0, len(df2.combined_genes)): 
        x = len(set(df2.combined_genes[i].split(','))) 
        if x > 1:
            drop_reads.append(i)
    if len(drop_reads) > 0:
        print ("STATUS: Dropping "+str(len(drop_reads))+" UMIs because of inconsistent gene annotations")
        df2 = df2.drop(drop_reads)

    df_umi = df2[["cell_barcode", "umi", "gene"]].groupby(["cell_barcode", "gene"]).count()
    df_umi = df_umi.reset_index()
    df_umi = df_umi.set_index('cell_barcode')
    df_umi.to_csv(os.path.join(output_path, "virus_al_gene_counts.csv"))
    return(df_umi)

def htseq_run(args):
    ## -- generate gtf for viral genes
    files_gtf = glob.glob(args.path2genome+"*.gtf")
    if len(files_gtf) > 1:
        print ("WARNING: too many gtf files in genome folder.")
        return False
    elif len(files_gtf) == 0:
        print ("STATUS: No gtf file not counting reads for individual viral genes")
        return False
    elif len(files_gtf) ==1:
        arg=["htseq-count", "--format=bam", "--idattr=gene_id", os.path.join(args.output_path, "virus_al_sort.bam"),
            files_gtf[0], "--samout="+os.path.join(args.output_path,"virus_genes_al_sort_gene_counts.sam")]
        ef._run_subprocesses(arg, "STATUS: running htseq for viral genes ", "running htseq for viral genes")
        dfumi = quantify_reads(args.output_path, os.path.join(args.output_path, "virus_al_sort_gene_counts.sam"))
        return dfumi

