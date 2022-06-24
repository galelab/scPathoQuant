__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "quantify viral genes"

import os
import re
import glob
import pandas as pd
import scqv.extra_functions as ef

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
                    reads_mapping[larray[0]]=re.sub("XF:\w+:", "", larray[-1])

    ## -- get umi and cell barcode information 
    if len(reads_mapping.keys()) > 0:
        tmpdf = pd.DataFrame.from_dict(reads_mapping,orient='index')   
        df = pd.read_csv(os.path.join(output_path,"_tmp", "barcode_umi_read_table.csv"))
        df1 = df[df['read'].isin(reads_mapping.keys())]
        df1 = df1.set_index('read')
        df1 = pd.concat([df1, tmpdf],  axis=1, join="inner")
        df1 = df1.reset_index() 
        # print (df1)
        df1 = df1.rename(columns={"index":"read", 0: "gene"})
        # print (df1)
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
        if df2.empty:
            print("STATUS: no reads able to be quantifie to viral genes.. will not add genes to 10x filtered files")
            return(False)
        else:
            df_umi = df2[["cell_barcode", "umi", "gene"]].groupby(["cell_barcode", "gene"]).count()
            df_umi = df_umi.reset_index()
            df_umi = df_umi.set_index('cell_barcode')
            df_umi.to_csv(os.path.join(output_path, "virus_al_gene_counts.csv"))
            return(df_umi)
            
    else: 
        print("STATUS: no reads mapping to viral genes.. will not add genes to 10x filtered files")
        df_umi = pd.DataFrame(columns=["cell_barcode", "gene", "umi"])
        df_umi = df_umi.set_index('cell_barcode')
        df_umi.to_csv(os.path.join(output_path, "virus_al_gene_counts.csv"))
        return(False)

def htseq_run(args):
    ## -- generate gtf for viral genes
    files_gtf = glob.glob(os.path.join(args.path2genome, "*.gtf"))
    if len(files_gtf) > 1:
        print ("WARNING: too many gtf files in genome folder.")
        return False
    elif len(files_gtf) == 0:
        print ("STATUS: No gtf file not counting reads for individual viral genes")
        return False
    elif len(files_gtf) ==1:
        arg=["htseq-count", "--format=bam", "--mode=intersection-nonempty", "--idattr=gene_id", os.path.join(args.output_path, "virus_al_sort.bam"),
            files_gtf[0], "--samout="+os.path.join(args.output_path,"virus_genes_al_sort_counts.sam")]
        ef._run_subprocesses(arg, "STATUS: running htseq for viral genes ", "running htseq for viral genes")
        dfumi = quantify_reads(args.output_path, os.path.join(args.output_path, "virus_genes_al_sort_counts.sam"))
        return dfumi

