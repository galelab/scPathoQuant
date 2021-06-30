__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "map reads with bowtie (currently only version)"

import re
import os
import glob
import subprocess
import scqv.extra_functions as ef
from sys import platform
import pandas as pd
from tqdm import tqdm
from scqv.variantcalling import sam2tsv as st

PATH = os.path.dirname(os.path.abspath(__file__))
PATH = re.sub("variantcalling", "", PATH)
if platform == "linux":
    samtoolspath = os.path.join(PATH, "extra_tools", "samtoolsv1.11_linux","bin")
    bcftoolspath =  os.path.join(PATH, "extra_tools", "bcftools1.10.2_linux", "bin")
    vcfutilspath =  os.path.join(PATH, "extra_tools", "vcfutils", "bin")


def variantcaller(args, dfgenes):
    '''Runs samtools/bcftools to call variants in the virus genome'''
    files_genome = glob.glob(os.path.join(args.path2genome, "*.fa"))
    genomefile = files_genome[0]
    header = open(genomefile, 'r').readline()
    header = re.sub(">", "", header).strip()

    if dfgenes is not False:
        gtf_coord = {}
        files_gtf = glob.glob(os.path.join(args.path2genome, "*.gtf"))
        with open(files_gtf[0]) as fin:
            for line in fin:
                if line.startswith("#"):
                    pass
                else:
                    larray = line.strip().split("\t")
                    if larray[2] == "gene":
                        tmp = larray[8].split(";")
                        for i in tmp:
                            if i.startswith("gene_name"):
                                gene = re.sub("gene_name ", "", i)
                                gene = re.sub("\"", "", gene)
                    for i in range(larray[3], larray[4]+1):
                        gtf_coord[i]=gene

    # arg=[os.path.join(samtoolspath, "samtools"), "faidx", genomefile, "-o", genomefile+".fai"]
    # ef._run_subprocesses(arg, "STATUS: indexing genome file", "indexing genome file")

    # arg=[os.path.join(bcftoolspath, "bcftools"), "mpileup", "--threads", args.processors, "-Q", "20", "-d", "250000000", "-L", "250000000", "-f", genomefile, os.path.join(args.output_path,"virus_al_sort.bam"), "-o", os.path.join(args.output_path,"virus.bcf")]
    # ef._run_subprocesses(arg, "STATUS: Running bcftools mpileup", "running bcftools mpileup")

    # arg=[os.path.join(bcftoolspath, "bcftools"), "call", "-v", "-c", os.path.join(args.output_path,"virus.bcf"), "-o", os.path.join(args.output_path,"virus_calls.vcf")]
    # ef._run_subprocesses(arg, "STATUS: Calling snps", "calling snps")

    variants = {}
    variants_ref = {}
    with open(os.path.join(args.output_path,"virus_calls.vcf")) as fin:
        for line in fin:
            if line.startswith("#"):
                pass
            else:
                larray=line.strip().split("\t")
                varis=larray[4].split(",")
                for v in varis:
                    variants.setdefault(larray[1], []).append(v)
                variants_ref.setdefault(larray[1], larray[3])
    print("STATUS: Number of variants found "+str(len(variants.keys())))
    print("STATUS: Loading barcode cell information ...")
    df = pd.read_csv(os.path.join(args.output_path,"_tmp", "barcode_umi_read_table.csv"))
    reads = df.read.to_list()
    dfvirus = pd.read_csv(os.path.join(args.output_path, "virus_al_counts.csv"),index_col=0)
    maxdepth = dfvirus.max()[0]

    # -- make new output folder 
    try:
        os.mkdir(os.path.join(args.output_path, "_tmp", "_variant_tmp"))
    except:
        # print("WARNING: variant tmp folder already exists")
        pass

    try:
        os.mkdir(os.path.join(args.output_path, "variant_calling_results"))
    except:
        pass 

    print("STATUS: Extract reads with variant  ...")
    outputfile = open(os.path.join(args.output_path, "variant_calling_results", "cellvariantcallstotal.txt"), 'w')
    outputfile.write("cell_barcode\tvariantinfo(variant_ref_base)\tpercentageofumi\tweight(virus umi count in cell/max virus umi count)\tweightedpercentage\tbase\n")
    for base in tqdm(variants):
        st.sam2tsv_function(os.path.join(args.output_path, "virus_al_sort.bam"), genomefile, int(base), outputfile=os.path.join(args.output_path, "_tmp", "_variant_tmp", "variant_"+str(base)+"_reads.txt"))
        # print("STATUS: Pull out variant information for "+base)
        outputfile1 = open(os.path.join(args.output_path, "variant_calling_results",  "cellvariantcalls_"+str(base)+".txt"), 'w')
        outputfile1.write("cell_barcode\tvariantinfo(variant_ref_base)\tpercentageofumi\tweight(virus umi count in cell/max virus umi count)\tweightedpercentage\tbase\n")
        dftmp = pd.read_table(os.path.join(args.output_path, "_tmp", "_variant_tmp", "variant_"+str(base)+"_reads.txt"), header=None)
        reads_tmp = set(reads).intersection(set(dftmp[0].to_list()))
        dftmpreads = dftmp[dftmp[0].isin(reads_tmp)]
        dftmpbarcode = df[df['read'].isin(reads_tmp)]
        dftmpreads = dftmpreads.set_index(0)
        dftmpbarcode = dftmpbarcode.set_index('read')
        dfc = pd.concat([dftmpreads, dftmpbarcode], axis=1, join="inner")
        # dfc.to_csv(os.path.join(args.output_path, "_tmp", "_variant_tmp", "variant_combined_df_"+str(base)+"_reads.txt"))
        dfc = dfc[['cell_barcode', 'umi', 4]]
        dfcbase = dfc.groupby(['cell_barcode',4]).count()
        dfcsum = dfcbase.groupby(['cell_barcode']).sum()
        for i in dfcsum.index.to_list():
            for vari in variants[str(base)]:
                try:
                    per = round(dfcbase.loc[i, vari]['umi']/dfcsum.loc[i]['umi'], 2)*100
                    depthweight = dfvirus.loc[i, "umi"]/maxdepth
                    weightedper = per*depthweight
                    if dfgenes is False:
                        outputfile.write(i+"\t"+vari+"_"+variants_ref[str(base)]+"_"+str(base)+"\t"+str(per)+"\t"+str(depthweight)+"\t"+str(weightedper)+"\t"+str(base)+"\n")
                        outputfile1.write(i+"\t"+vari+"_"+variants_ref[str(base)]+"_"+str(base)+"\t"+str(per)+"\t"+str(depthweight)+"\t"+str(weightedper)+"\t"+str(base)+"\n")
                    else:
                        try:
                            gene_variant=gtf_coord[str(base)]
                            outputfile.write(i+"\t"+vari+"_"+variants_ref[str(base)]+"_"+str(base)+"_"+gene_variant+"\t"+str(per)+"\t"+str(depthweight)+"\t"+str(weightedper)+"\t"+str(base)+"\n")
                            outputfile1.write(i+"\t"+vari+"_"+variants_ref[str(base)]+"_"+str(base)+"_"+gene_variant+"\t"+str(per)+"\t"+str(depthweight)+"\t"+str(weightedper)+"\t"+str(base)+"\n")
                        except KeyError:
                            gene_variant="OCR"
                            outputfile.write(i+"\t"+vari+"_"+variants_ref[str(base)]+"_"+str(base)+"_"+gene_variant+"\t"+str(per)+"\t"+str(depthweight)+"\t"+str(weightedper)+"\t"+str(base)+"\n")
                            outputfile1.write(i+"\t"+vari+"_"+variants_ref[str(base)]+"_"+str(base)+"_"+gene_variant+"\t"+str(per)+"\t"+str(depthweight)+"\t"+str(weightedper)+"\t"+str(base)+"\n")
                except KeyError:
                    pass
                    # print("WARNING: No variant "+vari+" in cell "+i)
        outputfile1.close()
    outputfile.close()