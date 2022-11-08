__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "integrate data in to matrix"

from scipy.sparse import coo_matrix
import scipy.io
import os
import shutil
import gzip
import csv
import numpy as np
import scqv.extra_functions as ef


def integrate_data_2_matrix(args, dfumi, gene_name):
    path10x = args.path10x

    if args.overwrite_feature_matrix:
        filter_folder="filtered_feature_bc_matrix"
        raw_folder="raw_feature_bc_matrix"
    else:
        filter_folder=args.output_filtered_folder+"_"+args.aligner
        raw_folder=args.output_raw_folder+"_"+args.aligner
        if os.path.isdir(os.path.join(path10x, "outs", filter_folder)) is False:
            shutil.copytree(os.path.join(path10x, "outs", args.input_filtered_folder), os.path.join(path10x, "outs", filter_folder), copy_function = shutil.copy)
        else:
            print("WARNING: filter folder exists so removing and recopying original")
            shutil.rmtree(os.path.join(path10x, "outs", filter_folder))
            shutil.copytree(os.path.join(path10x, "outs", args.input_filtered_folder), os.path.join(path10x, "outs", filter_folder), copy_function = shutil.copy)
        if os.path.isdir(os.path.join(path10x, "outs", raw_folder)) is False:
            shutil.copytree(os.path.join(path10x, "outs", args.input_raw_folder), os.path.join(path10x, "outs", raw_folder), copy_function = shutil.copy)
        else:
            print("WARNING: raw folder exists so removing and recopying original")
            shutil.rmtree(os.path.join(path10x, "outs", raw_folder))
            shutil.copytree(os.path.join(path10x, "outs", args.input_raw_folder), os.path.join(path10x, "outs", raw_folder), copy_function = shutil.copy)

    #####-- add data to filter folder files --#####

    # --Load in features and barcode matricies 
    print ("STATUS: Integrating viral copy counts into 10x matrix and feature files in filter folder")
    gene_names = [row for row in csv.reader(gzip.open(os.path.join(path10x,"outs", filter_folder, "features.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    barcodes = [row[0] for row in csv.reader(gzip.open(os.path.join(path10x,"outs", filter_folder, "barcodes.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]

    viralcopy = []

    # --Put counts in to a list 
    for i in dfumi.index.tolist():
        count = dfumi.loc[i,'umi']
        viralcopy.append([len(gene_names)+1, barcodes.index(i)+1, count])
    
    # -- Fill list with gene names (for new features.tsv file)
    gene_names.append((gene_name, gene_name, "Gene Expression"))

    # -- unzip matrix file so that it can be opened 
    arg=['gunzip', "-f", os.path.join(path10x,"outs",filter_folder,"matrix.mtx.gz")]
    ef._run_subprocesses(arg, "STATUS: unzipping matrix file", "unzipping matrix file")

    # -- Open matrix file 
    f_in = open(os.path.join(path10x,"outs",filter_folder,"matrix.mtx"), 'r')
    lines = f_in.readlines()

    # -- Edit header in the matrix file to have the correct new number of genes
    larray =lines[2].split(" ")
    larray[0]=len(gene_names)
    larray[2] = len(lines)+len(viralcopy)-3
    lines[2]=" ".join(str(v) for v in larray)
    lines[2]=lines[2]+"\n"

    # -- Load edited data into new matrix file 
    f_in = open(os.path.join(path10x,"outs",filter_folder,"matrix.mtx"), 'w')
    for line in lines:
        f_in.write(line)
    for i in viralcopy:
        f_in.write(str(i[0])+" "+str(i[1])+" "+str(i[2])+"\n")
    f_in.close()

    # -- re zip matrix file file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",filter_folder,"matrix.mtx")]
    ef._run_subprocesses(arg, "STATUS: zipping edited matrix file", "zipping edited matrix file")
    
    # -- load new information into new features file 
    f_in = open(os.path.join(path10x,"outs",filter_folder,"features.tsv"), 'w')
    for i in gene_names:
        f_in.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\n")
    f_in.close()

    # -- re zip features file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",filter_folder,"features.tsv")]
    ef._run_subprocesses(arg, "STATUS: zipping new features files filter folder", "zipping new features file filter folder")

    #####-- add data to raw folder files--#####
    print ("STATUS: Integrating viral copy counts into 10x matrix and feature files in raw folder")
    # --Load in features and barcode matricies 
    gene_names = [row for row in csv.reader(gzip.open(os.path.join(path10x,"outs", raw_folder, "features.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    barcodes = [row[0] for row in csv.reader(gzip.open(os.path.join(path10x,"outs", raw_folder, "barcodes.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    
    viralcopy = []

    # --Put counts in to a list 
    for i in dfumi.index.tolist():
        count = dfumi.loc[i,'umi']
        viralcopy.append([len(gene_names)+1, barcodes.index(i)+1, count])

    # -- Fill list with gene names (for new features.tsv file)
    gene_names.append((gene_name, gene_name, "Gene Expression"))

    # -- unzip matrix file so that it can be opened 
    arg=['gunzip', "-f", os.path.join(path10x,"outs",raw_folder,"matrix.mtx.gz")]
    ef._run_subprocesses(arg, "STATUS: unzipping matrix file", "unzipping matrix file")

    # -- Open matrix file 
    f_in = open(os.path.join(path10x,"outs",raw_folder,"matrix.mtx"), 'r')
    lines = f_in.readlines()

    # -- Edit header in the matrix file to have the correct new number of genes
    larray =lines[2].split(" ")
    larray[0]=len(gene_names)
    larray[2] = len(lines)+len(viralcopy)-3
    lines[2]=" ".join(str(v) for v in larray)
    lines[2]=lines[2]+"\n"
    # -- Load edited data into new matrix file 
    f_in = open(os.path.join(path10x,"outs",raw_folder,"matrix.mtx"), 'w')
    for line in lines:
        f_in.write(line)
    for i in viralcopy:
        f_in.write(str(i[0])+" "+str(i[1])+" "+str(i[2])+"\n")
    f_in.close()

    # -- re zip matrix file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",raw_folder,"matrix.mtx")]
    ef._run_subprocesses(arg, "STATUS: zipping edited matrix file", "zipping edited matrix file")

    # -- load new information into new features file 
    f_in = open(os.path.join(path10x,"outs",raw_folder,"features.tsv"), 'w')
    for i in gene_names:
        f_in.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\n")
    f_in.close()

    # -- re zip features file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",raw_folder,"features.tsv")]
    ef._run_subprocesses(arg, "STATUS: zipping new features files raw folder", "zipping new features file raw folder")

def integrate_viralgenes_data_2_matrix(args, dfumi):
    path10x = args.path10x

    if args.overwrite_feature_matrix:
        filter_folder="filtered_feature_bc_matrix"
        raw_folder="raw_feature_bc_matrix"
    else:
        filter_folder="filtered_feature_bc_matrix_"+args.aligner
        raw_folder="raw_feature_bc_matrix_"+args.aligner

        if os.path.isdir(os.path.join(path10x, "outs", filter_folder)) is False:
           shutil.copytree(os.path.join(path10x, "outs", "filtered_feature_bc_matrix"), os.path.join(path10x, "outs", filter_folder), copy_function = shutil.copy)
        if os.path.isdir(os.path.join(path10x, "outs", raw_folder)) is False:
           shutil.copytree(os.path.join(path10x, "outs", "raw_feature_bc_matrix"), os.path.join(path10x,  "outs", raw_folder), copy_function = shutil.copy)

    #####-- add data to filter folder files --#####
    print ("STATUS: Integrating viral gene counts into 10x matrix and feature files in filter folders")
    # --Load in features and barcode matricies 
    gene_names = [row for row in csv.reader(gzip.open(os.path.join(path10x,"outs",filter_folder,"features.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    barcodes = [row[0] for row in csv.reader(gzip.open(os.path.join(path10x,"outs",filter_folder,"barcodes.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]

    uniqgenes = set(dfumi['gene'].to_list())
    counts4genes = list()
    dfumi = dfumi.reset_index()
    counter=0

    # --Put counts in to a list 
    for gene in uniqgenes:
        dftemp =dfumi.loc[dfumi.gene==gene, ]  
        counter=counter+1
        for index, row in dftemp.iterrows(): 
            counts4genes.append([len(gene_names)+counter, barcodes.index(row["cell_barcode"])+1, row["umi"]])

    # -- Fill list with gene names (for new features.tsv file)
    for gene in uniqgenes:
        gene_names.append((gene, gene, "Gene Expression"))

    # -- unzip matrix file so that it can be opened 
    arg=['gunzip', "-f", os.path.join(path10x,"outs",filter_folder,"matrix.mtx.gz")]
    ef._run_subprocesses(arg, "STATUS: unzipping matrix file", "unzipping matrix file")

    # -- Open matrix file 
    f_in = open(os.path.join(path10x,"outs",filter_folder,"matrix.mtx"), 'r')
    lines = f_in.readlines()

    # -- Edit header in the matrix file to have the correct new number of genes
    larray =lines[2].split(" ")
    larray[0]=len(gene_names)
    larray[2] = len(lines)+len(counts4genes)-3
    lines[2]=" ".join(str(v) for v in larray)
    lines[2]=lines[2]+"\n"

    # -- Load edited data into new matrix file 
    f_in = open(os.path.join(path10x,"outs",filter_folder,"matrix.mtx"), 'w')
    for line in lines:
        f_in.write(line)
    for g in counts4genes:
        f_in.write(str(g[0])+" "+str(g[1])+" "+str(g[2])+"\n")
    f_in.close()

    # -- re zip matrix file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",filter_folder,"matrix.mtx")]
    ef._run_subprocesses(arg, "STATUS: zipping new matrix file", "zipping new matrix file")

    # -- load new information into new features file 
    f_in = open(os.path.join(path10x,"outs",filter_folder,"features.tsv"), 'w')
    for i in gene_names:
        f_in.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\n")
    f_in.close()

    # -- re zip features file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",filter_folder,"features.tsv")]
    ef._run_subprocesses(arg, "STATUS: zipping new features file in feature folder", "zipping new features file in feature folder")

    #####--add data to raw folder files--#####
    print ("STATUS: Integrating viral gene counts into 10x matrix and feature files in raw folder")
    
    # --Load in features and barcode matricies 
    gene_names = [row for row in csv.reader(gzip.open(os.path.join(path10x,"outs",raw_folder,"features.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    barcodes = [row[0] for row in csv.reader(gzip.open(os.path.join(path10x,"outs",raw_folder,"barcodes.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    
    uniqgenes = set(dfumi['gene'].to_list())
    counts4genes = list()
    dfumi = dfumi.reset_index()
    counter=0

    # --Put counts in to a list 
    for gene in uniqgenes:
        dftemp =dfumi.loc[dfumi.gene==gene, ]  
        counter=counter+1
        for index, row in dftemp.iterrows(): 
            counts4genes.append([len(gene_names)+counter, barcodes.index(row["cell_barcode"])+1, row["umi"]])

    # -- Fill list with gene names (for new features.tsv file)
    for gene in uniqgenes:
        gene_names.append((gene, gene, "Gene Expression"))

    # -- unzip matrix file so that it can be opened 
    arg=['gunzip', "-f", os.path.join(path10x,"outs",raw_folder,"matrix.mtx.gz")]
    ef._run_subprocesses(arg, "STATUS: unzipping matrix file", "unzipping matrix file")

    # -- Open matrix file 
    f_in = open(os.path.join(path10x,"outs",raw_folder,"matrix.mtx"), 'r')
    lines = f_in.readlines()

    # -- Edit header in the matrix file to have the correct new number of genes
    larray =lines[2].split(" ")
    larray[0]=len(gene_names)
    larray[2] = len(lines)+len(counts4genes)-3
    lines[2]=" ".join(str(v) for v in larray)
    lines[2]=lines[2]+"\n"

    # -- Load edited data into new matrix file 
    f_in = open(os.path.join(path10x,"outs",raw_folder,"matrix.mtx"), 'w')
    for line in lines:
        f_in.write(line)
    for g in counts4genes:
        f_in.write(str(g[0])+" "+str(g[1])+" "+str(g[2])+"\n")
    f_in.close()

    # -- re zip matrix file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",raw_folder,"matrix.mtx")]
    ef._run_subprocesses(arg, "STATUS: zipping new matrix file", "zipping new matrix file")

    # -- load new information into new features file 
    f_in = open(os.path.join(path10x,"outs",raw_folder,"features.tsv"), 'w')
    for i in gene_names:
        f_in.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\n")
    f_in.close()

    # -- re zip features file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",raw_folder,"features.tsv")]
    ef._run_subprocesses(arg, "STATUS: zipping new features file in raw folder", "zipping new features file in raw folder")    