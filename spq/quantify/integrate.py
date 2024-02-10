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
import spq.extra_functions as ef


def integrate_data_2_matrix(args, dfumidict, virus_names):
    path10x = args.path10x

    if args.overwrite_feature_matrix:
        filter_folder="filtered_feature_bc_matrix"
        raw_folder="raw_feature_bc_matrix"
    else:
        filter_folder=args.output_filtered_folder+"_scPathoQuant_"+args.aligner
        raw_folder=args.output_raw_folder+"_scPathoQuant_"+args.aligner
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
    print ("STATUS: Integrating pathogen copy counts into 10x matrix and feature files in filter folder")
    feature_names_filtered = [row for row in csv.reader(gzip.open(os.path.join(path10x,"outs", filter_folder, "features.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    barcodes_filtered = [row[0] for row in csv.reader(gzip.open(os.path.join(path10x,"outs", filter_folder, "barcodes.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    # --Load in features and barcode matricies 
    feature_names_raw = [row for row in csv.reader(gzip.open(os.path.join(path10x,"outs", raw_folder, "features.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    barcodes_raw = [row[0] for row in csv.reader(gzip.open(os.path.join(path10x,"outs", raw_folder, "barcodes.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    
    viralcopy_filtered = []
    viralcopy_raw = []
    # --Put counts in to a list 
    viruswithcounts=[]
    for virus_name, dfumi in dfumidict.items():
        if dfumi.shape[0] > 0:
            viruswithcounts.append(virus_name)
            for i in dfumi.index.tolist():
                count = dfumi.loc[i,'umi']
                viralcopy_filtered.append([len(feature_names_filtered)+1, barcodes_filtered.index(i)+1, count])
                viralcopy_raw.append([len(feature_names_raw)+1, barcodes_raw.index(i)+1, count])
        else: 
            print ("STATUS: Not adding counts for pathogen "+virus_name+" to matrix or feature files because there are no counts")

    # -- Fill list with gene names (for new features.tsv file)
    for virus_name in viruswithcounts:
        feature_names_filtered.append((virus_name, virus_name, "Gene Expression"))
        feature_names_raw.append((virus_name, virus_name, "Gene Expression"))

    # -- unzip matrix file so that it can be opened 
    arg=['gunzip', "-f", os.path.join(path10x,"outs",filter_folder,"matrix.mtx.gz")]
    ef._run_subprocesses(arg, "STATUS: unzipping matrix file", "unzipping matrix file")

    # -- Open matrix file 
    f_in = open(os.path.join(path10x,"outs",filter_folder,"matrix.mtx"), 'r')
    lines = f_in.readlines()

    # -- Edit header in the matrix file to have the correct new number of genes
    larray =lines[2].split(" ")
    larray[0]=len(feature_names_filtered)
    larray[2] = len(lines)+len(viralcopy_filtered)-3
    lines[2]=" ".join(str(v) for v in larray)
    lines[2]=lines[2]+"\n"

    # -- Load edited data into new matrix file 
    f_in = open(os.path.join(path10x,"outs",filter_folder,"matrix.mtx"), 'w')
    for line in lines:
        f_in.write(line)
    for i in viralcopy_filtered:
        f_in.write(str(i[0])+" "+str(i[1])+" "+str(i[2])+"\n")
    f_in.close()

    # -- re zip matrix file file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",filter_folder,"matrix.mtx")]
    ef._run_subprocesses(arg, "STATUS: zipping edited matrix file", "zipping edited matrix file")
    
    # -- load new information into new features file 
    f_in = open(os.path.join(path10x,"outs",filter_folder,"features.tsv"), 'w')
    for i in feature_names_filtered:
        f_in.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\n")
    f_in.close()

    # -- re zip features file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",filter_folder,"features.tsv")]
    ef._run_subprocesses(arg, "STATUS: zipping new features files filter folder", "zipping new features file filter folder")

    #####-- add data to raw folder files--#####
    print ("STATUS: Integrating pathogen copy counts into 10x matrix and feature files in raw folder")

    # -- unzip matrix file so that it can be opened 
    arg=['gunzip', "-f", os.path.join(path10x,"outs",raw_folder,"matrix.mtx.gz")]
    ef._run_subprocesses(arg, "STATUS: unzipping matrix file", "unzipping matrix file")

    # -- Open matrix file 
    f_in = open(os.path.join(path10x,"outs",raw_folder,"matrix.mtx"), 'r')
    lines = f_in.readlines()

    # -- Edit header in the matrix file to have the correct new number of genes
    larray =lines[2].split(" ")
    larray[0]=len(feature_names_raw)
    larray[2] = len(lines)+len(viralcopy_raw)-3
    lines[2]=" ".join(str(v) for v in larray)
    lines[2]=lines[2]+"\n"
    # -- Load edited data into new matrix file 
    f_in = open(os.path.join(path10x,"outs",raw_folder,"matrix.mtx"), 'w')
    for line in lines:
        f_in.write(line)
    for i in viralcopy_raw:
        f_in.write(str(i[0])+" "+str(i[1])+" "+str(i[2])+"\n")
    f_in.close()

    # -- re zip matrix file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",raw_folder,"matrix.mtx")]
    ef._run_subprocesses(arg, "STATUS: zipping edited matrix file", "zipping edited matrix file")

    # -- load new information into new features file 
    f_in = open(os.path.join(path10x,"outs",raw_folder,"features.tsv"), 'w')
    for i in feature_names_raw:
        f_in.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\n")
    f_in.close()

    # -- re zip features file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",raw_folder,"features.tsv")]
    ef._run_subprocesses(arg, "STATUS: zipping new features files raw folder", "zipping new features file raw folder")

def integrate_pathogenes_data_2_matrix(args, dfumidict):
    path10x = args.path10x

    if args.overwrite_feature_matrix:
        filter_folder="filtered_feature_bc_matrix"
        raw_folder="raw_feature_bc_matrix"
    else:
        filter_folder="filtered_feature_bc_matrix_"+"scPathoQuant_"+args.aligner
        raw_folder="raw_feature_bc_matrix_"+"scPathoQuant_"+args.aligner

        if os.path.isdir(os.path.join(path10x, "outs", filter_folder)) is False:
           shutil.copytree(os.path.join(path10x, "outs", "filtered_feature_bc_matrix"), os.path.join(path10x, "outs", filter_folder), copy_function = shutil.copy)
        if os.path.isdir(os.path.join(path10x, "outs", raw_folder)) is False:
           shutil.copytree(os.path.join(path10x, "outs", "raw_feature_bc_matrix"), os.path.join(path10x,  "outs", raw_folder), copy_function = shutil.copy)

    #####-- add data to filter folder files --#####
    print ("STATUS: Integrating pathogen gene counts into 10x matrix and feature files in filter folders")
    # --Load in features and barcode matricies 
    feature_names_filtered = [row for row in csv.reader(gzip.open(os.path.join(path10x,"outs",filter_folder,"features.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    barcodes_filtered = [row[0] for row in csv.reader(gzip.open(os.path.join(path10x,"outs",filter_folder,"barcodes.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    feature_names_raw = [row for row in csv.reader(gzip.open(os.path.join(path10x,"outs",raw_folder,"features.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    barcodes_raw = [row[0] for row in csv.reader(gzip.open(os.path.join(path10x,"outs",raw_folder,"barcodes.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    
    counts4genes_filtered = list()
    counts4genes_raw = list()
    viruswithcounts=[]
    uniqgenesfull=[]
    counter=0
    for virus_name, dfumi in dfumidict.items():
        if dfumi is not False:
            viruswithcounts.append(virus_name)
            uniqgenes = set(dfumi['gene'].to_list())
            uniqgenesfull.extend(uniqgenes)
            dfumi = dfumi.reset_index()
            # --Put counts in to a list 
            for gene in uniqgenes:
                dftemp =dfumi.loc[dfumi.gene==gene, ]  
                counter=counter+1
                for index, row in dftemp.iterrows(): 
                    counts4genes_filtered.append([len(feature_names_filtered)+counter, barcodes_filtered.index(row["cell_barcode"])+1, row["umi"]])
                    counts4genes_raw.append([len(feature_names_raw)+counter, barcodes_raw.index(row["cell_barcode"])+1, row["umi"]])
        else: 
            print("STATUS: no gene read counts for virus "+virus_name+" so not integrating info into matrix files")
    # -- Fill list with gene names (for new features.tsv file)
    for gene in uniqgenesfull:
        feature_names_filtered.append((gene, gene, "Gene Expression"))
        feature_names_raw.append((gene, gene, "Gene Expression"))

    # -- unzip matrix file so that it can be opened 
    arg=['gunzip', "-f", os.path.join(path10x,"outs",filter_folder,"matrix.mtx.gz")]
    ef._run_subprocesses(arg, "STATUS: unzipping matrix file", "unzipping matrix file")

    # -- Open matrix file 
    f_in = open(os.path.join(path10x,"outs",filter_folder,"matrix.mtx"), 'r')
    lines = f_in.readlines()

    # -- Edit header in the matrix file to have the correct new number of genes
    larray =lines[2].split(" ")
    larray[0]=len(feature_names_filtered)
    larray[2] = len(lines)+len(counts4genes_filtered)-3
    lines[2]=" ".join(str(v) for v in larray)
    lines[2]=lines[2]+"\n"

    # -- Load edited data into new matrix file 
    f_in = open(os.path.join(path10x,"outs",filter_folder,"matrix.mtx"), 'w')
    for line in lines:
        f_in.write(line)
    for g in counts4genes_filtered:
        f_in.write(str(g[0])+" "+str(g[1])+" "+str(g[2])+"\n")
    f_in.close()

    # -- re zip matrix file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",filter_folder,"matrix.mtx")]
    ef._run_subprocesses(arg, "STATUS: zipping new matrix file", "zipping new matrix file")

    # -- load new information into new features file 
    f_in = open(os.path.join(path10x,"outs",filter_folder,"features.tsv"), 'w')
    for i in feature_names_filtered:
        f_in.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\n")
    f_in.close()

    # -- re zip features file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",filter_folder,"features.tsv")]
    ef._run_subprocesses(arg, "STATUS: zipping new features file in feature folder", "zipping new features file in feature folder")

    #####--add data to raw folder files--#####
    print ("STATUS: Integrating pathogen gene counts into 10x matrix and feature files in raw folder")

    # -- unzip matrix file so that it can be opened 
    arg=['gunzip', "-f", os.path.join(path10x,"outs",raw_folder,"matrix.mtx.gz")]
    ef._run_subprocesses(arg, "STATUS: unzipping matrix file", "unzipping matrix file")

    # -- Open matrix file 
    f_in = open(os.path.join(path10x,"outs",raw_folder,"matrix.mtx"), 'r')
    lines = f_in.readlines()

    # -- Edit header in the matrix file to have the correct new number of genes
    larray =lines[2].split(" ")
    larray[0]=len(feature_names_raw)
    larray[2] = len(lines)+len(counts4genes_raw)-3
    lines[2]=" ".join(str(v) for v in larray)
    lines[2]=lines[2]+"\n"

    # -- Load edited data into new matrix file 
    f_in = open(os.path.join(path10x,"outs",raw_folder,"matrix.mtx"), 'w')
    for line in lines:
        f_in.write(line)
    for g in counts4genes_raw:
        f_in.write(str(g[0])+" "+str(g[1])+" "+str(g[2])+"\n")
    f_in.close()

    # -- re zip matrix file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",raw_folder,"matrix.mtx")]
    ef._run_subprocesses(arg, "STATUS: zipping new matrix file", "zipping new matrix file")

    # -- load new information into new features file 
    f_in = open(os.path.join(path10x,"outs",raw_folder,"features.tsv"), 'w')
    for i in feature_names_filtered:
        f_in.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\n")
    f_in.close()

    # -- re zip features file 
    arg=['gzip', "-f", os.path.join(path10x,"outs",raw_folder,"features.tsv")]
    ef._run_subprocesses(arg, "STATUS: zipping new features file in raw folder", "zipping new features file in raw folder")    