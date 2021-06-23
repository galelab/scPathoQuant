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
    viralcopy = []
    if args.overwrite_feature_matrix:
        filter_folder="filtered_feature_bc_matrix"
    else:
        filter_folder="filtered_feature_bc_matrix_"+args.aligner+"_"+args.alignment_type
        if os.path.isdir(os.path.join(path10x, "outs", filter_folder)) is False:
            shutil.copytree(os.path.join(path10x, "outs", "filtered_feature_bc_matrix"), os.path.join(path10x, "outs", filter_folder), copy_function = shutil.copy)
        else:
            print("WARNING: folder exists so removing and recopying original")
            shutil.rmtree(os.path.join(path10x, "outs", filter_folder))
            shutil.copytree(os.path.join(path10x, "outs", "filtered_feature_bc_matrix"), os.path.join(path10x, "outs", filter_folder), copy_function = shutil.copy)

    print ("STATUS: Integrating viral copy counts into 10x matrix and feature files")
    mat = scipy.io.mmread(os.path.join(path10x,"outs", filter_folder, "matrix.mtx.gz"))
    gene_names = [row for row in csv.reader(gzip.open(os.path.join(path10x,"outs", filter_folder, "features.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    barcodes = [row[0] for row in csv.reader(gzip.open(os.path.join(path10x,"outs", filter_folder, "barcodes.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    for i in barcodes:
        try:
            count = dfumi.loc[i,'umi']
            viralcopy.append(count)
        except KeyError:
            viralcopy.append(0)
    npviralcopy = np.array(viralcopy)
    mat = mat.todense()
    mat = np.vstack([mat,npviralcopy])
    mat = coo_matrix(mat)
    gene_names.append((gene_name, gene_name, "Gene Expression"))
    scipy.io.mmwrite(os.path.join(path10x,"outs",filter_folder,"matrix"), mat)
    arg=['gzip', "-f", os.path.join(path10x,"outs",filter_folder,"matrix.mtx")]
    ef._run_subprocesses(arg, "STATUS: zipping new matrix file", "zipping new matrix file")
    f_in = open(os.path.join(path10x,"outs",filter_folder,"features.tsv"), 'w')
    for i in gene_names:
        f_in.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\n")
    f_in.close()
    arg=['gzip', "-f", os.path.join(path10x,"outs",filter_folder,"features.tsv")]
    ef._run_subprocesses(arg, "STATUS: zipping new features file", "zipping new features file")


def integrate_viralgenes_data_2_matrix(args, dfumi):
    path10x = args.path10x
    viralcopy = []
    if args.overwrite_feature_matrix:
        filter_folder="filtered_feature_bc_matrix"
    else:
        filter_folder="filtered_feature_bc_matrix_"+args.aligner+"_"+args.alignment_type
        if os.path.isdir(os.path.join(path10x, "outs", filter_folder)) is False:
           shutil.copytree(os.path.join(path10x, "outs", "filtered_feature_bc_matrix"), os.path.join(path10x,  "outs", filter_folder), copy_function = shutil.copy)
    print ("STATUS: Integrating viral gene counts into 10x matrix and feature files")
    mat = scipy.io.mmread(os.path.join(path10x,"outs",filter_folder,"matrix.mtx.gz"))
    gene_names = [row for row in csv.reader(gzip.open(os.path.join(path10x,"outs",filter_folder,"features.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    barcodes = [row[0] for row in csv.reader(gzip.open(os.path.join(path10x,"outs",filter_folder,"barcodes.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    uniqgenes = set(dfumi['gene'].to_list())
    counts4genes = {}
    dfumi = dfumi.reset_index()
    for gene in uniqgenes:
        for i in barcodes:
            try:
                count = dfumi.loc[(dfumi['cell_barcode']==i) & (dfumi['gene'] == gene), 'umi'].iloc[0]
                counts4genes.setdefault(gene, []).append(count)
            except (KeyError, IndexError) as e:
                counts4genes.setdefault(gene, []).append(0)
    npgenes = np.array(list(counts4genes.values()))
    mat = mat.todense()
    mat = np.vstack([mat,npgenes])
    mat = coo_matrix(mat)
    for gene in list(counts4genes.keys()):
        gene_names.append((gene, gene, "Gene Expression"))
    scipy.io.mmwrite(os.path.join(path10x,"outs",filter_folder,"matrix"), mat)
    arg=['gzip', "-f", os.path.join(path10x,"outs",filter_folder,"matrix.mtx")]
    ef._run_subprocesses(arg, "STATUS: zipping new matrix file", "zipping new matrix file")
    f_in = open(os.path.join(path10x,"outs",filter_folder,"features.tsv"), 'w')
    for i in gene_names:
        f_in.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\n")
    f_in.close()
    arg=['gzip', "-f", os.path.join(path10x,"outs",filter_folder,"features.tsv")]
    ef._run_subprocesses(arg, "STATUS: zipping new features file", "zipping new features file")    