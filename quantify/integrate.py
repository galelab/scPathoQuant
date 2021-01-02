__author__ = "Leanne Whitmore"
__email__ = "leanne382@gmail.com"
__description__ = "integrate data in to matrix"

from scipy.sparse import coo_matrix
import scipy.io
import os
import gzip
import csv
import numpy as np
import extra_functions as ef


def integrate_data_2_matrix(path10x, dfumi, gene_name):
    viralcopy = []
    print ("STATUS: Integrating data")
    mat = scipy.io.mmread(os.path.join(path10x,"outs","filtered_feature_bc_matrix","matrix.mtx.gz"))
    gene_names = [row for row in csv.reader(gzip.open(os.path.join(path10x,"outs","filtered_feature_bc_matrix","features.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
    barcodes = [row[0] for row in csv.reader(gzip.open(os.path.join(path10x,"outs","filtered_feature_bc_matrix","barcodes.tsv.gz"), "rt", encoding="utf8"), delimiter="\t")]
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
    scipy.io.mmwrite(os.path.join(path10x,"outs","filtered_feature_bc_matrix","matrix"), mat)
    arg=['gzip', os.path.join(path10x,"outs","filtered_feature_bc_matrix","matrix.mtx")]
    ef._run_subprocesses(arg, "STATUS: zipping new matrix file", "zipping new matrix file")
    f_in = open(os.path.join(path10x,"outs","filtered_feature_bc_matrix","features.tsv"), 'w')
    for i in gene_names:
        f_in.write(i[0]+"\t"+i[1]+"\t"+i[2]+"\n")
    f_in.close()
    arg=['gzip', os.path.join(path10x,"outs","filtered_feature_bc_matrix","features.tsv")]
    ef._run_subprocesses(arg, "STATUS: zipping new features file", "zipping new features file")    