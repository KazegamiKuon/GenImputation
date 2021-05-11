import sys
sys.path.insert(1,'../../')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import genhelper.vcf_helper as vhelper
from sklearn.metrics import r2_score
import zarr

def af_to_maf(x):
    return x if x < 0.5 else 1 - x

def r2_score_by_bin(labels,data,x,y):
    r2_scores = []
    nb_var = []
    for label in labels:
        temp = df[df[bin_name] == label]
        xindexs = temp[xindex_col].values
        yindexs = temp[yindex_col].values
        r2_ = r2_score(ygt[yindexs],xgt[xindexs])
        r2_scores.append(r2_)
        nb_var.append(temp.shape[0])
    r2_scores = np.array(r2_scores)
    nb_var = np.array(nb_var)
    return r2_scores, nb_var