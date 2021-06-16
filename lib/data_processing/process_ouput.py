from lib.genhelper import vcf_helper as vhelper
import pandas as pd
import os
import numpy as np
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt

def plot_r2_by_maf(mafs,y_true,y_preds,bins=[0,0.001,0.005,0.01,0.05,0.2,0.5],labels=None):
    '''
    input:
        mafs: list maf data of each position
        y_true: groud truth data
        y_preds: list predict data or dict which value is predict data
        bins,label: see pandas.cut params to expland more
    '''

    if labels is None:
        labels = bins[1:]
    # convert ypreds to dict
    if type(y_preds) is not dict:
        keys = np.arange(len(y_preds))
        y_preds = dict(zip(keys,y_preds))
    
    nb_maf = len(mafs)
    maf_col_name = 'MAF'
    # stored index to map with x and y
    index_col_name = 'INDEX'
    bin_col_name = 'BIN'
    df_maf = pd.DataFrame({maf_col_name:mafs,index_col_name:np.arange(nb_maf)})
    df_maf[bin_col_name] = pd.cut(df_maf[maf_col_name],bins=bins,labels=labels)
    # df_maf[bin_col_name], labels = pd.qcut(df_maf[maf_col_name],q=nb_value_per_bin,labels=False,retbins=True)
    # labels = labels[1:]
    r2_dict = {}
    # nb_label = len(labels) # number label
    for label in labels:
        # get value from columns INDEX
        indexs = df_maf[df_maf[bin_col_name] == label][index_col_name].values.flatten()
        for key, y_pred in y_preds.items():
            temp_score = r2_score(y_true[indexs],y_pred[indexs])
            if key in r2_dict:                
                r2_dict[key].append(temp_score)
            else:
                r2_dict[key] = [temp_score]
    
    xaxis = np.arange(len(labels))  # x axis

    # draw plot and change r2_dict value to np array
    for key in r2_dict:
        r2_dict[key] = np.array(r2_dict[key])
        plt.plot(np.arange(len(labels)),r2_dict[key],label=''.join(str(key)))
    plt.xticks(ticks=xaxis,labels=labels)
    plt.legend()
    plt.show()
    return labels, r2_dict

