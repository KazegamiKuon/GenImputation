from lib.genhelper import vcf_helper as vhelper
import pandas as pd
import os
import numpy as np
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import zarr

def plot_r2_by_maf(mafs,y_true,y_preds,bins=[-1,0.001,0.005,0.01,0.05,0.2,0.5],labels=None,draw=True):
    '''
    input:
        mafs: list maf data of each position
        y_true: groud truth data
        y_preds: list predict data or dict which value is predict data
        bins,label: see pandas.cut params to expland more
        draw: draw plot or not
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
        if draw:
            plt.plot(np.arange(len(labels)),r2_dict[key],label=''.join(str(key)))
    if draw:
        plt.xticks(ticks=xaxis,labels=labels)
        plt.legend()
        plt.show()
    return labels, r2_dict
import types
def get_r2_score_minimac_result(true_callset:zarr.Group,pred_callset:zarr.Group,source_callset:zarr.Group,sample_func:types.FunctionType,dbname:str,**kwargs):
    '''
    input:
        true_callset, pred_callset, source_callset: zarr.Group type data.
            Run by vcf_helper.vcf_to_zarr to get it zarr path and use zarr.open_group to get that call set
            True_callset from ground truth data and pred_callset from minimac predict data
            Source_callset from minimac train data
        dbname: name of dataset
        sample_func: (type:["true","pred"],sample_name)=> return str
            use for get name of sample to mapping if data true and pred have different sample
    '''
    gt=dbname+'_gt'
    ds=dbname+'_ds'
    true_samples = [sample_func('true',sample) for sample in true_callset.samples[:]]
    pred_samples = [sample_func('pred',sample) for sample in pred_callset.samples[:]]
    pred_true_mask_samples = [sample in true_samples for sample in pred_samples]
    true_pred_mask_samples = [sample in pred_samples for sample in true_samples]
    intersection_variant_id = vhelper.get_dataframe_variant_id([true_callset.variants,pred_callset.variants])
    true_indexs=intersection_variant_id['index_0'].values
    pred_indexs=intersection_variant_id['index_1'].values
    afs = source_callset.variants.AF[:][:,0]
    mafs = [af if af <= 0.5 else 1-af for af in afs]
    
    y_true = true_callset.calldata.GT[:][true_indexs]
    y_true = y_true[:,true_pred_mask_samples,:]
    y_true = y_true.reshape((y_true.shape[0],y_true.shape[1]*2))
    y_pred = pred_callset.calldata.GT[:][pred_indexs]
    y_pred = y_pred[:,pred_true_mask_samples,:]
    y_pred = y_pred.reshape((y_pred.shape[0],y_pred.shape[1]*2))
    labels, gt_r2_dict = plot_r2_by_maf(mafs=mafs,y_true=y_true,y_preds={gt:y_pred},**kwargs)
    y_true = true_callset.calldata.GT[:][true_indexs]
    y_true = y_true[:,true_pred_mask_samples,:]
    y_true = np.sum(y_true,axis=2)
    y_pred = pred_callset.calldata.DS[:][pred_indexs]
    y_pred = y_pred[:,pred_true_mask_samples]
    labels, ds_r2_dict = plot_r2_by_maf(mafs=mafs,y_true=y_true,y_preds={ds:y_pred},**kwargs)
    r2_dict = dict({
        gt:gt_r2_dict[gt],
        ds:ds_r2_dict[ds]
    })
    return labels, r2_dict

def plot_label_and_dict(labels:list,data_dict:dict,title=""):
    for key in data_dict:
        data_dict[key] = np.array(data_dict[key])
        label = key if type(key) is str else str(key)
        plt.plot(np.arange(len(labels)),data_dict[key],label=label)
    xaxis = np.arange(len(labels))
    plt.xticks(ticks=xaxis,labels=labels)
    plt.title(title)
    plt.legend()
    plt.show()