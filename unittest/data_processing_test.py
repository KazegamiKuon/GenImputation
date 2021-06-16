import sys
import pandas as pd
from lib.utils import general as g
import os
import unittest
import numpy as np
from lib.data_processing import process_input as pi

class DataProcessingTest(unittest.TestCase):
    def __init__(self, methodName: str) -> None:
        super().__init__(methodName=methodName)
        self.vcf_file = '/home/cuong/VBDI/HungProject/GenImputation/data/raw/G1K_chr20_hg38.vcf.gz'
        self.manifest_file='/home/cuong/VBDI/HungProject/GenImputation/data/raw/infiniumomni2-5-8-v1-3-a2.csv.gz'
        self.hg_fasta_file='/home/cuong/VBDI/HungProject/GenImputation/data/raw/hg38.fa.gz'
        self.my_output_prefix = '/home/cuong/VBDI/HungProject/GenImputation/data/interim/G1K_chr22_hg38_manifest'
        self.chroms=['20']
        self.test_sample_list_file='/home/cuong/VBDI/HungProject/GenImputation/data/external/test_100_samples.txt'
        self.legend_file = '/home/cuong/VBDI/HungProject/GenImputation/data/interim/G1K_chr20_hg38_true.legend.gz'
        self.hap_file = '/home/cuong/VBDI/HungProject/GenImputation/data/interim/G1K_chr20_hg38_true.hap.gz'
        self.region_folder = '/home/cuong/VBDI/HungProject/GenImputation/data/interim/region'
        self.output_folder = '/home/cuong/VBDI/HungProject/GenImputation/data/train'
        self.default_config = '/home/cuong/VBDI/HungProject/GenImputation/data/external/region_default_config.json'

    # def test_plot_r2_by_maf(self):
    #     #params
    #     concat_datafame_params ={
    #         'header':None,
    #         'sep':' '
    #     }
    #     group_file_key = lambda x: x.split('_')[:-2]
    #     group_key_index = 0
    #     key_cols=[0,1,2,3,4]
    #     key_col_name = 'KEY'
    #     # get group data
    #     # predict
    #     group_list = g.get_group_file('/home/cuong/VBDI/HungProject/GenImputation/data/test/higher_model_result',key = group_file_key)
    #     group_predict_data = {}
    #     for k, paths in group_list[0]:
    #         key = ''.join(k)
    #         temp = list(paths)
    #         group_predict_data[key] = g.concat_dataframe(temp,**concat_datafame_params)
    #     # paper predict
    #     group_list = g.get_group_file('/home/cuong/VBDI/HungProject/GenImputation/data/test/paper_hybrid_result',key = group_file_key)
    #     group_paper_data = {}
    #     for k, paths in group_list[0]:
    #         key = ''.join(k)
    #         group_paper_data[key] = g.concat_dataframe(list(paths),**concat_datafame_params)
    #     # groud truth
    #     group_list = g.get_group_file('/home/cuong/VBDI/HungProject/GenImputation/data/test/gt_chr22_1_5',key = group_file_key)
    #     group_gtruth_data = {}
    #     for k, paths in group_list[0]:
    #         key = ''.join(k)
    #         group_gtruth_data[key] = g.concat_dataframe(list(paths),**concat_datafame_params)

    #     # get data predict
    #     keys = list(group_predict_data.keys())
    #     predict_data = group_predict_data[keys[group_key_index]]
    #     predict_data_keys = predict_data[key_cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)

    #     # get data predict
    #     keys = list(group_paper_data.keys())
    #     paper_data = group_paper_data[keys[group_key_index]]
    #     predict_paper_keys = paper_data[key_cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)

    #     # get data gt
    #     keys = list(group_gtruth_data.keys())
    #     gtruth_data = group_gtruth_data[keys[group_key_index]]
    #     gtruth_data_keys = gtruth_data[key_cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)

    #     assert np.all(predict_data_keys == gtruth_data_keys), 'not same key'
    #     y_paper_pred = paper_data.values[:,6:]
    #     y_me_pred = predict_data.values[:,5:]
    #     mafs = gtruth_data.values[:,5]
    #     y_true = gtruth_data.values[:,6:]
    #     labels, r2_dict = po.plot_r2_by_maf(mafs,y_true,{'paper':y_paper_pred,'own':y_me_pred})
    #     print('')
    
    def test_process_data(self):
        pi.process_data_to_legend(self.vcf_file,self.manifest_file,self.hg_fasta_file,self.chroms,self.my_output_prefix)
    
    def test_legend_to_region(self):
        # pi.legend_to_region(self.legend_file,self.hap_file,100,0.1,True,output_folder=self.output_folder)
        pi.legend_to_region(self.legend_file,self.hap_file,100,0.1,True)
    
    def test_region_config(self):
        pi.make_region_config(self.region_folder,self.default_config)
    # def test_genotyping_vcf(self):
    #     pi.genotyping_vcf(self.vcf_file,self.manifest_file,self.hg_fasta_file,self.my_output_prefix,self.chroms)

if __name__ == '__main__':
    unittest.main()