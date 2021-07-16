import sys
import pandas as pd
from lib.utils import general as g
import os
import unittest
import numpy as np
from lib.data_processing import process_input as pi
from lib.data_processing import GenNLPMaskedDataset
import sys
from transformers import ElectraTokenizer, ElectraConfig, Trainer, TrainingArguments
from transformers import EncoderDecoderModel, EncoderDecoderConfig
from lib.model.electra import ElectraEmbeddings, ElectraForMaskedLM
import torch
import numpy as np

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

        self.vcf_inter_file = '/client/user1/cuongdev/GenImputation/data/interim/vn_isec_1k_data_private/VN_20_hg38.vcf.gz'
        self.vcf_file_g1k_train = '/client/user1/cuongdev/GenImputation/data/train/G1K_chr20_biallelic_train.vcf.gz'
        self.nlp_train_prefix = '/client/user1/cuongdev/GenImputation/data/train/electra/temp/G1K_VN_chr20_biallelic_train'

        self.vocab_file = '/client/user1/cuongdev/GenImputation/data/train/electra/data_dir/vocab.txt'
        self.train_paths = [
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir/G1K_VN_chr20_biallelic_train.r0000.b0000.page.gz',
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir/G1K_VN_chr20_biallelic_train.r0000.b0001.page.gz',
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir/G1K_VN_chr20_biallelic_train.r0000.b0002.page.gz'
        ]
        self.eval_paths = [
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir/G1K_VN_chr20_biallelic_train.r0000.b0003.page.gz'
        ]

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
    
    # def test_process_data(self):
    #     pi.process_data_to_legend(self.vcf_file,self.manifest_file,self.hg_fasta_file,self.chroms,self.my_output_prefix)
    
    # def test_legend_to_region(self):
    #     # pi.legend_to_region(self.legend_file,self.hap_file,100,0.1,True,output_folder=self.output_folder)
    #     pi.process_legend_to_region(self.legend_file,self.hap_file,100,0.1,True)
    
    # def test_region_config(self):
    #     pi.process_region_config(self.region_folder,self.default_config)
    # # def test_genotyping_vcf(self):
    # #     pi.genotyping_vcf(self.vcf_file,self.manifest_file,self.hg_fasta_file,self.my_output_prefix,self.chroms)
    
    # def test_ssh_config(self):
    #     manifest_file = "/client/user1/data_imputation_vn/APMRA96_chr20_for_imputation_chr20_AF.vcf.gz"
    #     vcf_file = "/client/user1/data_imputation_vn/ref1014_chr20_for_imputation_chr20_AF.vcf.gz"
    #     fasta_file = "/client/user1/cuongdev/GenImputation/data/raw/hg38.fa.gz"
    #     chroms = ['chr20']
    #     ouput_prefix = '/client/user1/cuongdev/GenImputation/data/interim/ref1014_chr20_hg38'
    #     pi.process_data_to_legend(vcf_file,manifest_file,fasta_file,chroms,ouput_prefix)
    
    def test_page(self):
        pi.process_vcf_to_page_nlp(self.vcf_file_g1k_train,[self.vcf_inter_file],0,self.nlp_train_prefix,136,9)
    
    def test_nlp_dataset(self):
        train_dataset = GenNLPMaskedDataset(self.train_paths,self.vocab_file,seed=42)
        pass

    def test_electra(self):
        vocab_file = '/client/user1/cuongdev/GenImputation/data/train/electra/data_dir/vocab.txt'
        train_paths = [
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir/G1K_VN_chr20_biallelic_train.r0000.b0000.page.gz',
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir/G1K_VN_chr20_biallelic_train.r0000.b0001.page.gz',
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir/G1K_VN_chr20_biallelic_train.r0000.b0002.page.gz'
        ]
        eval_paths = [
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir/G1K_VN_chr20_biallelic_train.r0000.b0003.page.gz'
        ]
        output_dir = '/client/user1/cuongdev/GenImputation/data/train/electra/checkpoints'
        do_train = True
        do_eval = True
        do_predict = True
        evaluation_stategy = 'epoch'
        learning_rate=5e-5
        weight_decay = 0
        num_train_epochs = 20
        lr_scheduler_type = "linear"
        save_strategy = "epoch"
        no_cuda = True
        seed = 42
        run_name = "cuongdev_electra"
        per_device_train_batch_size=1
        per_device_eval_batch_size=1
        embedding_size = 16
        hidden_size = 64
        position_embedding_type =  'relative_key'
        num_hidden_layers = 1
        # data
        train_dataset = GenNLPMaskedDataset(train_paths,vocab_file,seed=seed)
        eval_dataset = GenNLPMaskedDataset(eval_paths,vocab_file,seed=seed)
        # config
        max_position_embeddings = train_dataset.max_position_embeddings()
        config = ElectraConfig(
            vocab_size=8,
            max_position_embeddings=max_position_embeddings,
            embedding_size = embedding_size ,
            hidden_size = hidden_size,
            position_embedding_type=position_embedding_type,
            num_hidden_layers = num_hidden_layers
        )
        # config.return_dict = True
        electra_model = ElectraForMaskedLM(config)
        ### data
        # data = train_dataset.__getitem__(0)
        # data = {key: torch.unsqueeze(val,dim=0) for key, val in data.items()}
        # temp_inputs = electra_model(**data)
        ###
        config = electra_model.config
        train_args = TrainingArguments(
            output_dir=output_dir,
            do_train=do_train,
            do_eval=do_eval,
            do_predict=do_predict,
            evaluation_strategy=evaluation_stategy,
            learning_rate=learning_rate,
            weight_decay=weight_decay,
            num_train_epochs=num_train_epochs,
            lr_scheduler_type=lr_scheduler_type,
            save_strategy = save_strategy,
            # save_steps=save_strategy,
            no_cuda=no_cuda,
            seed=seed,
            run_name=run_name,
            per_device_train_batch_size=per_device_train_batch_size,
        )
        trainer = Trainer(
            model = electra_model,
            args=train_args,
            train_dataset=train_dataset,
            eval_dataset = eval_dataset,
        )
        trainer.train()
        pass

    def test_how_model_run(self):
        vocab_file = '/client/user1/cuongdev/GenImputation/data/train/electra/data_dir/vocab.txt'
        train_paths = [
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir_2048/G1K_VN_chr20_biallelic_train.r0000.b0000.page.gz',
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir_2048/G1K_VN_chr20_biallelic_train.r0000.b0001.page.gz',
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir_2048/G1K_VN_chr20_biallelic_train.r0000.b0002.page.gz',
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir_2048/G1K_VN_chr20_biallelic_train.r0000.b0003.page.gz',
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir_2048/G1K_VN_chr20_biallelic_train.r0000.b0004.page.gz',
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir_2048/G1K_VN_chr20_biallelic_train.r0000.b0005.page.gz',
        ]
        eval_paths = [
            '/client/user1/cuongdev/GenImputation/data/train/electra/corpus_dir_2048/G1K_VN_chr20_biallelic_train.r0000.b0007.page.gz'
        ]
        tokenizer = ElectraTokenizer(vocab_file=vocab_file)
        eval_dataset = GenNLPMaskedDataset(eval_paths,tokenizer,seed=42,masked_per=0.15)
        model = ElectraForMaskedLM.from_pretrained('/client/user1/cuongdev/GenImputation/data/train/electra/checkpoints/small_2048/checkpoint-8568')
        output = model(**eval_dataset.__getitem__([0,1,2,3]))
        # output = model(**{key: torch.unsqueeze(val,0) for key, val in eval_dataset.__getitem__([0,1,2,3]).items()})
        pass
    
    def test_map_marker(self):
        manifest = '/client/user1/cuongdev/GenImputation/data/interim/vn_data_private/manifest_vn.vcf.gz'
        variant_paths = ['/client/user1/cuongdev/GenImputation/data/test/electra/corpus_dir_2048/G1K_VN_chr20_biallelic_test.r0000.b0000.variant.gz',
        '/client/user1/cuongdev/GenImputation/data/test/electra/corpus_dir_2048/G1K_VN_chr20_biallelic_test.r0001.b0000.variant.gz']
        hg_genome = '/client/user1/cuongdev/GenImputation/data/raw/hg38.fa.gz'
        pi.process_map_manifest_to_variant_nlp(variant_paths,manifest,hg_genome)
        pass

if __name__ == '__main__':
    unittest.main()