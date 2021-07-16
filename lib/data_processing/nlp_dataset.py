from keras.backend import shape
import torch
from torch.utils.data.dataset import Dataset
import typing
from transformers import PreTrainedTokenizer, ElectraModel
from tqdm.notebook import tqdm
from ..utils import general as g
import numpy as np
# import time
# import sys
from ..config.config_class import page_config, legend_config
import os
import json
import pandas as pd
from ..genhelper.config_class import vcf_zarr_config as vzconfig

def masked_token(rand,tokens:typing.List[int],masked_value:int,masked_per=0.9,masked_indexs=None)->str:
    nb_elements = len(tokens)-2
    temp = np.asarray(tokens.copy())
    # get masked index by percent nb token data exclude start/end token
    if masked_indexs is None:
        masked_indexs = rand.permutation(nb_elements)[:int(nb_elements*masked_per)]+1
    temp[masked_indexs] = masked_value
    temp = temp.tolist()
    return temp

def token_default_dict(tokenizer:PreTrainedTokenizer):
    return {key:[] for key, val in tokenizer("").items()}

class GenNLPMaskedDataset(Dataset):
    """Genotype data to nlp load"""
    def __init__(self,document_paths:typing.List[str],tokenizer:PreTrainedTokenizer,seed=42,masked_per=0.9,masked_by_flag=False) -> None:
        super().__init__()
        self.rand = np.random
        self.rand.seed(seed=seed)
        self.tokenizer = tokenizer
        # inputs = tokenizer(masked_data,return_tensors='pt')
        # labels = tokenizer(gtruth_data,return_tensors='pt')['input_ids']
        self.labels = []
        self.input_ids = "input_ids"
        self.token_type_ids = "token_type_ids"
        self.attention_mask = "attention_mask"
        self.maskeds=token_default_dict(self.tokenizer)
        # reading and save data
        for i, dpath in enumerate(tqdm(document_paths,desc="preprocess data from document")):
            token_file = page_config.get_file_path_from_page(dpath,page_config.token)
            masked_file = page_config.get_file_path_from_page(dpath,page_config.masked,'.s'+str(seed))
            labels = []
            maskeds = token_default_dict(self.tokenizer)
            masked_indexs = None
            if masked_by_flag:
                variant_path = page_config.get_file_path_from_page(dpath,page_config.variant)
                variant_df = pd.read_csv(variant_path)
                masked_indexs = np.where(variant_df[vzconfig.flag].values == int(legend_config.observe))[0]
            if os.path.isfile(token_file):
                with g.reading(token_file) as tokenf:
                    # use masked variable as temp token variable
                    maskeds = json.load(tokenf)
                    labels = maskeds[self.input_ids].copy()
                if os.path.isfile(masked_file):
                    with g.reading(masked_file) as maskedf:
                        maskeds = json.load(maskedf)
                else:
                    maskeds[self.input_ids] = [masked_token(self.rand,label,self.tokenizer.mask_token_id,masked_per,masked_indexs) for label in labels]
                    page_config.token_masked_to_json(maskeds,masked_file)
            else:
                with g.reading(dpath) as temp_doc:
                    for line in  tqdm(temp_doc,desc="reading document ~ {:03d}".format(i),leave=False):
                        if line == '\n':
                            continue
                        # tokenizer data
                        data = self.tokenizer(line)
                        # copy true data to label data
                        label = data[self.input_ids].copy()
                        labels.append(label)
                        # masked input data from data
                        data[self.input_ids] = masked_token(self.rand,label,self.tokenizer.mask_token_id,masked_per,masked_indexs)
                        # append data to dataset
                        [maskeds[key].append(val) for key, val in data.items()]
                page_config.token_masked_to_json(maskeds,masked_file)
                page_config.token_masked_to_json({**maskeds,**{self.input_ids:labels}},token_file)
            # append data to this data
            self.labels.extend(labels)
            [self.maskeds[key].extend(val) for key, val in maskeds.items()]                        
        # conver labels and masked data to tensor
        self.labels = torch.as_tensor(self.labels)
        for key, val in self.maskeds.items():
            self.maskeds[key] = torch.as_tensor(val)
    
    def max_position_embeddings(self):
        return self.labels.shape[1]

    def __len__(self):
        return len(self.labels)
    
    def __getitem__(self, index):
        item = {}
        item[self.input_ids] = self.maskeds[self.input_ids][index]
        # item = {key: val[index] for key, val in self.maskeds.items()}
        item["labels"] = self.labels[index]
        return item
    

class GenElectraModel(ElectraModel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)