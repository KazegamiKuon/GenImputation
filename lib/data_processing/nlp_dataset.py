import torch
from torch.utils.data.dataset import Dataset
import typing
from transformers import ElectraTokenizer
from tqdm.notebook import tqdm
from ..utils import general as g
import numpy as np
# import time
# import sys
from ..config.config_class import page_config
import os
import json

def masked_token(rand,tokens:typing.List[int],masked_value:int,masked_per=0.9)->str:
    nb_elements = len(tokens)-2
    temp = np.asarray(tokens.copy())
    # get masked index by percent nb token data exclude start/end token
    masked_indexs = rand.permutation(nb_elements)[:int(nb_elements*masked_per)]+1
    temp[masked_indexs] = masked_value
    return list(temp)

def token_default_dict(tokenizer:ElectraTokenizer):
    return {key:[] for key, val in tokenizer("").items()}

class GenNLPMaskedDataset(Dataset):
    """Genotype data to nlp load"""
    def __init__(self,document_paths:typing.List[str],vocab_file:str,nb_clone:int = 1,seed=42,masked_per=0.9) -> None:
        super().__init__()
        self.rand = np.random
        self.rand.seed(seed=seed)
        self.tokenizer = ElectraTokenizer(vocab_file=vocab_file,tokenize_chinese_chars=False)
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
            masked_file = page_config.get_file_path_from_page(dpath,page_config.masked,str(seed))
            labels = []
            maskeds = token_default_dict(self.tokenizer)

            if os.path.isfile(token_file):
                with g.reading(token_file) as tokenf:
                    # use masked variable as temp token variable
                    maskeds = json.load(tokenf)
                    labels = maskeds[self.input_ids].copy()
                if os.path.isfile(masked_file):
                    with g.reading(masked_file) as maskedf:
                        maskeds = json.load(maskedf)
                        continue
                else:
                    maskeds[self.input_ids] = [masked_token(self.rand,label,self.tokenizer.mask_token_id,masked_per) for label in labels]
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
                        data[self.input_ids] = masked_token(self.rand,label,self.tokenizer.mask_token_id,masked_per)
                        # append data to dataset
                        [maskeds[key].append(val) for key, val in data.items()]
                page_config.token_masked_to_json(maskeds,masked_file)
                page_config.token_masked_to_json({**maskeds,**{self.input_ids:labels}},token_file)
            # append data to this data
            self.labels.extend(labels)
            [self.maskeds[key].extend(val) for key, val in maskeds]                        
        # conver labels and masked data to tensor
        self.labels = torch.as_tensor(self.labels)
        for key, val in self.maskeds.items():
            self.maskeds[key] = torch.as_tensor(val)
    
    def __len__(self):
        return len(self.labels)
    
    def __getitem__(self, index):
        if torch.is_tensor(index):
            index = index.toList()
        item = {key: val[index] for key, val in self.maskeds.items()}
        item["labels"] = self.labels[index]
        return item
