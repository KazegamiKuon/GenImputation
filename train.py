import sys
import typing
sys.path.insert(0,'./')
from lib.data_processing import GenNLPMaskedDataset
from transformers import ElectraForMaskedLM, ElectraTokenizer, ElectraConfig, TrainingArguments
import pandas as pd
import numpy as np
from lib.utils import general as g
from lib.config.config_class import page_config
from lib.model.overwriter import OTrainingArguments, OTrainer
from lib.utils.metrics import evalpred_to_word, r2_score_transformers
import json
import os
from IPython.display import clear_output
from typing import List

# config_path = '/client/user1/cuongdev/GenImputation/data/train/electra_G1K_22_hs37d5/config.json'
# config_path = '/client/user1/cuongdev/GenImputation/data/train/electra_G1K_22_hs37d5/config_demo_r9.json'
config_path = '/client/user1/cuongdev/GenImputation/data/train/electra_G1K_22_hs37d5/config_demo_r9_same_length.json'
# config_path = '/client/user1/cuongdev/GenImputation/data/train/electra_G1K_22_hs37d5/config_demo_r3.json'

config = None
with g.reading(config_path) as cf:
    config = json.load(cf)
assert config is not None, "config can't none"

regions = config['regions']
batchs = config['batchs']
detail = config['model_name']
mode = config['mode']
masked_mode = None
if 'masked_mode' in config:
    masked_mode = config['masked_mode']
masked_flag_random = masked_mode == 'random'

train_region_paths = page_config.get_file_paths(config[page_config.file_train_prefix],page_config.page,regions,batchs)
test_region_paths = page_config.get_file_paths(config[page_config.file_test_prefix],page_config.page,regions,[0])
vocab_file = config[page_config.vocab_file]
save_dir = config[page_config.save_dir]

training_args = OTrainingArguments(**config[page_config.train_args])
output_dir = training_args.output_dir
logging_dir = training_args.logging_dir
modeling_args = ElectraConfig(**config[page_config.model_args])
tokenizer = ElectraTokenizer(vocab_file=vocab_file)
seed = training_args.seed

for i, region in enumerate(regions):
    clear_output(wait=True)
    print('Region {} trainning...'.format(region))
    print('Prevert region {} trainning...'.format(region-1))
    save_path = save_dir.format(region)+detail
    prevert_path = save_dir.format(region-1)+detail
    training_args.output_dir = output_dir.format(region)+detail
    training_args.logging_dir = logging_dir.format(region)+detail
    ## Train and eval data
    train_batch_paths = train_region_paths[i]
    train_dataset = GenNLPMaskedDataset(train_batch_paths, tokenizer, seed=seed, masked_by_flag=True, only_input=True,masked_flag_random = masked_flag_random)
    # eval_dataset = GenNLPMaskedDataset(train_batch_paths[-1:], tokenizer, seed=seed, masked_by_flag=True, only_input=True)
    ## test data
    test_batch_paths = test_region_paths[i]
    test_dataset = GenNLPMaskedDataset(test_batch_paths,tokenizer,seed=seed,masked_by_flag=True,only_input=True)
    ## model
    modeling_args.vocab_size = tokenizer.vocab_size
    if mode == 'pretrain':
        modeling_args.max_position_embeddings = 1300
    else:
        modeling_args.max_position_embeddings = train_dataset.max_position_embeddings()
    electra_model = ElectraForMaskedLM(modeling_args)
    if os.path.isdir(prevert_path):
        electra_model = ElectraForMaskedLM.from_pretrained(prevert_path)
    trainer = OTrainer(
        model = electra_model,
        args=training_args,
        train_dataset = train_dataset,
        eval_dataset = test_dataset,
        compute_metrics = r2_score_transformers,
    )
    trainer.train()
    trainer.save_model(save_path)
    output_test = trainer.predict(test_dataset)
    metrics = output_test.metrics
    test_result_path = os.path.join(save_path,'test_result.json')
    with g.writing(test_result_path) as trf:
        json.dump(metrics,trf)